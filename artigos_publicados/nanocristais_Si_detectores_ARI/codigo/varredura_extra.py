#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
varredura_extra.py -- complemento de varredura_gap.py.

(e1) Estruturas "min_viz = 2 + preenchimento" (sitios externos com >= 2 vizinhos
     internos sao incorporados; elimina pares de H apontando ao mesmo sitio vago):
     D = 1,2 nm -> Si39H40 (STO-3G, HF e PBE, integrais exatas).
     [D = 0,8 e 1,0 nm sao identicas a min_viz = 2; D = 0,9 nm e identica a D = 1,0 nm/min_viz = 2.]
(e2) Erro do density fitting (DF) em STO-3G: Si36H40 (D = 1,2 nm, min_viz = 2) com DF.
(e3) D = 1,5 nm com DF: Si72H64 (min_viz = 2), Si84H64 (min_viz = 2 + preenchimento)
     e, por ultimo, Si84H88 (min_viz <= 1).

Reutiliza as funcoes de varredura_gap.py e grava no MESMO ../dados/gaps.json
(chaves com 'mv2p' para o criterio com preenchimento; campo 'criterio' em todos os registros).
Uso: python3 varredura_extra.py [--etapas e1 e2 e3]
"""
import os, sys, time, json, argparse
import numpy as np
from scipy.spatial import cKDTree
from pyscf import gto, scf, dft
import varredura_gap as vg

def vizinhos_rede(rede):
    arv = cKDTree(rede)
    pares = arv.query_pairs(vg.CORTE, output_type="ndarray")
    viz = [[] for _ in range(len(rede))]
    for i, j in pares:
        viz[i].append(j); viz[j].append(i)
    return [np.array(v, dtype=int) for v in viz]

def construir_preencher(diametro_nm, n_celulas=8, min_viz=2):
    """Mesmo algoritmo de vg.construir + criterio de preenchimento (copia fiel de
    analise_geometria.construir_info(preencher=True))."""
    raio = diametro_nm * 10 / 2
    rede = vg.rede_diamante(n_celulas)
    rede = rede - rede.mean(axis=0)
    viz = vizinhos_rede(rede)
    dentro = np.linalg.norm(rede, axis=1) <= raio
    it = 0
    while True:
        it += 1
        nviz = np.array([np.sum(dentro[viz[i]]) if dentro[i] else -1 for i in range(len(rede))])
        ruins = np.where(dentro & (nviz < min_viz))[0]
        dentro[ruins] = False
        nviz_ext = np.array([np.sum(dentro[viz[i]]) if not dentro[i] else -1 for i in range(len(rede))])
        novos = np.where((~dentro) & (nviz_ext >= 2))[0]
        dentro[novos] = True
        if len(ruins) == 0 and len(novos) == 0:
            break
        if it > 500:
            raise RuntimeError("nao convergiu")
    si = rede[dentro]
    H = []
    for i in np.where(dentro)[0]:
        p = rede[i]
        for j in viz[i]:
            if not dentro[j]:
                u = (rede[j] - p) / np.linalg.norm(rede[j] - p)
                H.append(p + vg.D_SIH * u)
    return si, np.array(H).reshape(-1, 3)

def calcular_geom(db, etapa, D, tag, criterio, si, h, base, metodo, df=False):
    """Igual a vg.calcular, mas recebe a geometria pronta e uma etiqueta de criterio."""
    fm = vg.formula(si, h)
    chave = "D%.1f|%s|%s|%s|%s|%s" % (D, tag, fm, base, metodo, 'DF' if df else 'exato')
    if chave in db['calculos'] and db['calculos'][chave].get('converged'):
        print("  [ja feito] " + chave, flush=True); return db['calculos'][chave]
    ekey = "D%.1f|%s" % (D, tag)
    if ekey not in db['estruturas']:
        dg = vg.diagnostico(si, h); dg['formula'] = fm; dg['D_nominal_nm'] = D
        dg['min_viz'] = 2; dg['criterio'] = criterio
        dg['D_efetivo_nm'] = vg.diametro_efetivo_nm(len(si))
        db['estruturas'][ekey] = dg
    mol = gto.M(atom=vg.para_xyz(si, h), basis=base, unit='Angstrom', verbose=0,
                max_memory=vg.MAX_MEMORY, symmetry=False)
    Def = vg.diametro_efetivo_nm(len(si))
    reg = dict(etapa=etapa, D_nominal_nm=D, min_viz=(2 if tag == 'mv2' else (1 if tag == 'mv1' else '2p')),
               criterio=criterio, formula=fm, n_si=len(si), n_h=len(h), D_efetivo_nm=Def,
               base=base, metodo=metodo, nao=int(mol.nao), nelectron=int(mol.nelectron),
               threads=vg.N_THREADS, gap_brus_Dnominal_eV=float(vg.gap_brus(D)),
               gap_brus_Defetivo_eV=float(vg.gap_brus(Def)), gap_delerue_Defetivo_eV=float(vg.gap_delerue(Def)))
    if mol.nelectron % 2:
        reg.update(converged=False, obs='numero impar de eletrons: camada aberta, pulado')
        db['calculos'][chave] = reg; vg.salvar(db); return reg
    print("  -> %s  nao=%d  nel=%d" % (chave, mol.nao, mol.nelectron), flush=True)
    obs = []
    for ls in (0.0, 0.3):
        mf = scf.RHF(mol) if metodo.upper() == 'HF' else dft.RKS(mol)
        if metodo.upper() != 'HF':
            mf.xc = metodo.lower()
        if df:
            mf = mf.density_fit()
        mf.conv_tol = 1e-8; mf.max_cycle = 150; mf.level_shift = ls; mf.verbose = 0
        t0 = time.time()
        try:
            e = mf.kernel()
        except Exception as exc:
            reg.update(converged=False, obs='excecao: %r' % exc); db['calculos'][chave] = reg; vg.salvar(db)
            print("  !! excecao", exc, flush=True); return reg
        dt = time.time() - t0
        if mf.converged:
            break
        obs.append('nao convergiu com level_shift=%.1f' % ls)
    modo = ('DF (auxbase %s)' % getattr(mf.with_df, 'auxbasis', '?')) if df else \
           ('incore' if getattr(mf, '_eri', None) is not None else 'direto')
    no = mol.nelectron // 2
    homo, lumo = float(mf.mo_energy[no-1]), float(mf.mo_energy[no])
    reg.update(modo_eri=modo, converged=bool(mf.converged), level_shift=ls, e_tot_Ha=float(e),
               homo_eV=homo*vg.HARTREE_EV, lumo_eV=lumo*vg.HARTREE_EV, gap_eV=(lumo-homo)*vg.HARTREE_EV,
               tempo_s=round(dt, 1), obs='; '.join(obs))
    db['calculos'][chave] = reg; vg.salvar(db)
    print("     %s  gap=%.3f eV  E=%.6f Ha  conv=%s  %s  %.0f s" % (metodo, reg['gap_eV'], e, mf.converged, modo, dt), flush=True)
    return reg

CRIT_P = 'min_viz=2+preencher'

def etapa_e1(db):
    print("== Etapa (e1): min_viz=2 + preenchimento, D=1.2 nm, STO-3G exato ==", flush=True)
    si, h = construir_preencher(1.2)
    assert vg.formula(si, h) == 'Si39H40', vg.formula(si, h)
    for met in ('HF', 'PBE'):
        calcular_geom(db, 'e1', 1.2, 'mv2p', CRIT_P, si, h, 'sto-3g', met)
    vg.exportar_csv(db)

def etapa_e2(db):
    print("== Etapa (e2): erro do DF em STO-3G (Si36H40) ==", flush=True)
    si, h = vg.construir(1.2, 8, 2)
    for met in ('HF', 'PBE'):
        calcular_geom(db, 'e2', 1.2, 'mv2', 'min_viz=2', si, h, 'sto-3g', met, df=True)
    vg.exportar_csv(db)

def etapa_e3(db):
    print("== Etapa (e3): D=1.5 nm com DF, STO-3G ==", flush=True)
    si, h = vg.construir(1.5, 8, 2)
    for met in ('HF', 'PBE'):
        calcular_geom(db, 'e3', 1.5, 'mv2', 'min_viz=2', si, h, 'sto-3g', met, df=True)
    si, h = construir_preencher(1.5)
    assert vg.formula(si, h) == 'Si84H64', vg.formula(si, h)
    for met in ('HF', 'PBE'):
        calcular_geom(db, 'e3', 1.5, 'mv2p', CRIT_P, si, h, 'sto-3g', met, df=True)
    si, h = vg.construir(1.5, 8, 1)
    for met in ('HF', 'PBE'):
        calcular_geom(db, 'e3', 1.5, 'mv1', 'min_viz<=1', si, h, 'sto-3g', met, df=True)
    vg.exportar_csv(db)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--etapas', nargs='+', default=['e1', 'e2', 'e3'])
    args = ap.parse_args()
    db = vg.carregar()
    # etiqueta os registros antigos com o campo 'criterio'
    for k, r in db['calculos'].items():
        r.setdefault('criterio', {1: 'min_viz<=1', 2: 'min_viz=2', 3: 'min_viz=3 (vazio)'}.get(r.get('min_viz'), '?'))
    t0 = time.time()
    for et in args.etapas:
        {'e1': etapa_e1, 'e2': etapa_e2, 'e3': etapa_e3}[et](db)
        print("   tempo acumulado: %.0f s" % (time.time() - t0), flush=True)
    vg.exportar_csv(db)
    print("FIM. Total %.0f s" % (time.time() - t0), flush=True)
