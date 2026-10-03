#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
varredura_gap.py -- Resultado R3 do artigo (Secoes 4.3 e 4.5):
efeito do criterio de corte (min_viz) e da base sobre o gap HOMO-LUMO
de nanocristais de Si:H construidos a partir da rede diamante ideal.

Autocontido: copia as funcoes do gerador (codigo/nanocristal.py do livro).
Grava resultados incrementalmente em ../dados/gaps.json (um registro por
calculo, sobrevive a interrupcoes) e, ao final de cada etapa, ../dados/gaps.csv.

Uso:  python3 varredura_gap.py --etapas a b c d
      a = validacao (min_viz=2, STO-3G, HF/PBE, D=0.8-1.2 nm; reproduz o livro)
      b = sensibilidade ao corte (min_viz=1, STO-3G, HF/PBE, mesmos D)
      c = base (Si10H16 e Si22H28: STO-3G, 3-21G, def2-SVP; HF/PBE)
      d = D=1.5 nm (Si72H64 e Si84H88), STO-3G, HF/PBE
"""
import os, sys, time, json, argparse, csv, platform
import numpy as np
from pyscf import gto, scf, dft, lib
import pyscf

N_THREADS = 8
MAX_MEMORY = 2500          # MB por objeto Mole (maquina compartilhada)
HARTREE_EV = 27.2114
lib.num_threads(N_THREADS)

AQUI = os.path.dirname(os.path.abspath(__file__))
DADOS = os.path.join(AQUI, '..', 'dados')
os.makedirs(DADOS, exist_ok=True)
ARQ_JSON = os.path.join(DADOS, 'gaps.json')
ARQ_CSV = os.path.join(DADOS, 'gaps.csv')

# ----------------------------------------------------------------------
# Gerador de nanocristais (copiado de codigo/nanocristal.py do livro)
# ----------------------------------------------------------------------
A_SI, D_SIH, CORTE = 5.431, 1.48, 2.6

def rede_diamante(n_celulas):
    base = np.array([[0, 0, 0], [.25, .25, .25], [0, .5, .5], [.25, .75, .75],
                     [.5, 0, .5], [.75, .25, .75], [.5, .5, 0], [.75, .75, .25]])
    celulas = np.array([[i, j, k] for i in range(n_celulas)
                        for j in range(n_celulas) for k in range(n_celulas)])
    return (celulas[:, None, :] + base[None, :, :]).reshape(-1, 3) * A_SI

def construir(diametro_nm, n_celulas=8, min_viz=2):
    raio = diametro_nm * 10 / 2
    rede = rede_diamante(n_celulas)
    rede = rede - rede.mean(axis=0)
    dentro = np.linalg.norm(rede, axis=1) <= raio
    while True:
        idx = np.where(dentro)[0]
        pos = rede[idx]
        nviz = [np.sum((np.linalg.norm(pos - p, axis=1) > 0.1) &
                       (np.linalg.norm(pos - p, axis=1) < CORTE)) for p in pos]
        ruins = idx[np.array(nviz) < min_viz]
        if len(ruins) == 0:
            break
        dentro[ruins] = False
    si = rede[dentro]
    H = []
    for p in si:
        d = np.linalg.norm(rede - p, axis=1)
        for j in np.where((d > 0.1) & (d < CORTE))[0]:
            if not dentro[j]:
                u = (rede[j] - p) / d[j]
                H.append(p + D_SIH * u)
    return si, np.array(H).reshape(-1, 3)

def para_xyz(si, h, idx_p=None):
    L = []
    for i, p in enumerate(si):
        L.append("%-2s %.6f %.6f %.6f" % ("P" if i == idx_p else "Si", *p))
    for p in h:
        L.append("H %.6f %.6f %.6f" % tuple(p))
    return "\n".join(L)

# ----------------------------------------------------------------------
# Equacao de Brus (parametros do Cap. 3 do livro)
# ----------------------------------------------------------------------
HB2_2M = 3.80998   # hbar^2/(2 m0), em eV.A^2
E2_4PE = 14.3996   # e^2/(4 pi eps0), em eV.A
EG_BULK, ME, MH, EPS = 1.12, 0.26, 0.38, 11.7

def gap_brus(diametro_nm):
    R = diametro_nm * 10 / 2
    confinamento = HB2_2M * np.pi**2 / R**2 * (1/ME + 1/MH)
    coulomb = -1.786 * E2_4PE / (EPS * R)
    return EG_BULK + confinamento + coulomb

def gap_delerue(diametro_nm):
    """Ajuste de Delerue, Allan e Lannoo, PRB 48, 11024 (1993): E_g = 1.167 + 3.73/D^1.39 (D em nm)."""
    return 1.167 + 3.73 / diametro_nm**1.39

# ----------------------------------------------------------------------
# Utilidades estruturais
# ----------------------------------------------------------------------
def formula(si, h):
    return "Si%dH%d" % (len(si), len(h))

def diametro_efetivo_nm(n_si):
    """Diametro da esfera com o mesmo volume que n_si atomos ocupam no Si macico
    (8 atomos por celula cubica de aresta A_SI): D = a (3N/4pi)^(1/3)."""
    return A_SI * (3.0 * n_si / (4.0 * np.pi))**(1.0/3.0) / 10.0

def diagnostico(si, h):
    d = {}
    if len(si) == 0:
        return dict(n_si=0, n_h=0)
    dd = np.linalg.norm(si[:, None, :] - si[None, :, :], axis=2)
    nv = ((dd > 0.1) & (dd < CORTE)).sum(axis=1)
    d['n_si'] = int(len(si)); d['n_h'] = int(len(h))
    d['coord_si'] = {int(k): int((nv == k).sum()) for k in np.unique(nv)}
    d['n_SiH3'] = int((nv == 1).sum()); d['n_SiH2'] = int((nv == 2).sum()); d['n_SiH'] = int((nv == 3).sum())
    dd[dd < 0.1] = 99
    d['min_SiSi_A'] = float(dd.min())
    if len(h):
        dsh = np.linalg.norm(si[:, None, :] - h[None, :, :], axis=2)
        d['min_SiH_A'] = float(dsh.min())
        if len(h) > 1:
            dhh = np.linalg.norm(h[:, None, :] - h[None, :, :], axis=2)
            dhh[dhh < 0.1] = 99
            d['min_HH_A'] = float(dhh.min())
            d['n_pares_HH_lt_1p6A'] = int((dhh < 1.6).sum() // 2)
    return d

# ----------------------------------------------------------------------
# Persistencia
# ----------------------------------------------------------------------
def carregar():
    if os.path.exists(ARQ_JSON):
        with open(ARQ_JSON) as f:
            return json.load(f)
    return {"meta": {"pyscf": pyscf.__version__, "numpy": np.__version__,
                     "python": platform.python_version(), "threads": N_THREADS,
                     "max_memory_MB": MAX_MEMORY, "conv_tol": 1e-8,
                     "geometria": "ideal (rede diamante a=5.431 A, Si-H=1.48 A, nao relaxada)",
                     "gap_def": "(e_LUMO - e_HOMO) * 27.2114 eV, camada fechada (RHF/RKS)"},
            "estruturas": {}, "calculos": {}}

def salvar(db):
    tmp = ARQ_JSON + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(db, f, indent=1, ensure_ascii=False)
    os.replace(tmp, ARQ_JSON)

def exportar_csv(db):
    campos = ['chave', 'etapa', 'D_nominal_nm', 'min_viz', 'formula', 'n_si', 'n_h', 'D_efetivo_nm',
              'base', 'metodo', 'modo_eri', 'nao', 'nelectron', 'converged', 'ciclos',
              'level_shift', 'e_tot_Ha', 'homo_eV', 'lumo_eV', 'gap_eV', 'tempo_s', 'gap_brus_Dnominal_eV',
              'gap_brus_Defetivo_eV', 'gap_delerue_Defetivo_eV', 'obs']
    with open(ARQ_CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        w.writeheader()
        for k, r in db['calculos'].items():
            r2 = dict(r); r2['chave'] = k
            w.writerow(r2)

# ----------------------------------------------------------------------
# Calculo SCF
# ----------------------------------------------------------------------
def calcular(db, etapa, D, min_viz, base, metodo, df=False, forcar=False):
    si, h = construir(D, 8, min_viz)
    fm = formula(si, h)
    chave = "D%.1f|mv%d|%s|%s|%s|%s" % (D, min_viz, fm, base, metodo, 'DF' if df else 'exato')
    if chave in db['calculos'] and not forcar and db['calculos'][chave].get('converged'):
        print("  [ja feito] " + chave, flush=True)
        return db['calculos'][chave]
    ekey = "D%.1f|mv%d" % (D, min_viz)
    if ekey not in db['estruturas']:
        dg = diagnostico(si, h); dg['formula'] = fm; dg['D_nominal_nm'] = D; dg['min_viz'] = min_viz
        dg['D_efetivo_nm'] = diametro_efetivo_nm(len(si)) if len(si) else None
        db['estruturas'][ekey] = dg
    if len(si) == 0:
        reg = dict(etapa=etapa, D_nominal_nm=D, min_viz=min_viz, formula=fm, n_si=0, n_h=0,
                   base=base, metodo=metodo, converged=False, obs='estrutura vazia (min_viz>=3 erode todo o fragmento)')
        db['calculos'][chave] = reg; salvar(db); return reg
    mol = gto.M(atom=para_xyz(si, h), basis=base, unit='Angstrom', verbose=0,
                max_memory=MAX_MEMORY, symmetry=False)
    reg = dict(etapa=etapa, D_nominal_nm=D, min_viz=min_viz, formula=fm, n_si=len(si), n_h=len(h),
               D_efetivo_nm=diametro_efetivo_nm(len(si)), base=base, metodo=metodo,
               nao=int(mol.nao), nelectron=int(mol.nelectron), threads=N_THREADS,
               gap_brus_Dnominal_eV=float(gap_brus(D)),
               gap_brus_Defetivo_eV=float(gap_brus(diametro_efetivo_nm(len(si)))),
               gap_delerue_Defetivo_eV=float(gap_delerue(diametro_efetivo_nm(len(si)))))
    if mol.nelectron % 2:
        reg.update(converged=False, obs='numero impar de eletrons: camada aberta, pulado')
        db['calculos'][chave] = reg; salvar(db); return reg
    print("  -> %s  nao=%d  nel=%d" % (chave, mol.nao, mol.nelectron), flush=True)
    obs = []
    for tentativa, ls in enumerate((0.0, 0.3)):
        if metodo.upper() == 'HF':
            mf = scf.RHF(mol)
        else:
            mf = dft.RKS(mol); mf.xc = metodo.lower()
        if df:
            mf = mf.density_fit()
        mf.conv_tol = 1e-8; mf.max_cycle = 150; mf.level_shift = ls
        mf.verbose = 0
        t0 = time.time()
        try:
            e = mf.kernel()
        except Exception as exc:
            reg.update(converged=False, obs='excecao: %r' % exc); db['calculos'][chave] = reg; salvar(db)
            print("  !! excecao", exc, flush=True); return reg
        dt = time.time() - t0
        if mf.converged:
            break
        obs.append('nao convergiu com level_shift=%.1f' % ls)
    if df:
        modo = 'DF (auxbase %s)' % getattr(mf.with_df, 'auxbasis', '?')
    else:
        modo = 'incore' if getattr(mf, '_eri', None) is not None else 'direto'
    no = mol.nelectron // 2
    homo, lumo = float(mf.mo_energy[no-1]), float(mf.mo_energy[no])
    ciclos = None
    reg.update(modo_eri=modo, converged=bool(mf.converged), level_shift=ls, e_tot_Ha=float(e),
               homo_eV=homo*HARTREE_EV, lumo_eV=lumo*HARTREE_EV, gap_eV=(lumo-homo)*HARTREE_EV,
               tempo_s=round(dt, 1), obs='; '.join(obs))
    db['calculos'][chave] = reg; salvar(db)
    print("     %s  gap=%.3f eV  E=%.6f Ha  conv=%s  %s  %.0f s" %
          (metodo, reg['gap_eV'], e, mf.converged, modo, dt), flush=True)
    return reg

# ----------------------------------------------------------------------
# Etapas
# ----------------------------------------------------------------------
def etapa_a(db):
    print("== Etapa (a): validacao min_viz=2 / STO-3G ==", flush=True)
    for D in (0.8, 0.9, 1.0, 1.2):
        for met in ('HF', 'PBE'):
            calcular(db, 'a', D, 2, 'sto-3g', met)
    exportar_csv(db)

def etapa_b(db):
    print("== Etapa (b): sensibilidade ao corte (min_viz=1 e 3) / STO-3G ==", flush=True)
    for D in (0.8, 0.9, 1.0, 1.2):
        calcular(db, 'b', D, 3, 'sto-3g', 'HF')     # registra estrutura vazia
    for D in (0.8, 0.9, 1.0, 1.2):
        for met in ('HF', 'PBE'):
            calcular(db, 'b', D, 1, 'sto-3g', met)
    exportar_csv(db)

def etapa_c(db):
    print("== Etapa (c): base (Si10H16 e Si22H28, min_viz=2) ==", flush=True)
    # Si10H16 (D=0.8): tudo sem DF (barato); def2-SVP tambem com DF para medir o erro do DF
    for base in ('sto-3g', '3-21g', 'def2-svp'):
        for met in ('HF', 'PBE'):
            calcular(db, 'c', 0.8, 2, base, met, df=False)
    for met in ('HF', 'PBE'):
        calcular(db, 'c', 0.8, 2, 'def2-svp', met, df=True)
    # Si22H28 (D=1.0): STO-3G exato (livro); 3-21G e def2-SVP com DF
    for met in ('HF', 'PBE'):
        calcular(db, 'c', 1.0, 2, 'sto-3g', met, df=False)
    for base in ('3-21g', 'def2-svp'):
        for met in ('HF', 'PBE'):
            calcular(db, 'c', 1.0, 2, base, met, df=True)
    exportar_csv(db)

def etapa_c2(db):
    # complementar: Si22H28 3-21G sem DF (verificacao do erro de DF em cluster maior)
    print("== Etapa (c2): Si22H28 3-21G exato ==", flush=True)
    for met in ('HF', 'PBE'):
        calcular(db, 'c', 1.0, 2, '3-21g', met, df=False)
    exportar_csv(db)

def etapa_d(db):
    print("== Etapa (d): D=1.5 nm, STO-3G ==", flush=True)
    for met in ('HF', 'PBE'):
        calcular(db, 'd', 1.5, 2, 'sto-3g', met)
    for met in ('HF', 'PBE'):
        calcular(db, 'd', 1.5, 1, 'sto-3g', met)
    exportar_csv(db)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--etapas', nargs='+', default=['a', 'b', 'c', 'd'])
    args = ap.parse_args()
    db = carregar()
    t0 = time.time()
    for et in args.etapas:
        {'a': etapa_a, 'b': etapa_b, 'c': etapa_c, 'c2': etapa_c2, 'd': etapa_d}[et](db)
        print("   tempo acumulado: %.0f s" % (time.time() - t0), flush=True)
    exportar_csv(db)
    print("FIM. Total %.0f s" % (time.time() - t0), flush=True)
