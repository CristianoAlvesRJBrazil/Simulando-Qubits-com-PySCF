#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""varredura_gpu.py -- sistemas de 1,5 nm (Si72H64, Si84H64, Si84H88) em HF e PBE / STO-3G com
ajuste de densidade executado na GPU (gpu4pyscf). Mesmo formato de registro de varredura_gap.py;
grava em ../dados/gaps.json (chaves 'D1.5|<tag>|<formula>|sto-3g|<metodo>|DF').
Primeiro job: validacao em Si10H16 (compara com CPU: HF 16.823 / PBE 8.290 eV, integrais exatas).
Uso: python3 varredura_gpu.py"""
import os, sys, time, json
import numpy as np
import cupy
from pyscf import gto, lib
from gpu4pyscf import scf as gscf, dft as gdft
import varredura_gap as vg
from varredura_extra import construir_preencher
lib.num_threads(8)
GPU = cupy.cuda.runtime.getDeviceProperties(0)['name'].decode()

def calcular_gpu(db, etapa, D, tag, criterio, si, h, base, metodo, auxbasis=None, direto=False, dm0=None):
    fm = vg.formula(si, h)
    chave = "D%.1f|%s|%s|%s|%s|%s" % (D, tag, fm, base, metodo, "exato" if direto else "DF")
    if chave in db['calculos'] and db['calculos'][chave].get('converged'):
        print("  [ja feito] " + chave, flush=True); return db['calculos'][chave]
    ekey = "D%.1f|%s" % (D, tag)
    if ekey not in db['estruturas']:
        dg = vg.diagnostico(si, h); dg['formula'] = fm; dg['D_nominal_nm'] = D
        dg['min_viz'] = (2 if tag == 'mv2' else (1 if tag == 'mv1' else '2p')); dg['criterio'] = criterio
        dg['D_efetivo_nm'] = vg.diametro_efetivo_nm(len(si)); db['estruturas'][ekey] = dg
    mol = gto.M(atom=vg.para_xyz(si, h), basis=base, unit='Angstrom', verbose=0, max_memory=6000, symmetry=False)
    Def = vg.diametro_efetivo_nm(len(si))
    reg = dict(etapa=etapa, D_nominal_nm=D, min_viz=(2 if tag == 'mv2' else (1 if tag == 'mv1' else '2p')),
               criterio=criterio, formula=fm, n_si=len(si), n_h=len(h), D_efetivo_nm=Def, base=base, metodo=metodo,
               nao=int(mol.nao), nelectron=int(mol.nelectron), threads=8, gpu=GPU,
               gap_brus_Dnominal_eV=float(vg.gap_brus(D)), gap_brus_Defetivo_eV=float(vg.gap_brus(Def)),
               gap_delerue_Defetivo_eV=float(vg.gap_delerue(Def)))
    if mol.nelectron % 2:
        reg.update(converged=False, obs='numero impar de eletrons'); db['calculos'][chave] = reg; vg.salvar(db); return reg
    print("  -> %s  nao=%d  nel=%d  [GPU]" % (chave, mol.nao, mol.nelectron), flush=True)
    obs = []
    for ls in (0.0, 0.3):
        mf = gscf.RHF(mol) if metodo.upper() == 'HF' else gdft.RKS(mol, xc=metodo.lower())
        if direto:
            # SCF direto com screening 1e-10: o ruido da energia (~1e-5 Ha em 712+ funcoes) impede conv_tol=1e-7;
            # o gap depende do gradiente, que converge a ~1e-6. Criterios: dE < 2e-5 Ha e |g| < 1e-4.
            mf.direct_scf_tol = 1e-10; mf.init_guess = 'atom'; mf.conv_tol = 2e-5; mf.conv_tol_grad = 1e-4
        else:
            mf = mf.density_fit(auxbasis=auxbasis) if auxbasis else mf.density_fit(); mf.conv_tol = 1e-8
        mf.max_cycle = 150; mf.verbose = 4 if direto else 0
        if ls > 0: mf.level_shift = ls
        t0 = time.time()
        try:
            e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
        except Exception as exc:
            reg.update(converged=False, obs='excecao: %r' % exc); db['calculos'][chave] = reg; vg.salvar(db)
            print("  !! excecao", repr(exc)[:300], flush=True); return reg
        dt = time.time() - t0
        if mf.converged: break
        obs.append('nao convergiu com level_shift=%.1f' % ls)
    aux = '?' if direto else getattr(mf.with_df, 'auxbasis', '?')
    moe = cupy.asnumpy(mf.mo_energy) if hasattr(mf.mo_energy, 'get') else np.asarray(mf.mo_energy)
    no = mol.nelectron // 2
    homo, lumo = float(moe[no-1]), float(moe[no])
    reg.update(conv_tol=mf.conv_tol, conv_tol_grad=getattr(mf, 'conv_tol_grad', None), modo_eri=('direto-GPU (screening 1e-10)' if direto else 'DF-GPU (auxbase %s)' % aux), converged=bool(mf.converged), level_shift=ls, e_tot_Ha=float(e),
               homo_eV=homo*vg.HARTREE_EV, lumo_eV=lumo*vg.HARTREE_EV, gap_eV=(lumo-homo)*vg.HARTREE_EV,
               tempo_s=round(dt, 1), obs='; '.join(obs))
    db['calculos'][chave] = reg; vg.salvar(db)
    print("     %s  gap=%.3f eV  E=%.6f Ha  conv=%s  %s  %.0f s" % (metodo, reg['gap_eV'], e, mf.converged, reg['modo_eri'], dt), flush=True)
    return reg

if __name__ == '__main__':
    db = vg.carregar(); t0 = time.time()
    print("GPU:", GPU, flush=True)
    print("== (g0) validacao GPU: Si10H16 STO-3G DF ==", flush=True)
    si, h = vg.construir(0.8, 8, 2)
    for met in ('HF', 'PBE'):
        calcular_gpu(db, 'g0', 0.8, 'mv2', 'min_viz=2', si, h, 'sto-3g', met)
    print("== (g3) D=1.5 nm, STO-3G, SCF direto na GPU ==", flush=True)
    for tag, crit, (si, h) in (('mv2', 'min_viz=2', vg.construir(1.5, 8, 2)),
                               ('mv2p', 'min_viz=2+preencher', construir_preencher(1.5)),
                               ('mv1', 'min_viz<=1', vg.construir(1.5, 8, 1))):
        calcular_gpu(db, 'g3', 1.5, tag, crit, si, h, 'sto-3g', 'HF', direto=True)
        # PBE parte do palpite atomico: a densidade HF como ponto de partida tornou o SCF PBE muito lento
        calcular_gpu(db, 'g3', 1.5, tag, crit, si, h, 'sto-3g', 'PBE', direto=True)
        vg.exportar_csv(db)
    vg.exportar_csv(db)
    print("FIM. Total %.0f s" % (time.time() - t0), flush=True)
