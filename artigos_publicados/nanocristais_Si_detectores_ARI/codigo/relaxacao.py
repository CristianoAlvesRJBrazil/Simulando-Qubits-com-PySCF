#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
relaxacao.py -- Resultado R4 do artigo: o que muda quando a geometria ideal
(recortada da rede diamante perfeita) dos nanocristais Si_nH_m e relaxada.

Autocontido: reimplementa o gerador de nanocristais do Cap. 3 do livro
(rede diamante -> corte esferico -> remocao iterativa de Si com < min_viz
vizinhos Si -> passivacao com H a 1,48 A na direcao de cada vizinho removido).

Para cada sistema:
  1. geometria ideal: RHF/STO-3G e PBE/STO-3G (energia e gap HOMO-LUMO);
  2. otimizacao completa de geometria em PBE/STO-3G (geomeTRIC via PySCF);
  3. geometria relaxada: PBE (energia, gap) e RHF (energia, gap);
  4. metricas geometricas antes/depois: menor H-H, numero de pares H-H < 2,0 A,
     RMSD do nucleo de Si e da casca de H (alinhamento de Kabsch pelos Si),
     min/max Si-H e Si-Si, angulos H-Si-H dos grupos SiH2.
Resultados gravados de forma incremental em dados/relaxacao.json e os XYZ
(ideal e relaxado) em dados/.

Uso:  python3 relaxacao.py [rotulo ...]
      rotulos: Si10H16 Si19H28 Si22H28 Si36H40 (padrao: todos, nesta ordem)
"""
import os, sys, json, time
import numpy as np
from scipy.spatial.distance import pdist, squareform
from pyscf import gto, scf, dft, lib
from pyscf.geomopt import geometric_solver

lib.num_threads(6)
HA_EV = 27.211386
BASE = "sto-3g"
XC = "pbe"
MAX_MEMORY = 1500          # MB (maquina compartilhada)
RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DADOS = os.path.join(RAIZ, "dados")
os.makedirs(DADOS, exist_ok=True)
JSON_PATH = os.path.join(DADOS, "relaxacao.json")

# ------------------------------------------------------------------ gerador
A_SI, D_SIH, CORTE = 5.431, 1.48, 2.6

def rede_diamante(n_celulas):
    base = np.array([[0, 0, 0], [.25, .25, .25], [0, .5, .5], [.25, .75, .75],
                     [.5, 0, .5], [.75, .25, .75], [.5, .5, 0], [.75, .75, .25]])
    cel = np.array([[i, j, k] for i in range(n_celulas)
                    for j in range(n_celulas) for k in range(n_celulas)])
    return (cel[:, None, :] + base[None, :, :]).reshape(-1, 3) * A_SI

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

# ------------------------------------------------------------------ metricas
def kabsch(P, Q):
    """Rotacao/translacao que leva P sobre Q (minimos quadrados)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Pc, Qc

def alinhar(X, P, Q):
    """Aplica a X a transformacao que alinha P sobre Q."""
    R, Pc, Qc = kabsch(P, Q)
    return (X - Pc) @ R.T + Qc

def metricas(si, h):
    nsi, nh = len(si), len(h)
    dhh = pdist(h) if nh > 1 else np.array([np.inf])
    dsisi = pdist(si) if nsi > 1 else np.array([np.inf])
    dsih = np.linalg.norm(si[:, None, :] - h[None, :, :], axis=2)  # (nsi, nh)
    # ligacoes Si-H: H mais proximo de cada H (cada H tem um unico Si)
    lig_sih = dsih.min(axis=0)
    # ligacoes Si-Si: pares < 2.9 A (rede ideal: 2.352 A)
    lig_sisi = dsisi[dsisi < 2.9]
    # grupos SiHk e angulos H-Si-H
    grupos = {}
    angulos = []
    for i in range(nsi):
        js = np.where(dsih[i] < 1.9)[0]
        grupos[len(js)] = grupos.get(len(js), 0) + 1
        if len(js) == 2:
            v1 = h[js[0]] - si[i]; v2 = h[js[1]] - si[i]
            c = v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)
            angulos.append(float(np.degrees(np.arccos(np.clip(c, -1, 1)))))
    # numero de vizinhos Si de cada Si
    nviz = ((squareform(dsisi) < 2.9).sum(1) - 1) if nsi > 1 else np.array([0])
    return dict(
        formula="Si%dH%d" % (nsi, nh),
        min_HH=float(dhh.min()), n_HH_lt_2A=int((dhh < 2.0).sum()),
        n_HH_lt_1p5A=int((dhh < 1.5).sum()),
        min_SiH=float(lig_sih.min()), max_SiH=float(lig_sih.max()),
        min_SiSi=float(lig_sisi.min()), max_SiSi=float(lig_sisi.max()),
        n_lig_SiSi=int(len(lig_sisi)),
        coord_Si_total_min=int((nviz + (dsih < 1.9).sum(1)).min()),
        coord_Si_total_max=int((nviz + (dsih < 1.9).sum(1)).max()),
        grupos_SiHk={str(k): v for k, v in sorted(grupos.items())},
        ang_HSiH_min=float(min(angulos)) if angulos else None,
        ang_HSiH_max=float(max(angulos)) if angulos else None,
        ang_HSiH_media=float(np.mean(angulos)) if angulos else None,
        ang_HSiH=angulos,
        dist_HH_lt_3p5A=sorted(float(x) for x in dhh[dhh < 3.5]),
    )

# ------------------------------------------------------------------ PySCF
def montar_mol(si, h):
    atomos = [("Si", tuple(p)) for p in si] + [("H", tuple(p)) for p in h]
    mol = gto.M(atom=atomos, basis=BASE, unit="Angstrom", symmetry=False,
                verbose=0)
    mol.max_memory = MAX_MEMORY
    return mol

def gap_de(mf):
    e = mf.mo_energy; occ = mf.mo_occ
    homo = e[occ > 0].max(); lumo = e[occ == 0].min()
    return float(homo * HA_EV), float(lumo * HA_EV), float((lumo - homo) * HA_EV)

def calcular(mol, dm0=None):
    """RHF e PBE em STO-3G para a geometria de mol. Devolve dict e o mf PBE."""
    t0 = time.time()
    mf_hf = scf.RHF(mol); mf_hf.conv_tol = 1e-9
    mf_hf.kernel(dm0=dm0)
    t_hf = time.time() - t0
    h_hf, l_hf, g_hf = gap_de(mf_hf)
    t0 = time.time()
    mf = dft.RKS(mol); mf.xc = XC; mf.conv_tol = 1e-9
    mf.kernel(dm0=mf_hf.make_rdm1())
    t_pbe = time.time() - t0
    h, l, g = gap_de(mf)
    return dict(E_HF_Ha=float(mf_hf.e_tot), E_HF_eV=float(mf_hf.e_tot * HA_EV),
                HOMO_HF_eV=h_hf, LUMO_HF_eV=l_hf, gap_HF_eV=g_hf,
                conv_HF=bool(mf_hf.converged), t_HF_s=t_hf,
                E_PBE_Ha=float(mf.e_tot), E_PBE_eV=float(mf.e_tot * HA_EV),
                HOMO_PBE_eV=h, LUMO_PBE_eV=l, gap_PBE_eV=g,
                conv_PBE=bool(mf.converged), t_PBE_s=t_pbe), mf

def gravar_xyz(nome, si, h, comentario=""):
    L = ["%d" % (len(si) + len(h)), comentario]
    L += ["Si %12.6f %12.6f %12.6f" % tuple(p) for p in si]
    L += ["H  %12.6f %12.6f %12.6f" % tuple(p) for p in h]
    with open(os.path.join(DADOS, nome), "w") as f:
        f.write("\n".join(L) + "\n")

def carregar_json():
    if os.path.exists(JSON_PATH):
        with open(JSON_PATH) as f:
            return json.load(f)
    return {"_meta": dict(base=BASE, xc=XC, otimizador="geomeTRIC 1.1.1 (PySCF geomopt)",
                          pyscf=__import__("pyscf").__version__,
                          threads=6, max_memory_MB=MAX_MEMORY,
                          criterio="convergencia padrao do PySCF/geomeTRIC "
                                   "(dE 1e-6 Ha, grms 3e-4, gmax 4.5e-4 Ha/Bohr, "
                                   "drms 1.2e-3, dmax 1.8e-3 A)"),
            "sistemas": {}}

def salvar_json(dados):
    tmp = JSON_PATH + ".tmp"
    with open(tmp, "w") as f:
        json.dump(dados, f, indent=1, ensure_ascii=False)
    os.replace(tmp, JSON_PATH)

# ------------------------------------------------------------------ principal
SISTEMAS = {
    "Si10H16": dict(D=0.8, min_viz=2),
    "Si19H28": dict(D=0.9, min_viz=2),
    "Si22H28": dict(D=1.0, min_viz=2),
    "Si36H40": dict(D=1.2, min_viz=2),
}

def rodar(rotulo, maxsteps=100):
    par = SISTEMAS[rotulo]
    si0, h0 = construir(par["D"], min_viz=par["min_viz"])
    nsi, nh = len(si0), len(h0)
    assert "Si%dH%d" % (nsi, nh) == rotulo, (nsi, nh)
    print("=" * 70); print("%s  D=%.1f nm  min_viz=%d  (%d atomos)" %
                           (rotulo, par["D"], par["min_viz"], nsi + nh)); sys.stdout.flush()
    gravar_xyz(rotulo + "_ideal.xyz", si0, h0,
               "%s ideal: rede diamante a=%.3f A, corte D=%.1f nm, min_viz=%d, Si-H=%.2f A"
               % (rotulo, A_SI, par["D"], par["min_viz"], D_SIH))
    reg = dict(D_nm=par["D"], min_viz=par["min_viz"], n_Si=nsi, n_H=nh,
               n_atomos=nsi + nh)
    reg["geom_ideal"] = metricas(si0, h0)

    # --- geometria ideal
    mol = montar_mol(si0, h0)
    reg["n_AOs"] = mol.nao_nr(); reg["n_eletrons"] = mol.nelectron
    print("  AOs = %d" % mol.nao_nr()); sys.stdout.flush()
    res0, mf0 = calcular(mol)
    reg["ideal"] = res0
    print("  ideal : E_HF=%.6f Ha gap_HF=%.2f eV | E_PBE=%.6f Ha gap_PBE=%.2f eV (%.0f s + %.0f s)"
          % (res0["E_HF_Ha"], res0["gap_HF_eV"], res0["E_PBE_Ha"], res0["gap_PBE_eV"],
             res0["t_HF_s"], res0["t_PBE_s"])); sys.stdout.flush()

    # --- otimizacao PBE/STO-3G (todos os atomos livres)
    traj = []
    def cb(envs):
        traj.append(dict(passo=len(traj), E_Ha=float(envs["energy"]),
                         gmax=float(np.abs(envs["gradients"]).max())))
        print("    passo %3d  E=%.8f Ha  |g|max=%.2e  (%.0f s)"
              % (len(traj), envs["energy"], np.abs(envs["gradients"]).max(),
                 time.time() - t_opt)); sys.stdout.flush()
    t_opt = time.time()
    mf_opt = dft.RKS(mol); mf_opt.xc = XC; mf_opt.conv_tol = 1e-9
    opt = geometric_solver.GeometryOptimizer(mf_opt)
    opt.callback = cb
    opt.max_cycle = maxsteps
    opt.params = dict(assert_convergence=False)
    opt.kernel()
    t_opt = time.time() - t_opt
    mol_rel = opt.mol
    reg["otimizacao"] = dict(convergiu=bool(opt.converged), n_passos=len(traj),
                             tempo_s=t_opt, trajetoria=traj)
    print("  otimizacao: %s em %d passos, %.0f s" %
          ("convergiu" if opt.converged else "NAO convergiu", len(traj), t_opt))
    sys.stdout.flush()

    # --- geometria relaxada
    xyz = mol_rel.atom_coords(unit="Angstrom")
    si1, h1 = xyz[:nsi], xyz[nsi:]
    # alinhamento de Kabsch pelo nucleo de Si (remove translacao/rotacao global)
    si1a = alinhar(si1, si1, si0); h1a = alinhar(h1, si1, si0)
    reg["geom_relax"] = metricas(si1a, h1a)
    reg["geom_relax"]["RMSD_Si_A"] = float(np.sqrt(((si1a - si0) ** 2).sum(1).mean()))
    reg["geom_relax"]["RMSD_H_A"] = float(np.sqrt(((h1a - h0) ** 2).sum(1).mean()))
    reg["geom_relax"]["desloc_max_Si_A"] = float(np.linalg.norm(si1a - si0, axis=1).max())
    reg["geom_relax"]["desloc_max_H_A"] = float(np.linalg.norm(h1a - h0, axis=1).max())
    # raio de giro (medida do 'inchaco'/'encolhimento' do nucleo)
    reg["geom_ideal"]["Rg_Si_A"] = float(np.sqrt(((si0 - si0.mean(0)) ** 2).sum(1).mean()))
    reg["geom_relax"]["Rg_Si_A"] = float(np.sqrt(((si1a - si1a.mean(0)) ** 2).sum(1).mean()))
    gravar_xyz(rotulo + "_relax_pbe_sto3g.xyz", si1a, h1a,
               "%s relaxado PBE/STO-3G (geomeTRIC), alinhado ao nucleo de Si ideal" % rotulo)
    mol1 = montar_mol(si1a, h1a)
    res1, mf1 = calcular(mol1, dm0=None)
    reg["relax"] = res1
    reg["delta"] = dict(
        dE_PBE_Ha=res1["E_PBE_Ha"] - res0["E_PBE_Ha"],
        dE_PBE_eV=(res1["E_PBE_Ha"] - res0["E_PBE_Ha"]) * HA_EV,
        dE_PBE_meV_por_H=(res1["E_PBE_Ha"] - res0["E_PBE_Ha"]) * HA_EV * 1000 / nh,
        dE_HF_Ha=res1["E_HF_Ha"] - res0["E_HF_Ha"],
        dE_HF_eV=(res1["E_HF_Ha"] - res0["E_HF_Ha"]) * HA_EV,
        dgap_PBE_eV=res1["gap_PBE_eV"] - res0["gap_PBE_eV"],
        dgap_HF_eV=res1["gap_HF_eV"] - res0["gap_HF_eV"],
    )
    print("  relax : E_HF=%.6f Ha gap_HF=%.2f eV | E_PBE=%.6f Ha gap_PBE=%.2f eV"
          % (res1["E_HF_Ha"], res1["gap_HF_eV"], res1["E_PBE_Ha"], res1["gap_PBE_eV"]))
    print("  dE_PBE = %.4f eV ; dgap_PBE = %+.3f eV ; dgap_HF = %+.3f eV"
          % (reg["delta"]["dE_PBE_eV"], reg["delta"]["dgap_PBE_eV"], reg["delta"]["dgap_HF_eV"]))
    g0, g1 = reg["geom_ideal"], reg["geom_relax"]
    print("  minHH %.3f -> %.3f A ; pares<2A %d -> %d ; RMSD Si %.3f A, H %.3f A"
          % (g0["min_HH"], g1["min_HH"], g0["n_HH_lt_2A"], g1["n_HH_lt_2A"],
             g1["RMSD_Si_A"], g1["RMSD_H_A"]))
    print("  Si-H %.3f-%.3f ; Si-Si %.3f-%.3f ; H-Si-H %.1f-%.1f (media %.1f)"
          % (g1["min_SiH"], g1["max_SiH"], g1["min_SiSi"], g1["max_SiSi"],
             g1["ang_HSiH_min"] or 0, g1["ang_HSiH_max"] or 0, g1["ang_HSiH_media"] or 0))
    sys.stdout.flush()
    return reg

if __name__ == "__main__":
    rotulos = sys.argv[1:] or list(SISTEMAS)
    for r in rotulos:
        t0 = time.time()
        reg = rodar(r)
        reg["tempo_total_s"] = time.time() - t0
        dados = carregar_json()
        dados["sistemas"][r] = reg
        salvar_json(dados)
        print("  [%s gravado em %s; %.0f s]" % (r, JSON_PATH, reg["tempo_total_s"]))
