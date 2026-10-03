#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
defeito_ligacao_pendente.py -- ligacao pendente de superficie criada pela quebra de uma ligacao Si-H,
o defeito tipico que a radiacao ionizante produz em nanocristais de Si passivados com H
(dessorcao de H por excitacao eletronica). Resultado novo da versao ARI do artigo.

Para um aglomerado perfeito (Si10H16 ou Si22H28, min_viz = 2), em PBE/3-21G com integrais exatas:
  1. relaxa o aglomerado perfeito (RKS, geomeTRIC, criterios padrao -- os mesmos da Secao 4.6);
  2. agrupa os H em classes de equivalencia por simetria (impressao digital de distancias na
     geometria ideal, que tem a simetria exata do corte);
  3. para cada classe, remove um H da geometria relaxada e calcula o radical (dubleto) SEM relaxar:
     UKS-PBE e UHF, niveis de fronteira por canal de spin, niveis da ligacao pendente (identificados
     pela localizacao no Si que perdeu o H), populacao de spin de Mulliken e <S^2>;
  4. relaxa o radical (UKS-PBE/3-21G) para as classes escolhidas e repete a analise; calcula a
     energia de dissociacao Si-H, E(radical) + E(H) - E(perfeito);
  5. (opcional) pontos simples em def2-SVP nas geometrias relaxadas, como teste de base.
O gap de cada canal de spin, E_g^s = LUMO_s - HOMO_s, e a grandeza comparavel ao gap HOMO-LUMO de
camada fechada usado no restante do artigo; o gap efetivo do radical e min(E_g^alfa, E_g^beta).

Grava dados/defeito_<rotulo>.json de forma incremental e os XYZ relaxados em dados/defeito_xyz/.
Uso:  python3 codigo/defeito_ligacao_pendente.py Si10H16 [--gpu] [--relaxar SiH|todas|nenhuma] [--def2svp]
"""
import os, sys, json, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
AQUI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, AQUI)
import relaxacao as rx                                   # metricas() e alinhar() (Kabsch)
from pyscf import gto, lib
from pyscf.geomopt import geometric_solver

HA_EV = 27.211386
RAIZ = os.path.dirname(AQUI); DADOS = os.path.join(RAIZ, 'dados')
DIR_XYZ = os.path.join(DADOS, 'defeito_xyz'); os.makedirs(DIR_XYZ, exist_ok=True)
GPU = '--gpu' in sys.argv
BASE = '3-21g'
ORIGEM = {'Si10H16': 'nanocristal_D0.8nm_minviz2.xyz', 'Si22H28': 'nanocristal_D1.0nm_minviz2.xyz'}
lib.num_threads(8)

# ------------------------------------------------------------------ utilitarios
def ler_xyz(p):
    L = open(p).read().split('\n'); n = int(L[0])
    at = [l.split() for l in L[2:2 + n]]
    return [a[0] for a in at], np.array([[float(x) for x in a[1:4]] for a in at])

def gravar_xyz(nome, sim, X, comentario):
    with open(os.path.join(DIR_XYZ, nome), 'w') as f:
        f.write('%d\n%s\n' % (len(sim), comentario))
        f.writelines('%-2s %12.6f %12.6f %12.6f\n' % (s, *x) for s, x in zip(sim, X))

def montar(sim, X, spin, base=BASE):
    mol = gto.M(atom=[(s, tuple(x)) for s, x in zip(sim, X)], basis=base, unit='Angstrom',
                spin=spin, verbose=0)
    mol.max_memory = 3000
    return mol

def novo_mf(mol, metodo):
    aberto = mol.spin != 0
    if GPU:
        from gpu4pyscf import scf as gscf, dft as gdft
        mf = (gdft.UKS if aberto else gdft.RKS)(mol, xc='pbe') if metodo == 'PBE' else (gscf.UHF if aberto else gscf.RHF)(mol)
    else:
        from pyscf import scf, dft
        if metodo == 'PBE':
            mf = (dft.UKS if aberto else dft.RKS)(mol); mf.xc = 'pbe'
        else:
            mf = (scf.UHF if aberto else scf.RHF)(mol)
    mf.conv_tol = 1e-9; mf.max_cycle = 200
    return mf

def resolver(mol, metodo, dm0=None):
    """SCF com uma segunda tentativa amortecida; devolve o objeto na CPU (para analise) e o tempo."""
    t0 = time.time()
    mf = novo_mf(mol, metodo); mf.kernel(dm0=dm0)
    if not mf.converged:
        mf = novo_mf(mol, metodo); mf.level_shift = 0.3; mf.damp = 0.3; mf.max_cycle = 400; mf.kernel(dm0=dm0)
    t = time.time() - t0
    return (mf.to_cpu() if GPU else mf), t

def analisar(mf, idx_si=None):
    """Niveis de fronteira (eV). Para o radical: niveis por canal de spin, niveis da ligacao pendente
    identificados pela maior populacao de Mulliken no Si que perdeu o H, e populacao de spin nesse Si."""
    mol = mf.mol
    r = dict(E_Ha=float(mf.e_tot), convergiu=bool(mf.converged))
    if mol.spin == 0:
        e, o = mf.mo_energy, mf.mo_occ
        r.update(HOMO_eV=float(e[o > 0].max() * HA_EV), LUMO_eV=float(e[o == 0].min() * HA_EV))
        r['gap_eV'] = r['LUMO_eV'] - r['HOMO_eV']
        return r
    S = mol.intor('int1e_ovlp'); a0, a1 = mol.aoslice_by_atom()[idx_si][2:4]
    def pop_si(c):                                     # populacao de Mulliken do orbital no Si
        return float((c * (S @ c))[a0:a1].sum())
    (ea, eb), (oa, ob), (ca, cb) = mf.mo_energy, mf.mo_occ, mf.mo_coeff
    occa, vira = np.where(oa > 0)[0], np.where(oa == 0)[0]
    occb, virb = np.where(ob > 0)[0], np.where(ob == 0)[0]
    cand_o = occa[-6:]; pa = [pop_si(ca[:, i]) for i in cand_o]; i_db_o = cand_o[int(np.argmax(pa))]
    cand_v = virb[:6];  pb = [pop_si(cb[:, i]) for i in cand_v]; i_db_v = cand_v[int(np.argmax(pb))]
    dma, dmb = mf.make_rdm1()
    spin_pop = float((((dma - dmb) @ S).diagonal())[a0:a1].sum())
    s2 = float(mf.spin_square()[0])
    r.update(
        HOMO_a_eV=float(ea[occa[-1]] * HA_EV), LUMO_a_eV=float(ea[vira[0]] * HA_EV),
        HOMO_b_eV=float(eb[occb[-1]] * HA_EV), LUMO_b_eV=float(eb[virb[0]] * HA_EV),
        DB_ocupado_eV=float(ea[i_db_o] * HA_EV), DB_vazio_eV=float(eb[i_db_v] * HA_EV),
        DB_ocupado_e_HOMO_a=bool(i_db_o == occa[-1]), DB_vazio_e_LUMO_b=bool(i_db_v == virb[0]),
        pop_Si_DB_ocupado=float(max(pa)), pop_Si_DB_vazio=float(max(pb)),
        pop_spin_Si=spin_pop, S2=s2)
    r['gap_alfa_eV'] = r['LUMO_a_eV'] - r['HOMO_a_eV']
    r['gap_beta_eV'] = r['LUMO_b_eV'] - r['HOMO_b_eV']
    r['gap_efetivo_eV'] = min(r['gap_alfa_eV'], r['gap_beta_eV'])
    r['E_v_eV'] = r['HOMO_b_eV']; r['E_c_eV'] = r['LUMO_a_eV']      # bordas: canais sem o eletron da ligacao
    r['gap_bordas_eV'] = r['E_c_eV'] - r['E_v_eV']
    r['DB_ocupado_acima_Ev_eV'] = r['DB_ocupado_eV'] - r['E_v_eV']
    r['DB_vazio_abaixo_Ec_eV'] = r['E_c_eV'] - r['DB_vazio_eV']
    return r

def geometria_si(X, sim, i_si):
    """Vizinhanca do Si que perdeu o H: distancias e soma dos tres angulos de ligacao
    (360 graus = planar; 328,4 graus = tetraedrico ideal)."""
    X = np.asarray(X); p = X[i_si]
    d = np.linalg.norm(X - p, axis=1)
    viz = [j for j in range(len(X)) if j != i_si and ((sim[j] == 'Si' and d[j] < 2.9) or (sim[j] == 'H' and d[j] < 1.9))]
    ang = []
    for a in range(len(viz)):
        for b in range(a + 1, len(viz)):
            u, v = X[viz[a]] - p, X[viz[b]] - p
            ang.append(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
    return dict(vizinhos=[sim[j] for j in viz], dist_A=[float(d[j]) for j in viz],
                soma_angulos_graus=float(sum(ang)) if len(viz) == 3 else None)

def relaxar(mol, rotulo_log):
    """Otimizacao completa (todos os atomos livres), UKS/RKS-PBE, criterios padrao do geomeTRIC."""
    traj = []; t0 = time.time()
    def cb(envs):
        g = envs['gradients']; g = g.get() if hasattr(g, 'get') else np.asarray(g)      # cupy -> numpy
        traj.append(dict(passo=len(traj) + 1, E_Ha=float(envs['energy']), gmax=float(np.abs(g).max())))
        print('      %s passo %3d  E=%.8f Ha  |g|max=%.2e  (%.0f s)' % (rotulo_log, len(traj), envs['energy'], traj[-1]['gmax'], time.time() - t0), flush=True)
    opt = geometric_solver.GeometryOptimizer(novo_mf(mol, 'PBE'))
    opt.callback = cb; opt.max_cycle = 100; opt.params = dict(assert_convergence=False)
    opt.kernel()
    return opt.mol.atom_coords(unit='Angstrom'), dict(convergiu=bool(opt.converged), n_passos=len(traj), tempo_s=time.time() - t0, trajetoria=traj)

# ------------------------------------------------------------------ principal
def main():
    rotulo = sys.argv[1]
    modo = sys.argv[sys.argv.index('--relaxar') + 1] if '--relaxar' in sys.argv else 'SiH'
    arq_json = os.path.join(DADOS, 'defeito_%s.json' % rotulo)
    dados = json.load(open(arq_json)) if os.path.exists(arq_json) else {}
    def salvar():
        tmp = arq_json + '.tmp'; json.dump(dados, open(tmp, 'w'), indent=1, ensure_ascii=False); os.replace(tmp, arq_json)
    import pyscf, geometric
    dados['_meta'] = dict(base=BASE, xc='pbe', integrais='exatas (SCF direto na GPU)' if GPU else 'exatas (CPU)',
                          dispositivo='GPU (gpu4pyscf %s)' % __import__('gpu4pyscf').__version__ if GPU else 'CPU, 8 fios',
                          pyscf=pyscf.__version__, geometric=geometric.__version__, conv_tol_Ha=1e-9,
                          otimizacao='geomeTRIC, criterios padrao (dE 1e-6 Ha, grms 3e-4, gmax 4.5e-4 Ha/Bohr, drms 1.2e-3, dmax 1.8e-3 A)',
                          definicoes='E_v = HOMO beta, E_c = LUMO alfa; DB = orbital mais localizado no Si que perdeu o H; '
                                     'gap_efetivo = min(gap alfa, gap beta)')
    sim0, X0 = ler_xyz(os.path.join(DADOS, 'xyz', ORIGEM[rotulo]))
    nsi = sim0.count('Si'); assert 'Si%dH%d' % (nsi, len(sim0) - nsi) == rotulo
    assert all(x == 'Si' for x in sim0[:nsi]), 'o XYZ deve listar os Si antes dos H'
    print('=' * 72, '\n%s  (%s)' % (rotulo, dados['_meta']['dispositivo']), flush=True)

    # 1. aglomerado perfeito relaxado
    if 'perfeito' not in dados:
        mol = montar(sim0, X0, 0)
        Xr, info = relaxar(mol, 'perfeito')
        Xr = rx.alinhar(Xr, Xr[:nsi], X0[:nsi])
        gravar_xyz('%s_relax_pbe_321g.xyz' % rotulo, sim0, Xr, '%s perfeito relaxado PBE/3-21G' % rotulo)
        molr = montar(sim0, Xr, 0)
        mfp, tp = resolver(molr, 'PBE'); mfh, th = resolver(molr, 'HF')
        dados['perfeito'] = dict(otimizacao=info, X=Xr.tolist(), PBE=analisar(mfp), HF=analisar(mfh), t_s=tp + th)
        salvar()
        print('   perfeito: gap PBE %.3f eV, HF %.3f eV (%d passos)' % (dados['perfeito']['PBE']['gap_eV'], dados['perfeito']['HF']['gap_eV'], info['n_passos']), flush=True)
    Xr = np.array(dados['perfeito']['X'])
    Ep = dados['perfeito']['PBE']['E_Ha']

    # atomo de H isolado (para a energia de dissociacao)
    if 'atomo_H' not in dados:
        molH = gto.M(atom='H 0 0 0', basis=BASE, spin=1, verbose=0)
        dados['atomo_H'] = {m: float(resolver(molH, m)[0].e_tot) for m in ('PBE', 'HF')}; salvar()
    EH = dados['atomo_H']

    # 2. classes de equivalencia dos H (geometria ideal = simetria exata do corte)
    si_idx = [i for i, s in enumerate(sim0) if s == 'Si']; h_idx = [i for i, s in enumerate(sim0) if s == 'H']
    c = X0[si_idx].mean(0)
    classes = {}
    for k in h_idx:
        fp = tuple(np.round(np.sort(np.linalg.norm(X0 - X0[k], axis=1)), 3))
        classes.setdefault(fp, []).append(k)
    classes = sorted(classes.values(), key=lambda ks: np.linalg.norm(X0[ks[0]] - c))
    dados.setdefault('sitios', {})
    for n_cl, ks in enumerate(classes, 1):
        k = ks[0]
        dono = si_idx[int(np.argmin([np.linalg.norm(X0[k] - X0[i]) for i in si_idx]))]
        nH = sum(1 for j in h_idx if np.linalg.norm(X0[j] - X0[dono]) < 1.6)
        nviz = sum(1 for i in si_idx if i != dono and np.linalg.norm(X0[i] - X0[dono]) < 2.6)
        chave = 'classe%d' % n_cl
        s = dados['sitios'].setdefault(chave, dict(n_equivalentes=len(ks), H_removido=k, Si_dono=dono,
                                                   tipo='SiH%d' % nH, vizinhos_Si=nviz,
                                                   dist_Si_centro_A=float(np.linalg.norm(X0[dono] - c))))
        # indices no sistema sem o H removido
        sim1 = [x for j, x in enumerate(sim0) if j != k]
        i_si = dono - (1 if dono > k else 0)
        # 3. radical sem relaxar
        if 'vertical' not in s:
            X1 = np.delete(Xr, k, axis=0); mol1 = montar(sim1, X1, 1)
            mfp, tp = resolver(mol1, 'PBE'); mfh, th = resolver(mol1, 'HF')
            s['vertical'] = dict(PBE=analisar(mfp, i_si), HF=analisar(mfh, i_si), geometria=geometria_si(X1, sim1, i_si), t_s=tp + th)
            s['vertical']['BDE_PBE_eV'] = (s['vertical']['PBE']['E_Ha'] + EH['PBE'] - Ep) * HA_EV
            salvar()
            v = s['vertical']['PBE']
            print('   %s (%s, %d equiv.) vertical: gap_ef PBE %.3f eV (alfa %.3f, beta %.3f), spin no Si %.2f, S2 %.3f; HF gap_ef %.3f'
                  % (chave, s['tipo'], s['n_equivalentes'], v['gap_efetivo_eV'], v['gap_alfa_eV'], v['gap_beta_eV'], v['pop_spin_Si'], v['S2'],
                     s['vertical']['HF']['gap_efetivo_eV']), flush=True)
        # 4. radical relaxado
        quer = modo == 'todas' or (modo == 'SiH' and s['tipo'] == 'SiH1')
        if quer and 'relaxado' not in s:
            X1 = np.delete(Xr, k, axis=0); mol1 = montar(sim1, X1, 1)
            X2, info = relaxar(mol1, chave)
            si1 = [j for j, x in enumerate(sim1) if x == 'Si']
            X2 = rx.alinhar(X2, X2[si1], X1[si1])
            gravar_xyz('%s_%s_%s_relax_pbe_321g.xyz' % (rotulo, chave, s['tipo']), sim1, X2,
                       '%s sem um H (%s, %s): radical relaxado UKS-PBE/3-21G' % (rotulo, chave, s['tipo']))
            mol2 = montar(sim1, X2, 1)
            mfp, tp = resolver(mol2, 'PBE'); mfh, th = resolver(mol2, 'HF')
            s['relaxado'] = dict(otimizacao=info, PBE=analisar(mfp, i_si), HF=analisar(mfh, i_si),
                                 geometria=geometria_si(X2, sim1, i_si),
                                 desloc_Si_dono_A=float(np.linalg.norm(X2[i_si] - X1[i_si])), X=X2.tolist(), t_s=tp + th)
            s['relaxado']['BDE_PBE_eV'] = (s['relaxado']['PBE']['E_Ha'] + EH['PBE'] - Ep) * HA_EV
            s['relaxado']['E_relaxacao_eV'] = (s['vertical']['PBE']['E_Ha'] - s['relaxado']['PBE']['E_Ha']) * HA_EV
            salvar()
            v = s['relaxado']['PBE']
            print('   %s relaxado (%d passos): gap_ef PBE %.3f eV, BDE %.2f eV, E_relax %.2f eV, soma de angulos %s'
                  % (chave, info['n_passos'], v['gap_efetivo_eV'], s['relaxado']['BDE_PBE_eV'], s['relaxado']['E_relaxacao_eV'],
                     s['relaxado']['geometria']['soma_angulos_graus']), flush=True)
        # 5. teste de base (def2-SVP) nas geometrias relaxadas
        if '--def2svp' in sys.argv and 'relaxado' in s and 'def2svp' not in s:
            mol2 = montar(sim1, np.array(s['relaxado']['X']), 1, base='def2-svp')
            mfp, tp = resolver(mol2, 'PBE')
            s['def2svp'] = dict(PBE=analisar(mfp, i_si), t_s=tp); salvar()
    if '--def2svp' in sys.argv and 'def2svp' not in dados['perfeito']:
        mfp, tp = resolver(montar(sim0, Xr, 0, base='def2-svp'), 'PBE')
        dados['perfeito']['def2svp'] = dict(PBE=analisar(mfp), t_s=tp); salvar()
    print('concluido:', arq_json, flush=True)

if __name__ == '__main__':
    main()
