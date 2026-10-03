#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura R4 -- relaxacao de geometria (PBE, geomeTRIC).
(a) menor distancia H...H de cada hidrogenio, ordenada, na geometria ideal (circulos vazios) e
    relaxada em PBE/STO-3G (triangulos cheios), para cada sistema em dados/relaxacao.json
    (le os XYZ gravados por relaxacao.py);
(b) variacao do gap HOMO-LUMO (relaxada - ideal), HF e PBE, STO-3G e, se existir, controle 3-21G.
Saida: ../figuras/fig_relaxacao.png (300 dpi, 3,35 in, P&B)."""
import os, json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estilo as ST
AQUI = os.path.dirname(os.path.abspath(__file__)); RAIZ = os.path.join(AQUI, '..'); DADOS = os.path.join(RAIZ, 'dados')
def ler_xyz(p):
    L = open(p).read().split('\n'); n = int(L[0]); at = [l.split() for l in L[2:2 + n]]
    sim = np.array([a[0] for a in at]); xyz = np.array([[float(x) for x in a[1:4]] for a in at]); return sim, xyz
def dmin_por_H(p):
    sim, xyz = ler_xyz(p); H = xyz[sim == 'H']
    d = np.linalg.norm(H[:, None] - H[None], axis=2); d[d < 0.1] = 99
    return np.sort(d.min(1))
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'DejaVu Sans'], 'font.size': 7.5,
                     'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.5,
                     'axes.linewidth': 0.6, 'savefig.dpi': 300, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'legend.frameon': False})
sto = json.load(open(os.path.join(DADOS, 'relaxacao.json')))['sistemas']
p321 = os.path.join(DADOS, 'relaxacao_321g.json'); s321 = json.load(open(p321))['sistemas'] if os.path.exists(p321) else {}
sistemas = list(sto.keys())
fig, ax = plt.subplots(1, 2, figsize=(3.35, 2.3), gridspec_kw=dict(width_ratios=[1.25, 1]))
marc = ['o', 's', 'D', '^']
CORSIS = [ST.AZUL, ST.VERMELHO, ST.VERDE, ST.LARANJA]   # uma cor por sistema
for i, rot in enumerate(sistemas):
    d0 = dmin_por_H(os.path.join(DADOS, rot + '_ideal.xyz')); d1 = dmin_por_H(os.path.join(DADOS, rot + '_relax_pbe_sto3g.xyz'))
    x = np.arange(1, len(d0) + 1)
    nome = rot.replace('Si', 'Si$_{').replace('H', '}$H$_{') + '}$'
    cs = CORSIS[i % len(CORSIS)]
    ax[0].plot(x, d0, ls='none', marker=marc[i], mfc='white', mec=cs, ms=3, mew=0.8, label=nome + ' ideal')
    ax[0].plot(x, d1, ls='none', marker=marc[i], mfc=cs, mec=cs, ms=2.6, label=nome + ' relaxed')
ax[0].axhline(2.0, color='0.5', lw=0.6, ls='--'); ax[0].text(1.5, 2.03, '2.0 Å', fontsize=6, color='0.3', va='bottom')
ax[0].axhline(1.423, color='0.7', lw=0.5, ls=':'); ax[0].axhline(2.417, color='0.7', lw=0.5, ls=':')
ax[0].set_xlabel('hydrogen (ascending order)'); ax[0].set_ylabel('shortest H$\\cdots$H distance (Å)')
ax[0].set_ylim(1.3, 4.0); ax[0].grid(True, color='0.9', lw=0.5); ax[0].set_axisbelow(True)
ax[0].legend(loc='center right', bbox_to_anchor=(1.02, 0.60), handletextpad=0.3, labelspacing=0.25)
ax[0].set_title('(a)', loc='left', fontsize=8, fontweight='bold')
# (b) delta gap
rotulos, dHF, dPBE, dHF3, dPBE3 = [], [], [], [], []
for rot in sistemas:
    rotulos.append(rot.replace('Si', 'Si$_{').replace('H', '}$H$_{') + '}$')
    dHF.append(sto[rot]['delta']['dgap_HF_eV']); dPBE.append(sto[rot]['delta']['dgap_PBE_eV'])
    dHF3.append(s321[rot]['delta']['dgap_HF_eV'] if rot in s321 else np.nan); dPBE3.append(s321[rot]['delta']['dgap_PBE_eV'] if rot in s321 else np.nan)
x = np.arange(len(sistemas)); w = 0.2
ax[1].bar(x - 1.5 * w, dHF, w, color=ST.METODO['HF'], edgecolor='black', lw=0.6, label='HF, STO-3G')
ax[1].bar(x - 0.5 * w, dPBE, w, color=ST.METODO['PBE'], edgecolor='black', lw=0.6, label='PBE, STO-3G')
if s321:
    ax[1].bar(x + 0.5 * w, dHF3, w, color='white', edgecolor=ST.METODO['HF'], hatch='////', lw=0.8, label='HF, 3-21G')
    ax[1].bar(x + 1.5 * w, dPBE3, w, color='white', edgecolor=ST.METODO['PBE'], hatch='....', lw=0.8, label='PBE, 3-21G')
ax[1].axhline(0, color='black', lw=0.6)
ax[1].set_xticks(x); ax[1].set_xticklabels(rotulos, fontsize=6.5); ax[1].set_ylabel('$E_g$(relaxed) $-$ $E_g$(ideal) (eV)')
ax[1].legend(loc='upper right', handlelength=1.4, handletextpad=0.4, labelspacing=0.25)
ax[1].grid(True, axis='y', color='0.9', lw=0.5); ax[1].set_axisbelow(True)
ax[1].set_title('(b)', loc='left', fontsize=8, fontweight='bold')
for a in ax:
    for sp in ('top', 'right'): a.spines[sp].set_visible(False)
fig.tight_layout(pad=0.3, w_pad=0.8)
out = os.path.join(RAIZ, ST.SAIDA, 'fig_relaxacao.png'); ST.salvar(fig, out); print('figura gravada em', out)
