#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura: o efeito do criterio de corte sobrevive a troca de base?
Delta(gap) = gap(min_viz<=1) - gap(min_viz=2) no mesmo diametro nominal, em STO-3G e 3-21G,
comparado ao previsto apenas pela diferenca de tamanho (Eq. de Delerue nos dois D_efetivos).
Saida: ../figuras/fig_robustez.png (300 dpi, 3,35 in, P&B)."""
import os, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estilo as ST

AQUI = os.path.dirname(os.path.abspath(__file__)); RAIZ = os.path.join(AQUI, '..')
db = json.load(open(os.path.join(RAIZ, 'dados', 'gaps.json')))['calculos']

def busca(D, tag, base, met):
    c = [v for k, v in db.items() if v.get('converged') and abs(v['D_nominal_nm'] - D) < 1e-6
         and v['base'] == base and v['metodo'] == met and k.split('|')[1] == tag]
    c.sort(key=lambda v: 'DF' in v.get('modo_eri', ''))
    return c[0] if c else None

DIAMS = [0.8, 1.0, 1.2]
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'DejaVu Sans'],
                     'font.size': 7.5, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'savefig.dpi': 300,
                     'figure.facecolor': 'white', 'axes.facecolor': 'white', 'legend.frameon': False})
fig, ax = plt.subplots(figsize=(3.35, 2.5))
# cor distingue o metodo; preenchimento e traco distinguem a base
_HF, _PBE = ST.METODO['HF'], ST.METODO['PBE']
series = [('sto-3g', 'HF',  dict(marker='o', mfc=_HF,    ls='-',  color=_HF,  label='HF, STO-3G')),
          ('sto-3g', 'PBE', dict(marker='s', mfc=_PBE,   ls='-',  color=_PBE, label='PBE, STO-3G')),
          ('3-21g',  'HF',  dict(marker='o', mfc='white', ls='--', color=_HF,  label='HF, 3-21G')),
          ('3-21g',  'PBE', dict(marker='s', mfc='white', ls='--', color=_PBE, label='PBE, 3-21G'))]
for base, met, st in series:
    x, y = [], []
    for D in DIAMS:
        a, b = busca(D, 'mv1', base, met), busca(D, 'mv2', base, met)
        if a and b:
            x.append(D); y.append(a['gap_eV'] - b['gap_eV'])
    if x:
        ax.plot(x, y, ms=4, mew=0.8, mec=st['color'], lw=1.0, **{k: v for k, v in st.items() if k != 'color'} | {'color': st['color']})
# previsto so pelo tamanho
xs, ys = [], []
for D in DIAMS:
    a, b = busca(D, 'mv1', 'sto-3g', 'HF'), busca(D, 'mv2', 'sto-3g', 'HF')
    if a and b:
        xs.append(D); ys.append(a['gap_delerue_Defetivo_eV'] - b['gap_delerue_Defetivo_eV'])
ax.plot(xs, ys, ls=':', color=ST.CINZA, lw=1.4, marker='^', mfc=ST.CINZA, ms=3.5, label='size-only estimate (Eq. 2)')
ax.axhline(0, color='black', lw=0.5)
ax.set_xlabel('nominal diameter $D$ (nm)')
ax.set_ylabel('$\\Delta E_g = E_g(n_{\\min}\\leq1) - E_g(n_{\\min}=2)$ (eV)')
ax.set_xticks(DIAMS); ax.grid(True, color='0.9', lw=0.5); ax.set_axisbelow(True)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
ax.legend(loc='lower right', handlelength=2.4, labelspacing=0.3)
fig.tight_layout(pad=0.3)
out = os.path.join(RAIZ, ST.SAIDA, 'fig_robustez.png'); ST.salvar(fig, out)
print('figura gravada em', out)
