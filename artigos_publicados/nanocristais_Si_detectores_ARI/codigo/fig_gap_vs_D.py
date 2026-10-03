#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura: gap HOMO-LUMO vs diametro efetivo (HF e PBE, STO-3G; min_viz=2 vs 1),
com a equacao de Brus (parametros do Cap. 3) e o ajuste de Delerue-Allan-Lannoo
[PRB 48, 11024 (1993)] como referencias. Saida: ../figuras/fig_gap_vs_D.png (300 dpi, 3.35 in).
Legivel em preto-e-branco: identidade das series por marcador e estilo de traco."""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estilo as ST

AQUI = os.path.dirname(os.path.abspath(__file__))
ARQ_JSON = os.path.join(AQUI, '..', 'dados', 'gaps.json')
ARQ_FIG = os.path.join(AQUI, '..', ST.SAIDA, 'fig_gap_vs_D.png')
os.makedirs(os.path.dirname(ARQ_FIG), exist_ok=True)

HB2_2M, E2_4PE = 3.80998, 14.3996
EG_BULK, ME, MH, EPS = 1.12, 0.26, 0.38, 11.7
def gap_brus(d_nm):
    R = d_nm * 10 / 2
    return EG_BULK + HB2_2M * np.pi**2 / R**2 * (1/ME + 1/MH) - 1.786 * E2_4PE / (EPS * R)
def gap_delerue(d_nm):
    return 1.167 + 3.73 / d_nm**1.39

with open(ARQ_JSON) as f:
    db = json.load(f)

def serie(metodo, tag):
    """Pontos (D_ef, gap, formula) das series STO-3G por criterio (tag = 'mv1', 'mv2', 'mv2p').
    Prefere integrais exatas; aceita DF quando so ele existe (1,5 nm)."""
    pts = {}
    for k, r in db['calculos'].items():
        if not r.get('converged') or r['base'] != 'sto-3g' or r['metodo'] != metodo or k.split('|')[1] != tag:
            continue
        chave = (round(r['D_efetivo_nm'], 4), r['formula'])
        exato = r.get('modo_eri', '').startswith(('incore', 'direto'))
        if chave not in pts or exato:
            pts[chave] = r['gap_eV']
    pts = sorted((d, g, f) for (d, f), g in pts.items())
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts]), [p[2] for p in pts]

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7.5, 'axes.labelsize': 8, 'axes.titlesize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'lines.linewidth': 1.0, 'savefig.dpi': 300, 'figure.dpi': 300,
})

fig, ax = plt.subplots(figsize=(3.35, 3.3))
fig.patch.set_facecolor('white'); ax.set_facecolor('white')

# cor distingue o metodo; marcador e traco distinguem o criterio (leitura em P&B preservada)
_HF, _PBE = ST.METODO['HF'], ST.METODO['PBE']
estilos = {
    ('HF', 'mv2'):   dict(marker='o', mfc=_HF,    ls='-',  color=_HF,  label='HF, $n_{\\min}=2$'),
    ('HF', 'mv1'):   dict(marker='s', mfc='white', ls='--', color=_HF,  label='HF, $n_{\\min}\\leq1$'),
    ('HF', 'mv2p'):  dict(marker='D', mfc='white', ls='-.', color=_HF,  label='HF, $n_{\\min}=2$ + filling'),
    ('PBE', 'mv2'):  dict(marker='o', mfc=_PBE,   ls='-',  color=_PBE, label='PBE, $n_{\\min}=2$'),
    ('PBE', 'mv1'):  dict(marker='s', mfc='white', ls='--', color=_PBE, label='PBE, $n_{\\min}\\leq1$'),
    ('PBE', 'mv2p'): dict(marker='D', mfc='white', ls='-.', color=_PBE, label='PBE, $n_{\\min}=2$ + filling'),
}
xmax = 0.0
for (met, mv), st in estilos.items():
    x, y, fm = serie(met, mv)
    if len(x) == 0:
        continue
    xmax = max(xmax, x.max())
    ax.plot(x, y, marker=st['marker'], mfc=st['mfc'], mec=st['color'], ms=4.0, mew=0.8,
            ls=st['ls'], color=st['color'], lw=1.0, label=st['label'], zorder=3)

dd = np.linspace(0.6, 1.7, 300)
ax.plot(dd, gap_brus(dd), ls=':', color=ST.CINZA, lw=1.2, label='Brus, effective mass (Eq. 1)', zorder=2)
ax.plot(dd, gap_delerue(dd), ls=(0, (5, 1, 1, 1)), color=ST.ROXO, lw=1.2,
        label='Delerue et al. (1993) fit (Eq. 2)', zorder=2)
ax.axhline(1.12, color='0.6', lw=0.5, ls='-', zorder=1)
ax.text(1.68, 1.12 + 0.3, 'bulk Si (1.12 eV)', ha='right', va='bottom', fontsize=5.5, color='0.3')

ax.set_xlim(0.6, 1.7); ax.set_ylim(0, 20)
ax.set_xlabel('effective diameter $D_{\\mathrm{eff}}$ (nm)')
ax.set_ylabel('HOMO$-$LUMO gap (eV)')
ax.grid(True, color='0.9', lw=0.5); ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False, handlelength=2.4, borderaxespad=0.0, labelspacing=0.3, columnspacing=1.0, fontsize=6)
fig.tight_layout(pad=0.3)
ST.salvar(fig, ARQ_FIG)
print('figura gravada em', ARQ_FIG)
