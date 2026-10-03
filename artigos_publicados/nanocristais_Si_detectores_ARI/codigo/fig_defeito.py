#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_defeito.py -- diagrama de niveis de Kohn-Sham (PBE/3-21G) do aglomerado perfeito e dos radicais
com uma ligacao pendente de superficie (um H removido), para Si10H16 e Si22H28.

Para cada sitio: bordas E_v (HOMO beta) e E_c (LUMO alfa) e os niveis da ligacao pendente -- ocupado
(alfa) e vazio (beta). Usa a geometria relaxada quando existe; sitios so com calculo vertical aparecem
com simbolos vazados. Entrada: dados/defeito_<sistema>.json. Saida: figuras/fig_defeito.{png,pdf}."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
AQUI = os.path.dirname(os.path.abspath(__file__)); RAIZ = os.path.join(AQUI, '..')
sys.path.insert(0, AQUI)
import estilo as ST
plt.rcParams.update(ST.RC)

SIS = [s for s in ('Si10H16', 'Si22H28') if os.path.exists(os.path.join(RAIZ, 'dados', 'defeito_%s.json' % s))]
dados = {s: json.load(open(os.path.join(RAIZ, 'dados', 'defeito_%s.json' % s))) for s in SIS}
ncat = [1 + len(dados[s].get('sitios', {})) for s in SIS]
fig, axs = plt.subplots(1, len(SIS), figsize=(7.0, 2.7), sharey=True, gridspec_kw=dict(width_ratios=ncat), squeeze=False)
W = 0.34
for ax, s in zip(axs[0], SIS):
    d = dados[s]; per = d['perfeito']['PBE']
    ax.hlines([per['HOMO_eV'], per['LUMO_eV']], -W, W, color=ST.CINZA, lw=2.2)
    ax.text(0, per['LUMO_eV'] + 0.12, '%.2f eV' % per['gap_eV'], ha='center', va='bottom', fontsize=6.3, color=ST.CINZA)
    rot = ['pristine']
    sitios = sorted(d.get('sitios', {}).items(), key=lambda t: int(t[0][6:]))
    for i, (ch, st) in enumerate(sitios, 1):
        geo = 'relaxado' if 'relaxado' in st else 'vertical'; r = st[geo]['PBE']
        cheio = geo == 'relaxado'
        ax.hlines([r['E_v_eV'], r['E_c_eV']], i - W, i + W, color=ST.CINZA, lw=2.2, alpha=1.0 if cheio else 0.55)
        ax.hlines(r['DB_ocupado_eV'], i - W, i + W, color=ST.VERMELHO, lw=1.6, alpha=1.0 if cheio else 0.6)
        ax.hlines(r['DB_vazio_eV'], i - W, i + W, color=ST.AZUL, lw=1.6, ls=(0, (2, 1)), alpha=1.0 if cheio else 0.6)
        ax.plot(i, r['DB_ocupado_eV'], marker='^', ms=4, color=ST.VERMELHO, mfc=ST.VERMELHO if cheio else 'white', zorder=4)
        ax.plot(i, r['DB_vazio_eV'], marker='v', ms=4, color=ST.AZUL, mfc='white', zorder=4)
        ax.text(i, r['E_c_eV'] + 0.12, '%.2f' % r['gap_efetivo_eV'], ha='center', va='bottom', fontsize=6.3, color=ST.VERMELHO)
        grupo = 'SiH' if st['tipo'] == 'SiH1' else 'SiH$_2$'
        rot.append('S%s\n%s%s' % (ch[6:], grupo, '' if cheio else '*'))
    ax.set_xticks(range(len(rot))); ax.set_xticklabels(rot, fontsize=6.3)
    ax.set_xlim(-0.6, len(rot) - 0.4)
    ax.set_title('(%s) %s' % ('ab'[SIS.index(s)], s.replace('Si', 'Si$_{').replace('H', '}$H$_{') + '}$'), loc='left', fontsize=8, fontweight='bold')
    ax.grid(True, axis='y', color='0.9', lw=0.5); ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
axs[0][0].set_ylabel('Kohn–Sham level, PBE/3-21G (eV)')
from matplotlib.lines import Line2D
leg = [Line2D([], [], color=ST.CINZA, lw=2.2, label='band edges ($E_v$, $E_c$)'),
       Line2D([], [], color=ST.VERMELHO, lw=1.6, marker='^', ms=4, label='dangling bond, occupied ($\\alpha$)'),
       Line2D([], [], color=ST.AZUL, lw=1.6, ls=(0, (2, 1)), marker='v', mfc='white', ms=4, label='dangling bond, empty ($\\beta$)')]
fig.legend(handles=leg, loc='lower center', ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=0.6)
out = os.path.join(RAIZ, ST.SAIDA, 'fig_defeito.png'); ST.salvar(fig, out); print('gravado', out)
