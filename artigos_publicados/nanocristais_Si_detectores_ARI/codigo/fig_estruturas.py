#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figura: projecao ortografica dos nanocristais de D = 1,5 nm com min_viz = 2 (Si72H64) e
com preenchimento (Si84H64), destacando os pares H...H < 2,0 A. Saida: ../figuras/fig_estruturas.png
(300 dpi, 3.35 in, monocromatica)."""
import os, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import estilo as ST
AQUI = os.path.dirname(os.path.abspath(__file__)); RAIZ = os.path.join(AQUI, '..')
def ler_xyz(p):
    L = open(p).read().split('\n'); n = int(L[0]); at = []
    for l in L[2:2+n]:
        s, x, y, z = l.split()[:4]; at.append((s, float(x), float(y), float(z)))
    sim = np.array([a[0] for a in at]); xyz = np.array([a[1:] for a in at]); return sim, xyz
def rot(xyz, ang_deg, eixo):
    a = np.radians(ang_deg); c, s = np.cos(a), np.sin(a)
    R = {'x': [[1,0,0],[0,c,-s],[0,s,c]], 'y': [[c,0,s],[0,1,0],[-s,0,c]], 'z': [[c,-s,0],[s,c,0],[0,0,1]]}[eixo]
    return xyz @ np.array(R).T
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans','DejaVu Sans'], 'font.size': 7.5,
                     'savefig.dpi': 300, 'figure.facecolor': 'white', 'axes.facecolor': 'white'})
arqs = [(os.path.join(RAIZ,'dados','nanocristal_D1.5nm_minviz2.xyz'), '(a) $n_{\\min}=2$: Si$_{72}$H$_{64}$'),
        (os.path.join(RAIZ,'dados','xyz','nanocristal_D1.5nm_minviz2_preencher.xyz'), '(b) $n_{\\min}=2$ + filling: Si$_{84}$H$_{64}$')]
fig, axs = plt.subplots(1, 2, figsize=(3.35, 2.15))
for ax, (arq, tit) in zip(axs, arqs):
    sim, xyz = ler_xyz(arq); xyz = xyz - xyz[sim=='Si'].mean(0)
    xyz = rot(rot(xyz, 20, 'x'), 30, 'y')          # vista generica (evita superposicao total)
    si = sim == 'Si'; h = sim == 'H'
    # ligacoes Si-Si
    P = xyz[si]; d = np.linalg.norm(P[:,None]-P[None], axis=2)
    for i, j in zip(*np.where((d > 0.1) & (d < 2.6))):
        if i < j: ax.plot(P[[i,j],0], P[[i,j],1], color=ST.CINZACLARO, lw=0.7, zorder=1)
    # ligacoes Si-H
    Q = xyz[h]; dsh = np.linalg.norm(P[:,None]-Q[None], axis=2)
    for i, j in zip(*np.where(dsh < 1.6)):
        ax.plot([P[i,0],Q[j,0]],[P[i,1],Q[j,1]], color='0.75', lw=0.4, zorder=1)
    # atomos (mais profundos primeiro)
    ordem = np.argsort(xyz[:,2])
    for k in ordem:
        if sim[k] == 'Si': ax.plot(xyz[k,0], xyz[k,1], 'o', ms=3.2, mfc=ST.AZUL, mec='black', mew=0.3, zorder=2)
        else: ax.plot(xyz[k,0], xyz[k,1], 'o', ms=1.6, mfc='white', mec='0.4', mew=0.4, zorder=3)
    # pares H...H < 2,0 A
    dhh = np.linalg.norm(Q[:,None]-Q[None], axis=2); npar = 0
    for i, j in zip(*np.where((dhh > 0.1) & (dhh < 2.0))):
        if i < j:
            ax.plot(Q[[i,j],0], Q[[i,j],1], color=ST.VERMELHO, lw=1.8, zorder=4)
            ax.plot(Q[[i,j],0].mean(), Q[[i,j],1].mean(), 's', ms=3.5, mfc='none', mec=ST.VERMELHO, mew=1.0, zorder=5); npar += 1
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(tit + ('\n%d H$\\cdots$H pairs < 2.0 Å' % npar if npar else '\nno H$\\cdots$H pair < 2.0 Å'), fontsize=6.8, pad=3)
fig.tight_layout(pad=0.3, rect=[0, 0, 1, 0.97])
out = os.path.join(RAIZ, ST.SAIDA, 'fig_estruturas.png'); ST.salvar(fig, out); print('gravado', out)
