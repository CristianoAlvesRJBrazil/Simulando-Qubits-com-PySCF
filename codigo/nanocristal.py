"""Gerador de nanocristais de Si passivados (Cap. 3 do livro)."""
import numpy as np

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
    return si, np.array(H)

def para_xyz(si, h, idx_p=None):
    L = []
    for i, p in enumerate(si):
        L.append("%-2s %.6f %.6f %.6f" % ("P" if i == idx_p else "Si", *p))
    for p in h:
        L.append("H %.6f %.6f %.6f" % tuple(p))
    return "\n".join(L)
