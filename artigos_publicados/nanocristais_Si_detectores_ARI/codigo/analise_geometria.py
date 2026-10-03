#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analise_geometria.py -- Resultados R1 e R2 do artigo (Secoes 4.1 e 4.2).

Analise puramente geometrica dos nanocristais de Si passivados com H gerados
pelo algoritmo do Cap. 3 do livro "Simulando Qubits com PySCF":

    rede diamante -> corte esferico de raio R -> remocao iterativa de Si com
    menos de `min_viz` vizinhos Si -> passivacao com H (1,48 A) na direcao de
    cada vizinho Si ausente.

R1 -- sensibilidade ao criterio de coordenacao minima (min_viz).
R2 -- contatos H-H (origem do artefato de ~1,4 A registrado no livro).

Script AUTOCONTIDO: as funcoes `rede_diamante`, `construir` e `para_xyz` sao
copias fieis de codigo/nanocristal.py (nao importa nada do livro). Nao ha
nenhuma etapa aleatoria; a semente e fixada apenas por disciplina de
reprodutibilidade. Custo: segundos de CPU em 1 thread.

Uso:  OMP_NUM_THREADS=1 python3 analise_geometria.py
Saidas (relativas a pasta do artigo, um nivel acima de codigo/):
    dados/geometria_sensibilidade.csv, dados/geometria_sensibilidade.json,
    dados/nanocristal_D1.5nm_minviz{1,2}.xyz, dados/xyz/*.xyz,
    figuras/fig_contatos_HH.png, figuras/fig_distribuicao_HH.png,
    figuras/fig_estequiometria.png
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys, json, csv, time, math, itertools, platform
from collections import Counter, defaultdict
import numpy as np

np.random.seed(20260906)          # nenhuma etapa aleatoria; semente por disciplina

# ---------------------------------------------------------------------------
# 1. Gerador (copia fiel de codigo/nanocristal.py)
# ---------------------------------------------------------------------------
A_SI, D_SIH, CORTE = 5.431, 1.48, 2.6          # A: parametro de rede, Si-H, corte de vizinhanca
D_SISI = A_SI * math.sqrt(3) / 4               # 2,3517 A (ligacao Si-Si ideal)

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

# ---------------------------------------------------------------------------
# 2. Versao instrumentada (mesmo algoritmo, com historico) + criterio extra
# ---------------------------------------------------------------------------
def vizinhos_rede(rede):
    """Lista de vizinhos (indices) de cada sitio da rede completa, d < CORTE."""
    from scipy.spatial import cKDTree
    arv = cKDTree(rede)
    pares = arv.query_pairs(CORTE, output_type="ndarray")
    viz = [[] for _ in range(len(rede))]
    for i, j in pares:
        viz[i].append(j); viz[j].append(i)
    return [np.array(v, dtype=int) for v in viz]

def construir_info(diametro_nm, n_celulas=8, min_viz=2, preencher=False):
    """Igual a `construir`, mas devolve tambem o historico da filtragem.
    preencher=True ativa o criterio extra: sitios EXTERNOS com >= 2 vizinhos
    internos sao incorporados ao aglomerado (elimina, por construcao, pares de
    H apontando para o mesmo sitio vago)."""
    raio = diametro_nm * 10 / 2
    rede = rede_diamante(n_celulas)
    rede = rede - rede.mean(axis=0)
    viz = vizinhos_rede(rede)
    dentro = np.linalg.norm(rede, axis=1) <= raio
    n_corte = int(dentro.sum())
    hist_rem, hist_add, it = [], [], 0
    while True:
        it += 1
        nviz = np.array([np.sum(dentro[viz[i]]) if dentro[i] else -1 for i in range(len(rede))])
        ruins = np.where(dentro & (nviz < min_viz))[0]
        hist_rem.append(len(ruins))
        dentro[ruins] = False
        n_add = 0
        if preencher:
            nviz_ext = np.array([np.sum(dentro[viz[i]]) if not dentro[i] else -1 for i in range(len(rede))])
            novos = np.where((~dentro) & (nviz_ext >= 2))[0]
            n_add = len(novos)
            dentro[novos] = True
        hist_add.append(n_add)
        if len(ruins) == 0 and n_add == 0:
            break
        if it > 500:
            raise RuntimeError("nao convergiu")
    si = rede[dentro]
    H, pai, alvo = [], [], []
    for i in np.where(dentro)[0]:
        p = rede[i]
        for j in viz[i]:
            if not dentro[j]:
                u = (rede[j] - p) / np.linalg.norm(rede[j] - p)
                H.append(p + D_SIH * u); pai.append(i); alvo.append(j)
    H = np.array(H).reshape(-1, 3)
    info = dict(n_corte=n_corte, iteracoes=it, removidos_por_iteracao=hist_rem,
                adicionados_por_iteracao=hist_add,
                n_removidos=int(sum(hist_rem)), n_adicionados=int(sum(hist_add)))
    return si, H, info, np.array(pai, dtype=int), np.array(alvo, dtype=int), rede

# ---------------------------------------------------------------------------
# 3. Analise geometrica
# ---------------------------------------------------------------------------
def analisar(si, H, pai, alvo, rede):
    r = {}
    nSi, nH = len(si), len(H)
    r.update(n_Si=nSi, n_H=nH, n_atomos=nSi + nH, n_eletrons=14 * nSi + nH,
             formula=("Si%dH%d" % (nSi, nH)) if nSi else "(vazio)")
    if nSi == 0:
        return r
    r["razao_H_Si"] = nH / nSi
    # coordenacao
    dSS = np.linalg.norm(si[:, None] - si[None, :], axis=2)
    nvS = ((dSS > 0.1) & (dSS < CORTE)).sum(1)
    if nH:
        dSH = np.linalg.norm(si[:, None] - H[None, :], axis=2)
        nvH = (dSH < 1.6).sum(1)
    else:
        dSH = np.zeros((nSi, 0)); nvH = np.zeros(nSi, int)
    coord = nvS + nvH
    r.update(coord_min=int(coord.min()), coord_max=int(coord.max()),
             todos_tetracoordenados=bool(np.all(coord == 4)),
             min_SiSi=float(dSS[dSS > 0.1].min()) if nSi > 1 else None,
             min_SiH=float(dSH.min()) if nH else None)
    grupos = Counter(nvH.tolist())
    r.update(n_Si_nucleo=grupos.get(0, 0), n_SiH=grupos.get(1, 0), n_SiH2=grupos.get(2, 0),
             n_SiH3=grupos.get(3, 0), n_SiH4=grupos.get(4, 0))
    r["n_Si_superficie"] = int((nvH >= 1).sum())
    r["fracao_Si_superficie"] = r["n_Si_superficie"] / nSi
    r["vizinhos_Si_min"] = int(nvS.min()); r["vizinhos_Si_max"] = int(nvS.max())
    # atomo extremo em x (argumento do colapso para min_viz >= 3)
    ix = int(np.argmax(si[:, 0])); r["nviz_Si_do_atomo_extremo_x"] = int(nvS[ix])
    # raios
    rS = np.linalg.norm(si, axis=1)
    r.update(R_Si_max=float(rS.max()), D_ef_Si_nm=float(2 * rS.max() / 10),
             R_H_max=float(np.linalg.norm(H, axis=1).max()) if nH else None,
             R_Si_medio=float(rS.mean()))
    dens = 8 / A_SI**3                                    # Si/A^3 (diamante)
    r["D_eq_contagem_nm"] = float(2 * (3 * nSi / (4 * math.pi * dens))**(1 / 3) / 10)
    r["R_giro_Si"] = float(math.sqrt((rS**2).mean()))
    # ---- contatos H-H
    if nH < 2:
        r.update(min_HH=None); return r
    dHH = np.linalg.norm(H[:, None] - H[None, :], axis=2)
    np.fill_diagonal(dHH, np.inf)
    dmin = dHH.min(1)                                     # menor H-H de cada H
    r["min_HH"] = float(dmin.min())
    r["dmin_HH_por_H"] = [float(x) for x in np.round(dmin, 4)]
    # classificacao de cada par (i<j) com d < 4.0 A
    iu, ju = np.triu_indices(nH, 1)
    dp = dHH[iu, ju]
    sel = dp < 4.0
    classes, dists = [], []
    for i, j, d in zip(iu[sel], ju[sel], dp[sel]):
        pi, pj = pai[i], pai[j]
        if pi == pj:
            c = "geminal (mesmo Si)"
        elif alvo[i] == alvo[j]:
            c = "convergente (mesmo sitio vago; Si 2os vizinhos)"
        else:
            dS = np.linalg.norm(rede[pi] - rede[pj])
            if dS < CORTE:
                c = "vicinal (Si-Si ligados)"
            elif dS < 4.0:
                c = "2os vizinhos, sitios vagos distintos"
            elif dS < 4.7:
                c = "3os vizinhos (sitios vagos adjacentes)"
            else:
                c = "outros"
        classes.append(c); dists.append(float(d))
    classes = np.array(classes); dists = np.array(dists)
    r["pares_HH_por_limiar"] = {}
    for lim in (1.6, 2.0, 2.2, 2.4, 2.5):
        m = dists < lim
        r["pares_HH_por_limiar"][str(lim)] = dict(
            n_pares=int(m.sum()),
            n_H_envolvidos=int(len(set(iu[sel][m]) | set(ju[sel][m]))),
            fracao_H_envolvidos=float(len(set(iu[sel][m]) | set(ju[sel][m])) / nH),
            classes=dict(Counter(classes[m].tolist())))
    # espectro de distancias discretas (< 4 A)
    esp = Counter(zip(np.round(dists, 3).tolist(), classes.tolist()))
    r["espectro_HH"] = [dict(d=float(k[0]), classe=k[1], n=int(v))
                        for k, v in sorted(esp.items())]
    # multiplicidade dos sitios vagos (quantos H apontam para o mesmo sitio)
    mult = Counter(Counter(alvo.tolist()).values())
    r["sitios_vagos_por_multiplicidade"] = {str(k): int(v) for k, v in sorted(mult.items())}
    r["n_pares_convergentes_teorico"] = int(sum(v * k * (k - 1) // 2 for k, v in mult.items()))
    return r

# ---------------------------------------------------------------------------
# 4. Execucao
# ---------------------------------------------------------------------------
def main():
    base = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    ddir = os.path.join(base, "dados"); fdir = os.path.join(base, os.environ.get("FIG_DIR", "figuras"))
    os.makedirs(os.path.join(ddir, "xyz"), exist_ok=True); os.makedirs(fdir, exist_ok=True)

    DIAMS = [0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
    MINVIZ = [0, 1, 2, 3, 4]
    t_ini = time.time()

    # --- distancias teoricas (rede ideal)
    ang = math.degrees(math.acos(-1 / 3))                  # 109,471 graus
    d_gem = 2 * D_SIH * math.sin(math.radians(ang / 2))    # = D_SIH*sqrt(8/3)
    d_conv = 2 * (D_SISI - D_SIH) * math.sin(math.radians(ang / 2))
    teoria = dict(angulo_tetraedrico_graus=ang, Si_Si_ideal=D_SISI, Si_H=D_SIH,
                  HH_geminal=d_gem, HH_convergente_mesmo_sitio=d_conv,
                  HH_convergente_formula="2*(d_SiSi - d_SiH)*sin(109.47/2)",
                  distancia_H_ao_sitio_vago=D_SISI - D_SIH,
                  segundos_vizinhos_SiSi=A_SI / math.sqrt(2),
                  terceiros_vizinhos_SiSi=A_SI * math.sqrt(11) / 4,
                  vdW_HH_Bondi=2 * 1.20, vdW_HH_RowlandTaylor=2 * 1.10, H2_ligacao=0.741)
    print("== Distancias teoricas (rede ideal, a = %.3f A) ==" % A_SI)
    for k, v in teoria.items():
        print("  %-32s %s" % (k, ("%.4f" % v) if isinstance(v, float) else v))

    registros, tempos = [], {}
    xyz_guardados = {}
    for D in DIAMS:
        for mv in MINVIZ:
            t0 = time.time()
            si, H, info, pai, alvo, rede = construir_info(D, min_viz=mv)
            # verificacao: a versao instrumentada reproduz a original
            if mv >= 1:
                si0, H0 = construir(D, min_viz=mv)
                assert len(si0) == len(si) and (len(si) == 0 or np.allclose(np.sort(si0, 0), np.sort(si, 0)))
                assert len(H0) == len(H) and (len(H) == 0 or np.allclose(np.sort(H0, 0), np.sort(H, 0)))
            r = analisar(si, H, pai, alvo, rede)
            r.update(D_nm=D, min_viz=mv, criterio="min_viz=%d" % mv, tempo_s=time.time() - t0, **info)
            registros.append(r)
            if len(si):
                nome = "nanocristal_D%.1fnm_minviz%d.xyz" % (D, mv)
                txt = "%d\nSi nanocristal D=%.1f nm min_viz=%d %s (rede ideal, a=%.3f A, Si-H=%.2f A)\n%s\n" % (
                    len(si) + len(H), D, mv, r["formula"], A_SI, D_SIH, para_xyz(si, H))
                open(os.path.join(ddir, "xyz", nome), "w").write(txt)
                xyz_guardados[(D, mv)] = (txt, r)
        # criterio extra: min_viz=2 + preenchimento de sitios com >=2 vizinhos internos
        t0 = time.time()
        si, H, info, pai, alvo, rede = construir_info(D, min_viz=2, preencher=True)
        r = analisar(si, H, pai, alvo, rede)
        r.update(D_nm=D, min_viz=2, criterio="min_viz=2+preencher", tempo_s=time.time() - t0, **info)
        registros.append(r)
        if len(si):
            nome = "nanocristal_D%.1fnm_minviz2_preencher.xyz" % D
            open(os.path.join(ddir, "xyz", nome), "w").write(
                "%d\nSi nanocristal D=%.1f nm min_viz=2+preencher %s\n%s\n" % (len(si) + len(H), D, r["formula"], para_xyz(si, H)))

    # XYZ pedidos explicitamente (D=1,5 nm): min_viz=2 e o vizinho sobrevivente (min_viz=1);
    # min_viz=3 nao produz estrutura (colapsa) -- registrado no relatorio.
    for mv in (1, 2):
        open(os.path.join(ddir, "nanocristal_D1.5nm_minviz%d.xyz" % mv), "w").write(xyz_guardados[(1.5, mv)][0])
    open(os.path.join(ddir, "nanocristal_D1.5nm_minviz3.xyz"), "w").write(
        "0\nmin_viz=3: remocao iterativa colapsa o aglomerado (nenhum atomo sobrevive) -- ver geometria_sensibilidade.json\n")

    # ---------------- Tabelas no terminal ----------------
    def fmt(x, f="%.3f"):
        return "-" if x is None else (f % x if isinstance(x, float) else str(x))
    print("\n== R1: sensibilidade ao criterio min_viz ==")
    cab = ("D(nm) crit         formula     N_at N_el  H/Si  fS   SiH SiH2 SiH3 SiH4 nucleo "
           "R_Si_max D_ef  D_eq  cmin cmax minSiSi minSiH it rem add")
    print(cab)
    for r in registros:
        print("%4.1f  %-17s %-11s %4d %5s %5s %5s %4s %4s %4s %4s %5s   %6s %5s %5s  %4s %4s %6s %6s %2d %3d %3d" % (
            r["D_nm"], r["criterio"], r["formula"], r["n_atomos"], r["n_eletrons"],
            fmt(r.get("razao_H_Si"), "%.3f"), fmt(r.get("fracao_Si_superficie"), "%.3f"),
            fmt(r.get("n_SiH")), fmt(r.get("n_SiH2")), fmt(r.get("n_SiH3")), fmt(r.get("n_SiH4")), fmt(r.get("n_Si_nucleo")),
            fmt(r.get("R_Si_max"), "%.3f"), fmt(r.get("D_ef_Si_nm"), "%.3f"), fmt(r.get("D_eq_contagem_nm"), "%.3f"),
            fmt(r.get("coord_min")), fmt(r.get("coord_max")), fmt(r.get("min_SiSi"), "%.3f"), fmt(r.get("min_SiH"), "%.3f"),
            r["iteracoes"], r["n_removidos"], r["n_adicionados"]))
    print("\n== R2: contatos H-H ==")
    print("D(nm) crit         formula     N_H  minHH  n<1.6 n<2.0 n<2.2 n<2.4 n<2.5  H<2.0 fH<2.0  sitios(k=2) sitios(k=3)")
    for r in registros:
        if r.get("min_HH") is None:
            print("%4.1f  %-17s %-11s %4d   -" % (r["D_nm"], r["criterio"], r["formula"], r["n_H"])); continue
        p = r["pares_HH_por_limiar"]; m = r["sitios_vagos_por_multiplicidade"]
        print("%4.1f  %-17s %-11s %4d  %.3f  %5d %5d %5d %5d %5d  %5d  %.3f   %5s %5s" % (
            r["D_nm"], r["criterio"], r["formula"], r["n_H"], r["min_HH"],
            p["1.6"]["n_pares"], p["2.0"]["n_pares"], p["2.2"]["n_pares"], p["2.4"]["n_pares"], p["2.5"]["n_pares"],
            p["2.0"]["n_H_envolvidos"], p["2.0"]["fracao_H_envolvidos"], m.get("2", 0), m.get("3", 0)))
    print("\n== R2: espectro de distancias H-H (< 4 A), D = 1,5 nm ==")
    for r in registros:
        if r["D_nm"] == 1.5 and r.get("espectro_HH"):
            print("  [%s] %s" % (r["criterio"], r["formula"]))
            for e in r["espectro_HH"]:
                print("     d = %.3f A  n = %4d   %s" % (e["d"], e["n"], e["classe"]))
    print("\n== Colapso para min_viz >= 3 (removidos por iteracao) ==")
    for r in registros:
        if r["min_viz"] >= 3 and r["criterio"].startswith("min_viz=%d" % r["min_viz"]):
            print("  D=%.1f min_viz=%d: corte=%d Si -> %s -> sobram %d" % (
                r["D_nm"], r["min_viz"], r["n_corte"], r["removidos_por_iteracao"], r["n_Si"]))
    print("\n  Argumento geral: na rede diamante cada Si tem 2 vizinhos com +x e 2 com -x;")
    print("  o Si de maior x de QUALQUER aglomerado finito tem <= 2 vizinhos Si e e removido.")
    print("  nviz do Si extremo em x (min_viz=2):", [(r["D_nm"], r["nviz_Si_do_atomo_extremo_x"]) for r in registros if r["criterio"] == "min_viz=2"])

    # ---------------- CSV / JSON ----------------
    cols = ["D_nm", "min_viz", "criterio", "formula", "n_Si", "n_H", "n_atomos", "n_eletrons", "razao_H_Si",
            "n_Si_superficie", "fracao_Si_superficie", "n_Si_nucleo", "n_SiH", "n_SiH2", "n_SiH3", "n_SiH4",
            "vizinhos_Si_min", "coord_min", "coord_max", "todos_tetracoordenados", "min_SiSi", "min_SiH",
            "R_Si_max", "D_ef_Si_nm", "R_H_max", "R_giro_Si", "D_eq_contagem_nm",
            "min_HH", "n_pares_HH_lt1.6", "n_pares_HH_lt2.0", "n_pares_HH_lt2.2", "n_pares_HH_lt2.4", "n_pares_HH_lt2.5",
            "n_H_lt2.0", "fracao_H_lt2.0", "sitios_k2", "sitios_k3", "n_corte", "iteracoes", "n_removidos", "n_adicionados", "tempo_s"]
    with open(os.path.join(ddir, "geometria_sensibilidade.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for r in registros:
            p = r.get("pares_HH_por_limiar", {}); m = r.get("sitios_vagos_por_multiplicidade", {})
            linha = dict(r)
            for lim in ("1.6", "2.0", "2.2", "2.4", "2.5"):
                linha["n_pares_HH_lt" + lim] = p.get(lim, {}).get("n_pares")
            linha["n_H_lt2.0"] = p.get("2.0", {}).get("n_H_envolvidos")
            linha["fracao_H_lt2.0"] = p.get("2.0", {}).get("fracao_H_envolvidos")
            linha["sitios_k2"] = m.get("2", 0 if r["n_Si"] else None); linha["sitios_k3"] = m.get("3", 0 if r["n_Si"] else None)
            w.writerow([("" if linha.get(c) is None else (("%.6g" % linha[c]) if isinstance(linha[c], float) else linha[c])) for c in cols])
    meta = dict(script=os.path.abspath(__file__), data="2026-09-06", python=platform.python_version(),
                numpy=np.__version__, parametros=dict(A_SI=A_SI, D_SIH=D_SIH, CORTE=CORTE, n_celulas=8,
                diametros_nm=DIAMS, min_viz=MINVIZ, centro="media da rede completa = sitio intersticial tetraedrico (7/8,7/8,7/8)a"),
                distancias_teoricas=teoria, tempo_total_s=time.time() - t_ini, threads=os.environ["OMP_NUM_THREADS"])
    json.dump(dict(meta=meta, registros=registros), open(os.path.join(ddir, "geometria_sensibilidade.json"), "w"),
              indent=1, ensure_ascii=False)

    # ---------------- Figuras (P&B, 8,5 cm, 300 dpi) ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import estilo as ST
    from matplotlib import font_manager
    fam = "Liberation Sans" if any("Liberation Sans" == f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": fam, "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "lines.linewidth": 1.0, "lines.markersize": 4, "savefig.dpi": 300,
                         "figure.facecolor": "white", "axes.facecolor": "white", "legend.frameon": False})
    # min_viz = 0 e 1 produzem estruturas identicas (nenhum Si isolado sobrevive ao corte esferico);
    # por isso uma unica serie "min_viz <= 1" representa o corte sem filtro efetivo.
    estilo = {1: dict(**ST.CRITERIO[1], label="$n_{\\min}\\leq1$ (no effective filter)"),
              2: dict(**ST.CRITERIO[2], label="$n_{\\min}=2$")}
    EST_P = ST.CRITERIO['2p']
    import re
    def formula_mpl(fm):   # Si72H64 -> Si$_{72}$H$_{64}$
        return re.sub(r"(\d+)", r"$_{\1}$", fm)
    LARG = 3.35
    def serie(mv, chave, crit=None):
        crit = crit or ("min_viz=%d" % mv)
        rs = [r for r in registros if r["criterio"] == crit]
        return [r["D_nm"] for r in rs], [chave(r) for r in rs]

    # Fig 1: pares H-H < 2,0 A vs D
    fig, ax = plt.subplots(2, 1, figsize=(LARG, 4.3), sharex=True)
    for mv in (1, 2):
        x, y = serie(mv, lambda r: r["pares_HH_por_limiar"]["2.0"]["n_pares"]); ax[0].plot(x, y, **estilo[mv])
        x, y = serie(mv, lambda r: 100 * r["pares_HH_por_limiar"]["2.0"]["fracao_H_envolvidos"]); ax[1].plot(x, y, **estilo[mv])
    x, y = serie(2, lambda r: r["pares_HH_por_limiar"]["2.0"]["n_pares"], "min_viz=2+preencher")
    ax[0].plot(x, y, **EST_P, label="$n_{\\min}=2$ + filling")
    x, y = serie(2, lambda r: 100 * r["pares_HH_por_limiar"]["2.0"]["fracao_H_envolvidos"], "min_viz=2+preencher")
    ax[1].plot(x, y, **EST_P)
    ax[0].set_ylabel("H$\\cdots$H pairs with $d$ < 2.0 Å")
    ax[1].set_ylabel("H atoms involved (%)"); ax[1].set_xlabel("nominal diameter $D$ (nm)")
    ax[0].legend(loc="upper left")
    ax[0].text(-0.22, 1.0, "(a)", transform=ax[0].transAxes, va="bottom", ha="left", fontsize=8, fontweight="bold")
    ax[1].text(-0.22, 1.0, "(b)", transform=ax[1].transAxes, va="bottom", ha="left", fontsize=8, fontweight="bold")
    for a in ax:
        a.grid(True, color="0.85", lw=0.5); a.set_axisbelow(True); a.set_xticks(DIAMS)
    ax[1].set_ylim(bottom=-2)
    # anota valores da serie principal (min_viz=2)
    for xx, yy in zip(*serie(2, lambda r: r["pares_HH_por_limiar"]["2.0"]["n_pares"])):
        ax[0].annotate("%d" % yy, (xx, yy), textcoords="offset points", xytext=(0, 4), ha="center", fontsize=6.5)
    fig.tight_layout(h_pad=0.4); ST.salvar(fig, os.path.join(fdir, "fig_contatos_HH.png")); plt.close(fig)

    # Fig 2: distribuicao acumulada da menor distancia H-H por H, D = 1,5 nm
    fig, ax = plt.subplots(figsize=(LARG, 2.6))
    for mv in (1, 2):
        r = [r for r in registros if r["criterio"] == "min_viz=%d" % mv and r["D_nm"] == 1.5][0]
        d = np.sort(r["dmin_HH_por_H"]); y = np.arange(1, len(d) + 1) / len(d) * 100
        d2 = np.concatenate([[1.0], d, [5.0]]); y2 = np.concatenate([[0], y, [100]])
        st = dict(estilo[mv]); st.pop("marker"); st.pop("mfc", None)
        ax.step(d2, y2, where="post", label="%s, %s" % (st.pop("label"), formula_mpl(r["formula"])), **st)
    r = [r for r in registros if r["criterio"] == "min_viz=2+preencher" and r["D_nm"] == 1.5][0]
    d = np.sort(r["dmin_HH_por_H"]); y = np.arange(1, len(d) + 1) / len(d) * 100
    ax.step(np.concatenate([[1.0], d, [5.0]]), np.concatenate([[0], y, [100]]), where="post", color=ST.VERDE, ls="-.",
            label="$n_{\\min}=2$ + filling, %s" % formula_mpl(r["formula"]))
    for xv, txt in ((2.0, "2.0 Å"), (2.4, "vdW 2.4 Å")):
        ax.axvline(xv, color="0.4", lw=0.7, ls=(0, (1, 1.5))); ax.text(xv + 0.03, 6, txt, fontsize=6.5, rotation=90, va="bottom", color="0.2")
    ax.annotate("1.42 Å\n(shared\nvacant site)", xy=(1.424, 20), xytext=(1.03, 8), fontsize=6.5,
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    ax.annotate("2.42 Å (geminal SiH$_2$)", xy=(2.417, 56), xytext=(2.62, 40), fontsize=6.5,
                arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    ax.set_xlim(1.0, 4.2); ax.set_ylim(0, 102)
    ax.set_xlabel("shortest H$\\cdots$H distance of each H (Å)"); ax.set_ylabel("cumulative fraction of H (%)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.98)); ax.grid(True, color="0.85", lw=0.5); ax.set_axisbelow(True)
    ax.set_title("$D$ = 1.5 nm", fontsize=7.5)
    fig.tight_layout(); ST.salvar(fig, os.path.join(fdir, "fig_distribuicao_HH.png")); plt.close(fig)

    # Fig 3: estequiometria
    fig, ax = plt.subplots(2, 1, figsize=(LARG, 4.3), sharex=True)
    for mv in (1, 2):
        x, y = serie(mv, lambda r: r["razao_H_Si"]); ax[0].plot(x, y, **estilo[mv])
        x, y = serie(mv, lambda r: 100 * r["fracao_Si_superficie"]); ax[1].plot(x, y, **estilo[mv])
    x, y = serie(2, lambda r: r["razao_H_Si"], "min_viz=2+preencher"); ax[0].plot(x, y, **EST_P, label="$n_{\\min}=2$ + filling")
    x, y = serie(2, lambda r: 100 * r["fracao_Si_superficie"], "min_viz=2+preencher"); ax[1].plot(x, y, **EST_P)
    ax[0].set_ylabel("H/Si ratio"); ax[1].set_ylabel("surface Si (%)"); ax[1].set_xlabel("nominal diameter $D$ (nm)")
    ax[0].legend(loc="upper right")
    ax[0].text(-0.22, 1.0, "(a)", transform=ax[0].transAxes, va="bottom", ha="left", fontsize=8, fontweight="bold")
    ax[1].text(-0.22, 1.0, "(b)", transform=ax[1].transAxes, va="bottom", ha="left", fontsize=8, fontweight="bold")
    for a in ax:
        a.grid(True, color="0.85", lw=0.5); a.set_axisbelow(True); a.set_xticks(DIAMS)
    fig.tight_layout(h_pad=0.4); ST.salvar(fig, os.path.join(fdir, "fig_estequiometria.png")); plt.close(fig)

    print("\nTempo total: %.1f s (threads=%s). Fonte das figuras: %s" % (time.time() - t_ini, os.environ["OMP_NUM_THREADS"], fam))

if __name__ == "__main__":
    main()
