#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gera as tabelas LaTeX do manuscrito ARI (em ingles, ponto decimal) a partir dos dados.

Entradas: dados/gaps.json, dados/relaxacao.json, dados/relaxacao_321g.json.
Saidas (tabelas/):
    tab_gaps.tex, tab_robustez.tex, tab_base.tex, tab_relax.tex  -- texto principal
    tab_budget.tex                                              -- balanco de sensibilidade (novo na versao ARI)
    tab_df.tex, tab_custo.tex                                   -- material suplementar
e dados/orcamento.json com os numeros derivados citados no texto (Secao 4.7).
Valores ausentes aparecem como '---'. Uso: python3 codigo/gerar_tabelas.py
Com --pt, grava a versao em portugues (virgula decimal, termos traduzidos) em tabelas_pt/, para o
manuscrito em portugues da equipe, sem regravar os JSON derivados.
"""
import os, re, sys, json
AQUI = os.path.dirname(os.path.abspath(__file__)); RAIZ = os.path.join(AQUI, '..')
PT = '--pt' in sys.argv
TAB = os.path.join(RAIZ, 'tabelas_pt' if PT else 'tabelas'); os.makedirs(TAB, exist_ok=True)
TXT = {'$\\nviz=2$+fill': '$\\nviz=2$+pr.', 'in-core': 'mem.', 'direct (GPU)': 'direto (GPU)', 'direct': 'direto',
       'Cutting criterion': 'Critério de corte', 'Filling': 'Preenchimento', 'Relaxation': 'Relaxação',
       'Basis set': 'Base', 'Density fitting': 'Ajuste de densidade', 'three systems': 'três sistemas', 'all': 'todas',
       'pristine': 'perfeito', 'unrelaxed': 'não relaxado', 'relaxed': 'relaxado'}
def t(x):
    return TXT.get(x, x) if PT else x
g = json.load(open(os.path.join(RAIZ, 'dados', 'gaps.json')))['calculos']
HC = 1239.84198   # eV nm

def f(x, nd=2):
    if not isinstance(x, (int, float)): return '---'
    s = '%.*f' % (nd, x)
    return s.replace('.', ',') if PT else s
def sinal(x):
    s = '%+.2f' % x
    return s.replace('.', ',') if PT else s
def formula_tex(fm):
    m = re.match(r'Si(\d+)H(\d+)', fm); return 'Si$_{%s}$H$_{%s}$' % m.groups()
def busca(D, tag, base, met, df=None):
    """Registro convergido para (D, tag, base, metodo); df=None: prefere exato, aceita DF."""
    cands = [(k, v) for k, v in g.items() if v.get('converged') and abs(v['D_nominal_nm'] - D) < 1e-6
             and v['base'] == base and v['metodo'] == met and k.split('|')[1] == tag]
    if df is True:  cands = [c for c in cands if c[0].endswith('|DF')]
    if df is False: cands = [c for c in cands if c[0].endswith('|exato')]
    if not cands: return None
    cands.sort(key=lambda c: c[0].endswith('|DF'))   # exato primeiro
    return cands[0][1] | {'_df': cands[0][0].endswith('|DF')}
def grava(nome, linhas):
    # sinal de menos tipografico: '-0.25' -> '$-$0.25' (nao afeta '3-21G' nem os intervalos '6.0--10.6')
    txt = re.sub(r'(?<=[\s(&])-(?=\d)', '$-$', '\n'.join(linhas))
    open(os.path.join(TAB, nome + '.tex'), 'w', encoding='utf-8').write(txt + '\n    \\bottomrule\n')
BASES = {'sto-3g': 'STO-3G', '3-21g': '3-21G', 'def2-svp': 'def2-SVP'}

# ---------------- Tabela: gaps por criterio (STO-3G)
L = []
CRIT = [('mv1', '$\\nviz\\le1$'), ('mv2', '$\\nviz=2$'), ('mv2p', t('$\\nviz=2$+fill'))]
for D in (0.8, 0.9, 1.0, 1.2, 1.5):
    first = True
    for tag, rot in CRIT:
        if tag == 'mv2p' and D == 0.9:      # o preenchimento em 0,9 nm produz a mesma estrutura que D=1,0/min_viz=2
            hf = busca(1.0, 'mv2', 'sto-3g', 'HF'); pbe = busca(1.0, 'mv2', 'sto-3g', 'PBE')
        else:
            hf = busca(D, tag, 'sto-3g', 'HF'); pbe = busca(D, tag, 'sto-3g', 'PBE')
        r = hf or pbe
        if r is None:
            if tag == 'mv2p' and D in (0.8, 1.0): continue
            L.append('    %s & %s & --- & --- & --- & --- & --- & --- \\\\' % (f(D, 1) if first else '', rot)); first = False; continue
        dfm = '$^{a}$' if (hf and hf['_df']) or (pbe and pbe['_df']) else ''
        L.append('    %s & %s & %s & %d & %s%s & %s%s & %s & %s \\\\' % (
            f(D, 1) if first else '', rot, formula_tex(r['formula']), r['nao'],
            f(hf['gap_eV']) if hf else '---', dfm if hf else '', f(pbe['gap_eV']) if pbe else '---', dfm if pbe else '',
            f(r['gap_brus_Defetivo_eV']), f(r['gap_delerue_Defetivo_eV'])))
        first = False
grava('tab_gaps', L)

# ---------------- Tabela: base (Si10H16 e Si22H28, min_viz=2)
L = []
for D, fm in ((0.8, 'Si$_{10}$H$_{16}$'), (1.0, 'Si$_{22}$H$_{28}$')):
    first = True
    for base in ('sto-3g', '3-21g', 'def2-svp'):
        hf = busca(D, 'mv2', base, 'HF'); pbe = busca(D, 'mv2', base, 'PBE'); r = hf or pbe
        if r is None:
            L.append('    %s & %s & --- & --- & --- & --- & --- \\\\' % (fm if first else '', BASES[base])); first = False; continue
        dfm = '$^{a}$' if r['_df'] else ''
        L.append('    %s & %s%s & %d & %s & %s & %s & %s \\\\' % (fm if first else '', BASES[base], dfm, r['nao'],
                 f(hf['gap_eV']) if hf else '---', f(pbe['gap_eV']) if pbe else '---',
                 f(hf['tempo_s'], 0) if hf else '---', f(pbe['tempo_s'], 0) if pbe else '---'))
        first = False
grava('tab_base', L)

# ---------------- Tabela S1: erro do ajuste de densidade
L = []; erro_df_max = 0.0
for D, base, fm in ((0.8, 'def2-svp', 'Si$_{10}$H$_{16}$'), (1.0, '3-21g', 'Si$_{22}$H$_{28}$'), (1.2, 'sto-3g', 'Si$_{36}$H$_{40}$')):
    for met in ('HF', 'PBE'):
        ex = busca(D, 'mv2', base, met, df=False); dfr = busca(D, 'mv2', base, met, df=True)
        if ex and dfr:
            erro_df_max = max(erro_df_max, abs(dfr['gap_eV'] - ex['gap_eV']))
            L.append('    %s & %s & %s & %s & %s & %s & %s \\\\' % (fm, BASES[base], met,
                     f(ex['gap_eV'], 3), f(dfr['gap_eV'], 3), f((dfr['gap_eV'] - ex['gap_eV']) * 1000, 1), f((dfr['e_tot_Ha'] - ex['e_tot_Ha']) * 1000, 2)))
        else:
            L.append('    %s & %s & %s & --- & --- & --- & --- \\\\' % (fm, BASES[base], met))
grava('tab_df', L)

# ---------------- Tabela S2: custo (min_viz=2, STO-3G)
L = []
for D in (0.8, 0.9, 1.0, 1.2, 1.5):
    hf = busca(D, 'mv2', 'sto-3g', 'HF'); pbe = busca(D, 'mv2', 'sto-3g', 'PBE'); r = hf or pbe
    if r is None: L.append('    %s & --- & --- & --- & --- & --- & --- & --- \\\\' % f(D, 1)); continue
    modo = 'DF' if r['_df'] else t('in-core' if r.get('modo_eri', '').startswith('incore') else ('direct (GPU)' if 'GPU' in r.get('modo_eri', '') else 'direct'))
    L.append('    %s & %s & %d & %d & %d & %s & %s & %s \\\\' % (f(D, 1), formula_tex(r['formula']), r['n_si'] + r['n_h'], r['nelectron'], r['nao'],
             f(hf['tempo_s'], 0) if hf else '---', f(pbe['tempo_s'], 0) if pbe else '---', modo))
grava('tab_custo', L)

# ---------------- Tabela: robustez do efeito do criterio (STO-3G vs 3-21G)
# Delta = gap(min_viz<=1) - gap(min_viz=2) no mesmo diametro nominal; compara-se com o previsto
# so pelo tamanho, via Eq. de Delerue avaliada nos dois D_efetivos.
L = []
for D in (0.8, 1.0, 1.2):
    linha_feita = False
    for base in ('sto-3g', '3-21g'):
        cel = []
        for met in ('HF', 'PBE'):
            a = busca(D, 'mv1', base, met); b = busca(D, 'mv2', base, met)
            cel.append(f(a['gap_eV'] - b['gap_eV']) if (a and b) else '---')
        a = busca(D, 'mv1', base, 'HF'); b = busca(D, 'mv2', base, 'HF')
        if a and b:
            prev = a['gap_delerue_Defetivo_eV'] - b['gap_delerue_Defetivo_eV']
            par = '%s$\\to$%s' % (formula_tex(a['formula']), formula_tex(b['formula']))
        else:
            prev, par = None, '---'
        L.append('    %s & %s & %s & %s & %s & %s \\\\' % (
            f(D, 1) if not linha_feita else '', par if not linha_feita else '', BASES[base],
            cel[0], cel[1], f(prev) if prev is not None else '---'))
        linha_feita = True
grava('tab_robustez', L)

# ---------------- Tabela: relaxacao
L = []; RELAX = {}
for arq, base in (('relaxacao.json', 'STO-3G'), ('relaxacao_321g.json', '3-21G')):
    p = os.path.join(RAIZ, 'dados', arq)
    if not os.path.exists(p): continue
    S = json.load(open(p))['sistemas']; RELAX[base] = S
    for rot, s in S.items():
        gi, gr, d, o = s['geom_ideal'], s['geom_relax'], s['delta'], s['otimizacao']
        L.append('    %s & %s & %d & %s & %s & %s & %s & %s & %s$\\to$%s & %d$\\to$%d & %s/%s & %s (%s) & %s (%s) \\\\' % (
            formula_tex(rot), base, o['n_passos'], f(d['dE_PBE_eV'], 2), f(d['dE_PBE_meV_por_H'], 0),
            f(gr['min_SiSi']) + '--' + f(gr['max_SiSi']), f(gr['min_SiH']) + '--' + f(gr['max_SiH']),
            f(gr['ang_HSiH_media'], 1) if gr['ang_HSiH_media'] else '---',
            f(gi['min_HH']), f(gr['min_HH']), gi['n_HH_lt_2A'], gr['n_HH_lt_2A'],
            f(gr['RMSD_Si_A']), f(gr['RMSD_H_A']),
            f(s['relax']['gap_HF_eV']), sinal(d['dgap_HF_eV']),
            f(s['relax']['gap_PBE_eV']), sinal(d['dgap_PBE_eV'])))
grava('tab_relax', L)

# ---------------- Tabela: balanco de sensibilidade (budget) e propagacao para o detector
# Para cada fonte: |dEg| em HF e PBE e o efeito relativo |dEg|/Eg_ref, com Eg_ref = Eq. (2) de
# Delerue avaliada no D_ef da estrutura de referencia. Em primeira ordem, esse efeito relativo e o
# mesmo no limite de rendimento luminoso Y = 1e6/(beta Eg), na energia por par (Klein) e no
# comprimento de onda de emissao (lambda = hc/E).
linhas_b, orc = [], []
def item(fonte, sistema, base, dHF, dPBE, ref, nota=''):
    rel = [abs(x) / ref * 100 for x in (dHF, dPBE) if x is not None]
    orc.append(dict(fonte=fonte, sistema=re.sub(r'[\$\{\}_]', '', sistema), base=base, dHF_eV=dHF, dPBE_eV=dPBE,
                    Eg_ref_eV=ref, rel_min_pct=min(rel), rel_max_pct=max(rel),
                    lambda_ref_nm=HC / ref, dlambda_max_nm=HC / ref * max(rel) / 100))
    faixa = f(min(rel), 1) if abs(max(rel) - min(rel)) < 0.05 else '%s--%s' % (f(min(rel), 1), f(max(rel), 1))
    def cel(x):
        if x is None: return '---'
        return '$<$' + f(0.005, 3) if abs(x) < 0.005 else f(abs(x))
    if max(rel) < 0.1: faixa = '$<$' + f(0.1, 1)
    linhas_b.append('    %s & %s & %s%s & %s & %s & %s & %s \\\\' % (t(fonte), t(sistema), t(base), nota, cel(dHF), cel(dPBE), f(ref), faixa))
def delta(D1, t1, D2, t2, base):
    out = []
    for met in ('HF', 'PBE'):
        a = busca(D1, t1, base, met); b = busca(D2, t2, base, met)
        out.append(a['gap_eV'] - b['gap_eV'] if (a and b) else None)
    ref = busca(D2, t2, base, 'HF') or busca(D2, t2, base, 'PBE')
    return out, ref
SUP = '$^{b}$'
primeiro = True
for D, base in ((0.8, '3-21g'), (1.0, '3-21g'), (1.2, '3-21g'), (1.5, 'sto-3g')):
    (dh, dp), r = delta(D, 'mv1', D, 'mv2', base)
    item('Cutting criterion' if primeiro else '', '%s (%s~nm)' % (formula_tex(r['formula']), f(D, 1)), BASES[base], dh, dp,
         r['gap_delerue_Defetivo_eV'], SUP if base == 'sto-3g' else ''); primeiro = False
linhas_b.append('    \\addlinespace')
primeiro = True
for D, base in ((1.2, '3-21g'), (1.5, 'sto-3g')):
    (dh, dp), r = delta(D, 'mv2p', D, 'mv2', base)
    item('Filling' if primeiro else '', '%s (%s~nm)' % (formula_tex(r['formula']), f(D, 1)), BASES[base], dh, dp,
         r['gap_delerue_Defetivo_eV'], SUP if base == 'sto-3g' else ''); primeiro = False
linhas_b.append('    \\addlinespace')
primeiro = True
for base in ('3-21G', 'STO-3G'):
    for rot, s in RELAX.get(base, {}).items():
        Dn = s['D_nm']; r = busca(Dn, 'mv2', 'sto-3g', 'HF')   # Eg_ref depende so de D_ef
        item('Relaxation' if primeiro else '', formula_tex(rot), base, s['delta']['dgap_HF_eV'], s['delta']['dgap_PBE_eV'],
             r['gap_delerue_Defetivo_eV'], SUP if base == 'STO-3G' else ''); primeiro = False
linhas_b.append('    \\addlinespace')
primeiro = True
for (b1, b2) in (('sto-3g', '3-21g'), ('3-21g', 'def2-svp')):
    for D in (0.8, 1.0):
        (dh, dp), r = delta(D, 'mv2', D, 'mv2', b1)   # so para obter a referencia
        vals = []
        for met in ('HF', 'PBE'):
            a = busca(D, 'mv2', b1, met); b = busca(D, 'mv2', b2, met)
            vals.append(b['gap_eV'] - a['gap_eV'] if (a and b) else None)
        item('Basis set' if primeiro else '', formula_tex(r['formula']), '%s$\\to$%s' % (BASES[b1], BASES[b2]), vals[0], vals[1],
             r['gap_delerue_Defetivo_eV']); primeiro = False
linhas_b.append('    \\addlinespace')
r = busca(0.8, 'mv2', 'sto-3g', 'HF')
item('Density fitting', 'three systems', 'all', erro_df_max, None, r['gap_delerue_Defetivo_eV'], '$^{c}$')
grava('tab_budget', linhas_b)
if not PT:
    json.dump(dict(nota='Numeros derivados da tabela tab_budget (gerar_tabelas.py). rel = |dEg|/Eg_ref*100; '
                        'Eg_ref = Eq. de Delerue no D_ef da estrutura de referencia; lambda = hc/Eg_ref.',
                   itens=orc), open(os.path.join(RAIZ, 'dados', 'orcamento.json'), 'w'), indent=1, ensure_ascii=False)

print('tabelas geradas em', TAB)
for n in ('tab_gaps', 'tab_base', 'tab_robustez', 'tab_relax', 'tab_budget', 'tab_df', 'tab_custo'):
    print('--', n); print(open(os.path.join(TAB, n + '.tex')).read())

# ---------------- Tabela: ligacao pendente de superficie (defeito induzido por radiacao)
# Le dados/defeito_<sistema>.json (codigo/defeito_ligacao_pendente.py). Gap efetivo do radical =
# min(gap alfa, gap beta); Delta = gap efetivo - gap do aglomerado perfeito relaxado (mesmo metodo).
L = []; DEF = {}
for rot in ('Si10H16', 'Si22H28'):
    p = os.path.join(RAIZ, 'dados', 'defeito_%s.json' % rot)
    if not os.path.exists(p): continue
    d = json.load(open(p)); DEF[rot] = d
    per = d['perfeito']
    if L: L.append('    \\addlinespace')
    L.append('    %s & %s & & & & & %s & --- & --- & %s & --- & --- & --- \\\\' % (
        formula_tex(rot), t('pristine'), f(per['PBE']['gap_eV']), f(per['HF']['gap_eV'])))
    for ch, s in sorted(d.get('sitios', {}).items(), key=lambda t: int(t[0][6:])):
        grupo = 'SiH' if s['tipo'] == 'SiH1' else 'SiH$_2$'
        primeira = True
        for geo in ('vertical', 'relaxado'):
            if geo not in s: continue
            r = s[geo]; pb, hf = r['PBE'], r['HF']
            cab = ('S%s & %s & %d & %s' % (ch[6:], grupo, s['n_equivalentes'], f(s['dist_Si_centro_A']))) if primeira else ' & & & '
            L.append('    & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\' % (
                cab, t('unrelaxed' if geo == 'vertical' else 'relaxed'),
                f(pb['gap_efetivo_eV']), sinal(pb['gap_efetivo_eV'] - per['PBE']['gap_eV']),
                f(pb['DB_ocupado_acima_Ev_eV']), f(hf['gap_efetivo_eV']), sinal(hf['gap_efetivo_eV'] - per['HF']['gap_eV']),
                f(pb['pop_spin_Si']), f(r['BDE_PBE_eV']) if geo == 'relaxado' else '---'))
            primeira = False
if L:
    grava('tab_defeito', L)
    # numeros derivados para o texto
    resumo = {}
    for rot, d in DEF.items():
        per = d['perfeito']; vs = [s for s in d['sitios'].values()]
        def faixa(geo, met, chave):
            v = [s[geo][met][chave] for s in vs if geo in s]
            return (min(v), max(v)) if v else None
        resumo[rot] = dict(
            gap_perfeito_PBE=per['PBE']['gap_eV'], gap_perfeito_HF=per['HF']['gap_eV'],
            gap_ef_vertical_PBE=faixa('vertical', 'PBE', 'gap_efetivo_eV'), gap_ef_relaxado_PBE=faixa('relaxado', 'PBE', 'gap_efetivo_eV'),
            gap_ef_vertical_HF=faixa('vertical', 'HF', 'gap_efetivo_eV'), gap_ef_relaxado_HF=faixa('relaxado', 'HF', 'gap_efetivo_eV'),
            spin_Si=faixa('vertical', 'PBE', 'pop_spin_Si'),
            BDE_relaxado=[s['relaxado']['BDE_PBE_eV'] for s in vs if 'relaxado' in s],
            E_relax=[s['relaxado']['E_relaxacao_eV'] for s in vs if 'relaxado' in s])
    if not PT: json.dump(resumo, open(os.path.join(RAIZ, 'dados', 'defeito_resumo.json'), 'w'), indent=1)
    print('-- tab_defeito'); print(open(os.path.join(TAB, 'tab_defeito.tex')).read())
