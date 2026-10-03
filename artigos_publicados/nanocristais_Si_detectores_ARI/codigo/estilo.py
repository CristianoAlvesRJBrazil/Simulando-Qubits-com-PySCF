#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""estilo.py -- paleta e estilos compartilhados pelas figuras do artigo.

Versao ARI: figuras coloridas (a Elsevier publica cor online sem custo), mas a identidade de
cada serie continua redundante -- cor E marcador E estilo de traco --, de modo que a figura
segue legivel em impressao monocromatica e para leitores com daltonismo.

Paleta de Okabe-Ito, segura para daltonismo (deuteranopia e protanopia).
"""

import os

# --- modo de saida ----------------------------------------------------------
# FIG_MONO=1  -> paleta convertida para tons de cinza (versao impressa em preto e branco)
# FIG_DIR=... -> diretorio de saida das figuras (padrao: figuras)
MONO  = os.environ.get('FIG_MONO', '0') == '1'
SAIDA = os.environ.get('FIG_DIR', 'figuras')
# FIG_LANG=pt -> rotulos em portugues e virgula decimal (versao em portugues para a equipe).
# Os scripts escrevem os rotulos em ingles; salvar() os traduz pelo dicionario abaixo antes de gravar.
LANG  = os.environ.get('FIG_LANG', 'en')
if LANG == 'pt':
    import locale, matplotlib
    try:
        locale.setlocale(locale.LC_NUMERIC, 'pt_BR.UTF-8')
        matplotlib.rcParams['axes.formatter.use_locale'] = True     # 0,8 nos eixos
    except locale.Error:
        pass

# No modo monocromatico nao se converte por luminancia: cores distintas podem cair em
# cinzas proximos. Usa-se uma rampa explicita, bem separada, reproduzindo a versao impressa.
_RAMPA_MONO = {
    '#0072B2': '#000000',   # AZUL     -> preto
    '#D55E00': '#595959',   # VERMELHO -> cinza medio
    '#009E73': '#000000',   # VERDE    -> preto (distinto por traco e marcador)
    '#E69F00': '#737373',   # LARANJA  -> cinza claro
    '#CC79A7': '#595959',   # ROXO     -> cinza medio
    '#56B4E9': '#8C8C8C',   # CEU      -> cinza
    '#4D4D4D': '#333333',   # CINZA    -> cinza escuro
    '#BFBFBF': '#BFBFBF',   # CINZACLARO
}

def _c(hexcor):
    return _RAMPA_MONO[hexcor] if MONO else hexcor

# --- paleta base ------------------------------------------------------------
AZUL     = _c('#0072B2')   # metodo HF; criterio min_viz = 2
VERMELHO = _c('#D55E00')   # metodo PBE
VERDE    = _c('#009E73')   # criterio min_viz = 2 + preenchimento
LARANJA  = _c('#E69F00')   # criterio min_viz <= 1
ROXO     = _c('#CC79A7')   # ajuste de Delerue
CEU      = _c('#56B4E9')   # auxiliar
CINZA    = _c('#4D4D4D')   # equacao de Brus, eixos auxiliares
CINZACLARO = _c('#BFBFBF') # grades e linhas de referencia

# --- por criterio de corte --------------------------------------------------
CRITERIO = {
    1:    dict(color=LARANJA, ls='--', marker='s', mfc='white'),
    2:    dict(color=AZUL,    ls='-',  marker='o'),
    '2p': dict(color=VERDE,   ls='-.', marker='D', mfc='white'),
}

# --- por metodo -------------------------------------------------------------
METODO = {'HF': AZUL, 'PBE': VERMELHO}

# --- rcParams comuns --------------------------------------------------------
RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'lines.linewidth': 1.2, 'lines.markersize': 4,
    'savefig.dpi': 300, 'figure.dpi': 300,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'legend.frameon': False,
}


# --- gravacao ---------------------------------------------------------------
# A Elsevier pede >= 1000 dpi para graficos de linha em bitmap; por isso cada figura e gravada
# tambem em PDF vetorial (usado pelo LaTeX e enviado a revista), com fontes TrueType embutidas
# (pdf.fonttype = 42), alem do PNG de 600 dpi para visualizacao rapida.
def salvar(fig, caminho_png, **kw):
    import matplotlib
    kw.setdefault('facecolor', 'white')
    if LANG == 'pt':
        _traduzir_figura(fig)
    with matplotlib.rc_context({'pdf.fonttype': 42, 'ps.fonttype': 42}):
        fig.savefig(caminho_png, dpi=600, **kw)
        fig.savefig(caminho_png[:-4] + '.pdf', **kw)


# --- traducao para a versao em portugues ------------------------------------
# Trechos (substituidos do mais longo para o mais curto) e depois virgula decimal nos numeros.
_PT = {
    '(no effective filter)': '(sem filtro efetivo)', '+ filling': '+ preench.', 'filling': 'preenchimento',
    'H$\\cdots$H pairs with $d$ < 2.0 Å': 'pares H$\\cdots$H com $d$ < 2,0 Å', 'H atoms involved (%)': 'H envolvidos (%)',
    'nominal diameter $D$ (nm)': 'diâmetro nominal $D$ (nm)', '(shared\nvacant site)': '(sítio vago\ncomum)',
    'shortest H$\\cdots$H distance of each H (Å)': 'menor distância H$\\cdots$H de cada H (Å)',
    'cumulative fraction of H (%)': 'fração acumulada de H (%)', 'H/Si ratio': 'razão H/Si', 'surface Si (%)': 'Si de superfície (%)',
    'H$\\cdots$H pairs <': 'pares H$\\cdots$H <', 'no H$\\cdots$H pair <': 'nenhum par H$\\cdots$H <',
    'Brus, effective mass (Eq. 1)': 'Brus, massa efetiva (Eq. 1)', 'Delerue et al. (1993) fit (Eq. 2)': 'ajuste de Delerue et al. (1993) (Eq. 2)',
    'bulk Si': 'Si maciço', 'effective diameter $D_{\\mathrm{eff}}$ (nm)': 'diâmetro efetivo $D_{\\mathrm{ef}}$ (nm)',
    'HOMO$-$LUMO gap (eV)': 'gap HOMO$-$LUMO (eV)', 'size-only estimate (Eq. 2)': 'previsto só pelo tamanho (Eq. 2)',
    ' relaxed': ' relaxado', 'hydrogen (ascending order)': 'hidrogênio (ordem crescente)',
    'shortest H$\\cdots$H distance (Å)': 'menor distância H$\\cdots$H (Å)',
    '$E_g$(relaxed) $-$ $E_g$(ideal) (eV)': '$E_g$(relaxada) $-$ $E_g$(ideal) (eV)',
    'pristine': 'perfeito', 'Kohn–Sham level, PBE/3-21G (eV)': 'nível de Kohn–Sham, PBE/3-21G (eV)',
    'band edges': 'bordas de banda', 'dangling bond, occupied': 'ligação pendente, ocupado', 'dangling bond, empty': 'ligação pendente, vazio',
}
def traduzir(texto):
    if LANG != 'pt' or not texto:
        return texto
    import re
    for en in sorted(_PT, key=len, reverse=True):
        texto = texto.replace(en, _PT[en])
    return re.sub(r'(?<=\d)\.(?=\d)', ',', texto)          # 2.0 -> 2,0 (so entre digitos)

def _traduzir_figura(fig):
    import matplotlib.text as mtext
    import functools
    from matplotlib.ticker import FixedFormatter, FuncFormatter
    fig.canvas.draw()                                     # materializa rotulos e legendas
    for t in fig.findobj(mtext.Text):
        t.set_text(traduzir(t.get_text()))
    for ax in fig.axes:
        for eixo in (ax.xaxis, ax.yaxis):
            fmt = eixo.get_major_formatter()
            if isinstance(fmt, FixedFormatter):          # rotulos de categoria (matplotlib antigo)
                fmt.seq = [traduzir(x) for x in fmt.seq]
            elif isinstance(fmt, FuncFormatter) and isinstance(getattr(fmt, 'func', None), functools.partial) \
                    and fmt.func.args and isinstance(fmt.func.args[0], dict):   # set_xticklabels no matplotlib >= 3.8
                rot = fmt.func.args[0]
                for k in rot:
                    rot[k] = traduzir(rot[k])
