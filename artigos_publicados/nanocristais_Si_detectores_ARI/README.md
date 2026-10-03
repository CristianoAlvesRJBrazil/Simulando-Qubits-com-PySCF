# Modelos atomísticos de pontos quânticos de silício para a simulação de detectores de radiação ionizante

Códigos Python/PySCF da pesquisa do artigo:

> C. da Costa Alves, L. Oliveira de Sá, N. A. Barbosa, F. M. Araújo-Moreira.
> *Atomistic models of silicon quantum dots for simulating ionizing radiation detectors: an open construction
> protocol and the sensitivity of the electronic gap to geometric choices.*
> Submetido à *Applied Radiation and Isotopes* (Elsevier).

**Situação:** submetido. O link para o artigo publicado entra aqui quando ele sair.

*English summary: Python/PySCF scripts behind the article above. They build hydrogen-passivated silicon
nanocrystals from the diamond lattice, validate the structures, and compute HOMO–LUMO gaps, geometry relaxations,
and surface dangling-bond defects. Script names, comments, and data keys are in Portuguese; the glossary below
translates the terms needed to read the data.*

## O que a pesquisa faz

O artigo apresenta um protocolo aberto para construir aglomerados Si<sub>n</sub>H<sub>m</sub> a partir da rede
diamante:

- corte esférico;
- filtro de coordenação mínima;
- passivação direcional com hidrogênio;
- critério de preenchimento, que elimina pares H–H espúrios.

O artigo também mede quanto o *gap* HOMO–LUMO depende dessas escolhas e compara essa sensibilidade com a
do método, da base, do ajuste de densidade e da relaxação de geometria. O efeito é propagado às grandezas
usadas na simulação de detectores. Por fim, os scripts constroem os modelos do defeito típico induzido por
radiação: a ligação pendente que fica quando um hidrogênio é removido.

O gerador de nanocristais amplia o do Capítulo 3 do livro ([`codigo/nanocristal.py`](../../codigo/nanocristal.py)).

## Scripts

| Script | O que faz | Grava |
|---|---|---|
| `analise_geometria.py` | Gera os nanocristais (0,8–2,0 nm; n<sub>min</sub> = 0–4 e com preenchimento), valida as estruturas e analisa os contatos H–H. Leva segundos | `dados/geometria_sensibilidade.*`, `dados/xyz/`, figuras de contatos e estequiometria |
| `varredura_gap.py` | *Gap* HF e PBE em função do critério de corte (STO-3G) e da base (STO-3G, 3-21G, def2-SVP), na CPU | `dados/gaps.json`, `dados/gaps.csv` |
| `varredura_extra.py` | Critério de preenchimento a 1,2 nm e erro do ajuste de densidade em STO-3G | `dados/gaps.json` |
| `varredura_gpu.py` | Nanocristais de 1,5 nm na GPU (gpu4pyscf), por SCF direto | `dados/gaps.json` |
| `varredura_321g.py` | Efeito do critério de corte repetido na base 3-21G (GPU) | `dados/gaps.json` |
| `relaxacao.py` | Relaxação de geometria PBE/STO-3G (geomeTRIC) e métricas antes/depois | `dados/relaxacao.json`, XYZ |
| `relaxacao_321g.py` | Relaxação de controle em PBE/3-21G | `dados/relaxacao_321g.json`, XYZ |
| `defeito_ligacao_pendente.py` | Ligação pendente de superfície: remove um H de cada sítio distinto por simetria e calcula o radical (UKS-PBE e UHF), sem e com relaxação | `dados/defeito_*.json`, `dados/defeito_xyz/` |
| `gerar_tabelas.py` | Tabelas do artigo a partir dos dados | `tabelas/*.tex` |
| `fig_*.py`, `estilo.py` | Figuras do artigo a partir dos dados (PDF vetorial e PNG de 600 dpi) | `figuras/` |

## Requisitos

- Python 3.11, NumPy 2.4, SciPy 1.16, PySCF 2.14, geomeTRIC 1.1 e Matplotlib 3.10.
- Para os scripts de GPU: gpu4pyscf 1.8.1 com CUDA. Os cálculos do artigo usaram uma placa de 4 GB.

```bash
pip install numpy scipy matplotlib pyscf geometric
pip install gpu4pyscf-cuda12x      # opcional, so para os scripts de GPU
```

## Como reproduzir

Rode os comandos a partir desta pasta. Os scripts de cálculo gravam um registro por cálculo, de forma
incremental, e uma execução interrompida pode ser retomada. Nenhuma etapa é aleatória: com as mesmas
versões de software, os resultados numéricos se repetem bit a bit.

```bash
OMP_NUM_THREADS=1 python3 codigo/analise_geometria.py           # estruturas e contatos H-H (segundos)
python3 codigo/varredura_gap.py --etapas a b c c2               # gaps na CPU
python3 codigo/varredura_extra.py --etapas e1 e2
python3 codigo/varredura_gpu.py                                 # 1,5 nm na GPU
python3 codigo/varredura_321g.py                                # criterio de corte em 3-21G (GPU)
python3 codigo/relaxacao.py Si10H16 Si19H28
python3 codigo/relaxacao_321g.py Si10H16 Si19H28
python3 codigo/defeito_ligacao_pendente.py Si10H16 --gpu --relaxar todas --def2svp
python3 codigo/defeito_ligacao_pendente.py Si22H28 --gpu --relaxar SiH
python3 codigo/gerar_tabelas.py                                 # tabelas
for f in fig_estruturas fig_gap_vs_D fig_robustez fig_relaxacao fig_defeito; do python3 codigo/$f.py; done
```

O conjunto completo de dados e geometrias acompanha o artigo como material suplementar e será depositado em
repositório público com DOI após a publicação.

## Glossário

| Português | Inglês |
|---|---|
| `min_viz`, `mv2` | critério de coordenação mínima n<sub>min</sub>; n<sub>min</sub> = 2 |
| `preencher`, `mv2p` | critério de preenchimento (*filling*); n<sub>min</sub> = 2 com preenchimento |
| `rede_diamante`, `vizinhos` | rede diamante (*diamond lattice*), vizinhos (*neighbors*) |
| `D_nominal_nm`, `D_efetivo_nm` | diâmetro nominal e diâmetro efetivo D<sub>eff</sub> |
| `base`, `metodo`, `modo_eri` | base (*basis set*), método, modo das integrais (*in core*, direto, ajuste de densidade) |
| `perfeito` | aglomerado sem defeito (*pristine*) |
| `vertical`, `relaxado` | radical sem relaxação (*unrelaxed*) e relaxado (*relaxed*) |
| `sitios`, `classe` | sítios, classe de simetria |
| `gap_efetivo_eV`, `DB_ocupado_eV`, `DB_vazio_eV` | *gap* efetivo; níveis ocupado e vazio da ligação pendente |
| `orcamento` | balanço de sensibilidade (*sensitivity budget*) |
