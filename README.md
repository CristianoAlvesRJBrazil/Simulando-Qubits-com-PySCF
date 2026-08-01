<div align="center">

# Simulando Qubits com PySCF

### Uma Jornada do Zero aos Pontos Quânticos para a Computação Quântica

**Cristiano da Costa Alves** · **Nilseia Aparecida Barbosa** · **Fernando Manuel Araújo Moreira**

*Do primeiro `import pyscf` ao qubit de spin em Si:P*

[![Amostra](https://img.shields.io/badge/📥_Amostra_gratuita-63_páginas-1f6feb?style=for-the-badge)](livro/Livro_PySCF_Amostra.pdf)
[![Notebooks](https://img.shields.io/badge/📓_Notebooks-executáveis-f97316?style=for-the-badge)](notebooks/)

![Python](https://img.shields.io/badge/Python-3.10+-3776ab?logo=python&logoColor=white)
![PySCF](https://img.shields.io/badge/PySCF-2.3+-0b7285)
![Licença](https://img.shields.io/badge/livro-todos_os_direitos_reservados-6b7280)

</div>

---

## Sobre a obra

Modelar *qubits* deixou de ser curiosidade acadêmica para virar competência
estratégica — mas quem tenta entrar na área esbarra em pacotes de estrutura
eletrônica com curvas de aprendizado íngremes e licenças caras.

Este livro-texto ensina, **do zero absoluto**, a simular *qubits* de spin em
pontos quânticos com o [PySCF](https://pyscf.org): aberto, gratuito e
integralmente programável em Python. Não se pressupõe nenhuma experiência prévia
com química quântica computacional.

O diferencial é um **projeto longitudinal único**. Em vez de exemplos soltos, o
leitor constrói — capítulo a capítulo — um *pipeline* completo de simulação de um
*qubit* de spin em **silício dopado com fósforo (Si:P)**, do arquivo de geometria
até a estimativa do tempo de coerência.

> **Cada número deste livro saiu de um cálculo real.** Não há resultados
> ilustrativos: energias, geometrias, constantes de acoplamento e figuras foram
> gerados com o PySCF e podem ser reproduzidos pelo leitor. Quando o modelo
> encontra seus limites, o texto **diz isso com todas as letras** em vez de
> maquiar o resultado.

---

## 📥 Amostra gratuita

**[Baixe a amostra — 63 páginas, Partes 0 e I completas](livro/Livro_PySCF_Amostra.pdf)**

Cobre os capítulos 1 a 5: do ambiente Python até a Teoria do Funcional da
Densidade, incluindo a construção completa do nanocristal de Si:P. Ao terminar a
amostra você terá um ponto quântico atomístico pronto para calcular.

Este repositório hospeda a **amostra, os códigos e os notebooks**. O manuscrito
completo não é distribuído aqui.

---

## 📓 Notebooks executáveis

Um por capítulo, com saídas reais salvas. Rodam no **Google Colab** sem instalar
nada — clique no *badge*:

| Capítulo | Notebook | |
|---|---|---|
| 1 — Python Essencial para Química Quântica | [`ipynb`](notebooks/Capitulo_01_Python_Essencial.ipynb) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CristianoAlvesRJBrazil/Simulando-Qubits-com-PySCF/blob/main/notebooks/Capitulo_01_Python_Essencial.ipynb) |
| 2 — O Primeiro Cálculo com PySCF | [`ipynb`](notebooks/Capitulo_02_Primeiro_Calculo_PySCF.ipynb) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CristianoAlvesRJBrazil/Simulando-Qubits-com-PySCF/blob/main/notebooks/Capitulo_02_Primeiro_Calculo_PySCF.ipynb) |
| 3 — Do Átomo ao Nanocristal | [`ipynb`](notebooks/Capitulo_03_Do_Atomo_ao_Nanocristal.ipynb) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CristianoAlvesRJBrazil/Simulando-Qubits-com-PySCF/blob/main/notebooks/Capitulo_03_Do_Atomo_ao_Nanocristal.ipynb) |
| 4 — Hartree-Fock e o Diagrama de Orbitais ⏱️ | [`ipynb`](notebooks/Capitulo_04_Hartree_Fock_Orbitais.ipynb) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CristianoAlvesRJBrazil/Simulando-Qubits-com-PySCF/blob/main/notebooks/Capitulo_04_Hartree_Fock_Orbitais.ipynb) |
| 5 — Teoria do Funcional da Densidade ⏱️ | [`ipynb`](notebooks/Capitulo_05_DFT_Corrigindo_o_Gap.ipynb) | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CristianoAlvesRJBrazil/Simulando-Qubits-com-PySCF/blob/main/notebooks/Capitulo_05_DFT_Corrigindo_o_Gap.ipynb) |

<sub>⏱️ executa cálculos Hartree-Fock reais; a varredura do *gap* leva alguns minutos.</sub>


### 🎓 Aulas em PowerPoint

Apresentações prontas para sala de aula, com a identidade visual do livro e os
mesmos números calculados. Ver [`apresentacoes/`](apresentacoes/).

| Aula | Capítulo | |
|---|---|---|
| [`Cap. 2 — O Primeiro Cálculo`](apresentacoes/Aula_Capitulo_02_Primeiro_Calculo_PySCF.pptx) | 2 | 17 slides |

---

## 🔷 O Projeto Âncora

O fio condutor da obra. A cada capítulo, uma peça é acrescentada ao *pipeline*:

| Etapa | Marco | Cap. |
|:---:|---|:---:|
| 0 | Ambiente funcional + primeiro cálculo do átomo de silício | 1–2 |
| 1 | Geometria do nanocristal Si:P e o orbital do doador | 3–4 |
| 2 | *Gap* corrigido por DFT e densidade de spin do doador | 5 |
| 3 | Acoplamento hiperfino elétron–núcleo de ³¹P | 6 |
| 4 | Tensor *g* e a frequência de ressonância | 7 |
| 5 | Acoplamento de troca *J*: a porta de dois *qubits* | 8 |
| 6 | Tempo de coerência e purificação isotópica | 9 |
| 7 | Controle elétrico por eletrodo de porta | 10 |
| 8 | Registrador de dez *qubits* e a porta √SWAP | 11 |
| 9 | *Pipeline* reprodutível e relatório final | 12 |

### Resultados da amostra

Obtidos nos capítulos disponíveis publicamente:

- **Nanocristal Si₇₁P₁H₆₄** de 1,5 nm — 136 átomos, 1073 elétrons, coordenação 4
  em todos os silícios, sem ligações pendentes.
- **1073 elétrons — número ímpar.** O fósforo quebra a paridade, e o elétron
  desemparelhado resultante *é* o *qubit*.
- **O orbital do doador é localizado:** 83 % da densidade do SOMO num raio de
  3 Å do fósforo — a confirmação, por primeiros princípios, do "hidrogênio
  artificial" previsto pela intuição química.
- **A correlação eletrônica importa:** o *gap* HF/STO-3G de ~17 eV cai para ~8 eV
  com PBE (e ~10 eV com B3LYP), aproximando-se da faixa experimental sem perder
  a tendência de confinamento.

---

## 📖 Progresso da obra

**Parte 0 — Preparando o Laboratório Computacional**
- ✅ 1 — Python Essencial para Química Quântica
- ✅ 2 — O Primeiro Cálculo com PySCF

**Parte I — Estrutura Eletrônica de Pontos Quânticos**
- ✅ 3 — Do Átomo ao Nanocristal: Construindo um Ponto Quântico
- ✅ 4 — Hartree-Fock e o Diagrama de Orbitais do Ponto Quântico
- ✅ 5 — Teoria do Funcional da Densidade: Corrigindo o *Gap*

**Parte II — Os Parâmetros Magnéticos do Qubit**
- 📝 6 — O Acoplamento Hiperfino: o Diálogo Elétron–Núcleo
- 📝 7 — O Tensor *g*: a Frequência de Ressonância do Qubit

**Parte III — Dois Qubits e a Coerência**
- 📝 8 — O Acoplamento de Troca: Portas de Dois Qubits
- 📝 9 — Decoerência: Quanto Tempo o Qubit Sobrevive

**Parte IV — Controlando o Qubit**
- 📝 10 — O Eletrodo de Porta: Controle Elétrico do Qubit

**Parte V — Do Qubit ao Dispositivo**
- 📝 11 — Além do PySCF: do Parâmetro ao Dispositivo
- 📝 12 — O Projeto Âncora Consolidado

<sub>✅ disponível na amostra · 📝 redigido, não distribuído aqui · plano completo
na [proposta editorial](proposta/Livro_PySCF_Quantico.pdf)</sub>

---

## 📂 Estrutura do repositório

| Pasta | Conteúdo |
|---|---|
| [`livro/`](livro/) | Amostra pública: `Livro_PySCF_Amostra.pdf` e suas fontes LaTeX em [`livro/tex/`](livro/tex/). |
| [`notebooks/`](notebooks/) | *Notebooks* Jupyter, um por capítulo, com saídas reais. |
| [`codigo/`](codigo/) | Módulos Python reutilizáveis (ex.: gerador de nanocristais). |
| [`apresentacoes/`](apresentacoes/) | Aulas em PowerPoint (16:9) que acompanham os capítulos. |
| [`proposta/`](proposta/) | Proposta editorial: justificativa, público-alvo e sumário detalhado. |

---

## 🚀 Como usar

### Notebooks — o caminho mais rápido

Abra qualquer notebook pelo *badge* do Colab acima: nada a instalar. Localmente:

```bash
conda create -n pyscf-livro python=3.11
conda activate pyscf-livro
pip install pyscf numpy scipy matplotlib jupyter py3Dmol
jupyter notebook
```

### Do zero a um ponto quântico calculável

```python
import numpy as np
from codigo.nanocristal import construir, para_xyz
from pyscf import gto

si, h = construir(1.5)                              # nanocristal de 1,5 nm passivado
idx_p = int(np.argmin(np.linalg.norm(si, axis=1)))  # doador de P no centro

mol = gto.M(atom=para_xyz(si, h, idx_p), basis="sto-3g", spin=1)
print(f"Si{len(si)-1}P1H{len(h)} — {mol.nelectron} elétrons, {mol.nao} funções de base")
```

```
Si71P1H64 — 1073 elétrons, 712 funções de base
```

O `spin=1` não é detalhe técnico: são **1073 elétrons**, um número ímpar. O
elétron desemparelhado que sobra é o *qubit*.

### Compilar a amostra

```bash
cd livro
latexmk -pdf Livro_PySCF_Amostra.tex
```

Requer uma distribuição LaTeX (TeX Live) com `tcolorbox`, `listings`, `fncychap`,
`multirow` e `newtxtt`.

---

## 🛠️ Requisitos

Python 3.10+, [PySCF](https://pyscf.org) 2.3+, NumPy, SciPy, Matplotlib e
Jupyter. Opcionais: `py3Dmol` (visualização 3D interativa dos nanocristais) e
Avogadro/VMD/VESTA para inspecionar os arquivos `.xyz`.

O PySCF tem suporte pleno a **Linux** e **macOS**. No **Windows**, use o WSL ou
o Google Colab.

---

<div align="center">

*Obra em desenvolvimento · Todos os direitos reservados aos autores*

</div>
