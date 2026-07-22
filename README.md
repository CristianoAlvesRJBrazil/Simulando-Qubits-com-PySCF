# Simulando Qubits com PySCF

### Uma Jornada do Zero aos Pontos Quânticos para a Computação Quântica

**Autor:** Cristiano Alves

Livro-texto que ensina, **do zero absoluto**, a modelar *qubits* de spin baseados
em pontos quânticos usando o [PySCF](https://pyscf.org) — uma plataforma de
química quântica aberta, gratuita e programável em Python. Partindo do primeiro
`import pyscf`, o leitor constrói, capítulo a capítulo, um *pipeline* completo de
simulação de um *qubit* de spin em **silício dopado com fósforo (Si:P)**.

> Todo o código do livro é **executável e reprodutível**. As energias, geometrias
> e figuras apresentadas foram geradas por cálculos reais com o PySCF.

---

## 📂 Estrutura do repositório

| Pasta | Conteúdo |
|-------|----------|
| [`livro/`](livro/) | Fonte LaTeX do livro (`Livro_PySCF.tex`), o PDF compilado e as figuras. |
| [`notebooks/`](notebooks/) | *Notebooks* Jupyter, um por capítulo, prontos para rodar no Colab ou localmente. |
| [`codigo/`](codigo/) | Módulos Python reutilizáveis (ex.: gerador de nanocristais). |
| [`proposta/`](proposta/) | Proposta editorial original da obra (sumário detalhado e plano). |

---

## 📖 Progresso da obra

**Parte 0 — Preparando o Laboratório Computacional**
- ✅ Capítulo 1 — Python Essencial para Química Quântica
- ✅ Capítulo 2 — O Primeiro Cálculo com PySCF

**Parte I — Estrutura Eletrônica de Pontos Quânticos**
- ✅ Capítulo 3 — Do Átomo ao Nanocristal: Construindo um Ponto Quântico
- ✅ Capítulo 4 — Hartree-Fock e o Diagrama de Orbitais do Ponto Quântico
- ⏳ Capítulo 5 — Teoria do Funcional da Densidade para Semicondutores

*(demais capítulos em elaboração — ver [`proposta/`](proposta/) para o plano completo)*

---

## 🔷 O Projeto Âncora

Um projeto longitudinal atravessa todo o livro: a construção passo a passo de um
*qubit* de spin em Si:P. A cada capítulo, uma etapa é acrescentada.

| Etapa | Marco | Capítulo |
|-------|-------|----------|
| 0 | Ambiente funcional + átomo de silício isolado | 1–2 |
| 1 | Geometria do nanocristal Si:P + orbital do doador (SOMO localizado) | 3–4 |
| 2+ | Parâmetros magnéticos, acoplamento de troca, descoerência | 5+ |

Resultado já obtido: o orbital que hospeda o *qubit* tem **83 % de sua densidade
num raio de 3 Å do fósforo** — a confirmação de primeiros princípios do
"hidrogênio artificial".

---

## 🚀 Como usar

### Notebooks (recomendado para começar)

Abra qualquer notebook de [`notebooks/`](notebooks/) no
[Google Colab](https://colab.research.google.com) (sem instalar nada) ou
localmente:

```bash
conda create -n pyscf-livro python=3.11
conda activate pyscf-livro
pip install pyscf numpy scipy matplotlib jupyter
jupyter notebook
```

### Compilar o livro

```bash
cd livro
latexmk -pdf Livro_PySCF.tex
```

Requer uma distribuição LaTeX (TeX Live) com os pacotes `tcolorbox`, `listings`,
`fncychap` e `newtxtt`.

---

## 🛠️ Software

Python 3.10+, PySCF 2.3+, NumPy, SciPy, Matplotlib, Jupyter; opcionalmente
`py3Dmol` para visualização 3D.

---

*Obra em desenvolvimento. Todos os direitos reservados ao autor.*
