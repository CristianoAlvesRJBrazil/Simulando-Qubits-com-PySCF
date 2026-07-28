# Notebooks — *Simulando Qubits com PySCF*

*Notebooks* Jupyter que acompanham o livro, capítulo a capítulo. Cada um é
autocontido e pode ser executado no **Google Colab** (sem instalação) ou
localmente no **Jupyter** (via Anaconda).

## Conteúdo

| Notebook | Capítulo |
|----------|----------|
| `Capitulo_01_Python_Essencial.ipynb` | 1 — Python Essencial para Química Quântica |
| `Capitulo_02_Primeiro_Calculo_PySCF.ipynb` | 2 — O Primeiro Cálculo com PySCF |
| `Capitulo_03_Do_Atomo_ao_Nanocristal.ipynb` | 3 — Do Átomo ao Nanocristal: Construindo um Ponto Quântico |
| `Capitulo_04_Hartree_Fock_Orbitais.ipynb` | 4 — Hartree-Fock e o Diagrama de Orbitais do Ponto Quântico |
| `Capitulo_05_DFT_Corrigindo_o_Gap.ipynb` | 5 — Teoria do Funcional da Densidade: Corrigindo o *Gap* |

*(novos capítulos serão adicionados aqui)*

⏱️ Os *notebooks* dos Capítulos 4 e 5 executam cálculos pesados; o do Cap. 5 é o
mais longo (três métodos — HF, PBE e B3LYP — em quatro tamanhos). Em ambos, para
um teste rápido, remova o `1.2` da lista de diâmetros.

⏱️ O *notebook* do Capítulo 4 executa cálculos Hartree-Fock reais e é o mais
demorado até aqui: a varredura do *gap* leva alguns minutos (o ponto de 1,2 nm,
com 364 funções de base, responde pela maior parte). Para um teste rápido, remova
o `1.2` da lista de diâmetros. Ele regenera a figura `gap_vs_tamanho.png` usada
no livro.

O *notebook* do Capítulo 3 usa também **matplotlib** (gráfico do confinamento) e
**py3Dmol** (visualização 3D interativa do nanocristal). Ambos são opcionais: a
célula de visualização degrada com elegância caso o py3Dmol não esteja instalado.
Ele gera o arquivo de geometria `SiP_1.5nm.xyz`, reutilizado nos capítulos
seguintes.

## Como executar

**No Colab:** faça upload do `.ipynb` (ou abra a partir do GitHub) e execute as
células na ordem com `Shift+Enter`. Descomente a linha `!pip install pyscf` na
primeira vez.

**Localmente (Anaconda):**

```bash
conda create -n pyscf-livro python=3.11
conda activate pyscf-livro
pip install pyscf numpy scipy matplotlib jupyter
jupyter notebook
```

## Projeto Âncora

As células marcadas com 🔷 **Projeto Âncora** constroem, ao longo do livro, um
*pipeline* reprodutível de simulação de um *qubit* de spin em silício dopado com
fósforo (Si:P) — do primeiro `import pyscf` até a estimativa do tempo de
descoerência.
