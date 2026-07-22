# Notebooks — *Simulando Qubits com PySCF*

*Notebooks* Jupyter que acompanham o livro, capítulo a capítulo. Cada um é
autocontido e pode ser executado no **Google Colab** (sem instalação) ou
localmente no **Jupyter** (via Anaconda).

## Conteúdo

| Notebook | Capítulo |
|----------|----------|
| `Capitulo_01_Python_Essencial.ipynb` | 1 — Python Essencial para Química Quântica |
| `Capitulo_02_Primeiro_Calculo_PySCF.ipynb` | 2 — O Primeiro Cálculo com PySCF |

*(novos capítulos serão adicionados aqui)*

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
