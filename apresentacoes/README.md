# Apresentações — *Simulando Qubits com PySCF*

Aulas em PowerPoint (16:9) que acompanham os capítulos do livro. Seguem a
identidade visual da obra: azul institucional, caixas didáticas coloridas
(💡 nota · ⚠️ atenção · ✓ conclusão · 🔷 Projeto Âncora) e blocos de código em
fonte monoespaçada.

## Conteúdo

| Arquivo | Capítulo | Slides |
|---|---|---|
| `Aula_Capitulo_02_Primeiro_Calculo_PySCF.pptx` | 2 — O Primeiro Cálculo com PySCF | 17 |

## Estrutura de uma aula

Cada apresentação segue o mesmo roteiro:

1. **Capa** e **roteiro** com os objetivos
2. **Desenvolvimento** — uma ideia por slide, alternando conceito e código
3. **Caixas de atenção** para os erros clássicos do capítulo
4. **Mãos à obra** — o exercício guiado
5. **🔷 Projeto Âncora** — o avanço do fio condutor do livro
6. **Resumo** e gancho para a aula seguinte

Todos os números exibidos (energias, funções de base, saídas de console) vêm dos
cálculos reais do livro — não são valores ilustrativos.

## Como regerar

As apresentações são construídas por *script* com
[`python-pptx`](https://python-pptx.readthedocs.io), o que mantém o padrão visual
consistente e permite reconstruí-las quando o texto do capítulo mudar.

```bash
pip install python-pptx
python build_ppt.py          # gera o .pptx
```

Para conferir o resultado sem abrir o PowerPoint:

```bash
libreoffice --headless --convert-to pdf Aula_Capitulo_02_Primeiro_Calculo_PySCF.pptx
```

## Aviso sobre os capítulos 6 a 12

Assim como o manuscrito, **aulas de capítulos não divulgados não devem ir para o
repositório público**. O `.gitignore` já protege o padrão
`apresentacoes/Aula_Capitulo_0[6-9]*` e `apresentacoes/Aula_Capitulo_1[0-2]*`.
