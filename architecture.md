# Arquitetura — QSAR de Ecotoxicidade para Formulações Cosméticas Amazônicas

## Objetivo

Este projeto tem como objetivo prever a **ecotoxicidade** (inibição do crescimento de algas) de ativos amazônicos e, futuramente, de formulações cosméticas, utilizando modelos **QSAR (Quantitative Structure–Activity Relationship)** baseados em descritores moleculares.

O projeto integra o **TCC/PJC 2026 da Universidade Federal do Amapá (UNIFAP)** e foi desenvolvido integralmente em Python.

---

# Visão Geral

```text
                 Base Pública ECOTOX (EPA)
                           │
                           ▼
                 CAS → PubChem → SMILES
                           │
                           ▼
                  Descritores (RDKit)
                           │
                           ▼
              Treinamento do modelo QSAR
          Random Forest + Validação cruzada
                           │
                  Modelo congelado
                           │
         ┌─────────────────┴─────────────────┐
         ▼                                   ▼
 Validação externa                 Validação do modelo
 (ativos amazônicos)          (qualidade, Y-scrambling,
                                    leverage e GHS)
         │
         ▼
 Futuras formulações (Fase 3)
```

---

# Por que duas fases?

Os dados públicos e os dados experimentais respondem à mesma pergunta científica, porém em escalas diferentes.

- **Base pública:** uma molécula → um efeito biológico;
- **Base experimental:** um ingrediente → um efeito biológico.

Ainda não existe uma tabela contendo a composição completa das formulações cosméticas. Portanto, não é possível treinar diretamente um modelo para misturas.

A estratégia adotada segue a prática recomendada para modelos QSAR:

1. treinar utilizando uma base pública ampla;
2. congelar o modelo treinado;
3. validar externamente utilizando os dados experimentais.

Essa abordagem reduz o risco de **overfitting** em uma base pequena (6–17 amostras) e segue os princípios de validação da **OECD**.

---

# Fase 1 — Treinamento com Base Pública

## Objetivo

Construir um modelo QSAR utilizando dados públicos da **ECOTOX Knowledgebase (US EPA)**.

| Etapa | Ferramenta / Fonte |
|--------|--------------------|
| Dados de toxicidade | ECOTOX Knowledgebase |
| Conversão CAS → SMILES | PubChem PUG REST |
| Descritores | RDKit |
| Variável-resposta | pEC50 = −log10(EC50 em mol/L) |
| Modelo | Random Forest |
| Validação | Cross Validation (k = 5) |

### Funcionamento do Random Forest Regressor

![Diagrama de funcionamento do Random Forest Regressor para predição QSAR de ecotoxicidade](/home/paula/.gemini/antigravity-ide/brain/9c0943d7-ab88-45bf-9400-b187ec99615c/random_forest_qsar_diagram.png)

O **Random Forest Regressor** é um algoritmo de aprendizado de máquina supervisionado baseado no princípio de métodos ensemble (*ensemble learning*), que constrói e combina previsões de múltiplas árvores de decisão. O processo de funcionamento do algoritmo é descrito pelas seguintes etapas:

1. **Amostragem com Reposição (Bootstrap):**
   * A partir da base filtrada do ECOTOX contendo $N$ amostras, o algoritmo gera subconjuntos de dados aleatórios e com reposição (*bootstrap samples*). Cada árvore de decisão individual é treinada em um desses subconjuntos. Isso significa que uma mesma amostra pode ser usada múltiplas vezes para treinar uma árvore, enquanto outras não são utilizadas nessa árvore específica.

2. **Construção de Árvores de Decisão Descorrelacionadas:**
   * Uma árvore de regressão tenta subdividir o espaço de dados em regiões homogêneas em relação à variável-resposta ($pEC50$).
   * **Subespaço Aleatório (Feature Randomness):** Em cada nó da árvore, em vez de avaliar todos os descritores moleculares para decidir onde dividir os dados, o algoritmo sorteia aleatoriamente um subconjunto dos descritores disponíveis. A árvore escolhe o melhor divisor somente a partir desse subconjunto sorteado, o que diminui a correlação entre as árvores criadas e aumenta a diversidade do modelo.

3. **Predição por Agregação (Bagging):**
   * Uma vez treinadas as 300 árvores independentes (configuradas via `n_estimators=300`), o modelo final calcula a predição para um novo composto químico.
   * Cada árvore percorre seus nós de decisão até estimar um valor de $pEC50$ específico. A predição final do Random Forest é a **média aritmética simples** das estimativas de todas as 300 árvores individuais:
     $$\hat{y} = \frac{1}{B} \sum_{b=1}^{B} f_b(x)$$
     onde $B$ é o número de árvores e $f_b(x)$ é a predição da árvore $b$. Este processo reduz significativamente a variância geral sem elevar o viés do modelo.

4. **Importância dos Descritores (MDI):**
   * O algoritmo rastreia o quanto cada descritor molecular reduz a impureza dos nós (neste caso, a variância dos erros de regressão) ao longo de todas as árvores em que é selecionado. A média ponderada dessas reduções fornece a importância relativa de cada característica na tomada de decisão final do modelo.

---

## Fluxo

```text
ECOTOX
     │
     ▼
Filtragem dos testes
     │
     ▼
CAS Number
     │
     ▼
PubChem
     │
     ▼
SMILES
     │
     ▼
RDKit
     │
     ▼
Descritores
     │
     ▼
Random Forest
     │
     ▼
Modelo QSAR
```

O download da base ECOTOX deve ser realizado localmente seguindo as instruções presentes em **README.md**.

---

# Fase 2 — Validação Externa

Após o treinamento, o modelo permanece **congelado**, ou seja, não é ajustado novamente utilizando os dados experimentais.

Isso evita que uma base pequena sobrescreva o conhecimento aprendido na base pública.

| Etapa | Fonte |
|--------|-------|
| Ingredientes | `Dados_QSAR_Saile_PJC2026.xlsx` |
| Aba | `5_SMILES_Ingredientes` |
| Descritores | Mesma função RDKit utilizada na Fase 1 |
| Comparação | Predição do modelo × % de inibição observado |

---

## Fluxo

```text
Planilha experimental
          │
          ▼
       SMILES
          │
          ▼
RDKit
          │
          ▼
Mesmo conjunto de descritores
          │
          ▼
Modelo treinado
          │
          ▼
Predição de pEC50
          │
          ▼
Comparação com os resultados experimentais
```

---

# Validação do Modelo

A validação não constitui uma nova fase do projeto.

Ela funciona como uma **auditoria** realizada após o treinamento do modelo (Fase 1) e antes da interpretação dos resultados obtidos na validação externa (Fase 2).

O objetivo é verificar se o modelo realmente aprendeu relações entre estrutura molecular e atividade biológica, em vez de apenas memorizar padrões da base de treinamento.

| Verificação | Objetivo | Implementação |
|-------------|----------|---------------|
| Qualidade dos dados | Detectar duplicatas, outliers e descritores sem variação | `validacao.py` |
| Y-scrambling | Avaliar se o desempenho decorre de sinal químico real ou apenas de ruído | `validacao.py` |
| Domínio de aplicabilidade (Leverage) | Verificar se um ativo amazônico pertence ao espaço químico conhecido pelo modelo | `validacao.py` |
| Matriz de confusão (GHS) | Avaliar em quais categorias de toxicidade o modelo acerta ou erra | `validacao.py` |

---

## Fluxo de Validação

```text
Modelo treinado
        │
        ├────────► Qualidade dos dados
        │
        ├────────► Y-scrambling
        │
        ├────────► Domínio de aplicabilidade
        │
        └────────► Matriz de confusão (GHS)
```

Essas verificações **não modificam o modelo**.

Elas apenas fornecem evidências sobre sua confiabilidade.

### Matriz de Confusão

O modelo produz uma predição contínua em **pEC50**.

Para facilitar a interpretação ambiental, esses valores são convertidos em categorias de toxicidade aguda aquática do **Sistema Globalmente Harmonizado (GHS)**.

A matriz compara:

```text
Categoria observada
        ×
Categoria prevista
```

permitindo identificar se o modelo tende a:

- superestimar a toxicidade;
- subestimar a toxicidade;
- ou classificar corretamente cada composto.

---

# Fase 3 — Formulações Cosméticas (Pendente)

Atualmente não existe uma tabela relacionando cada formulação aos ingredientes que a compõem e às respectivas concentrações.

Sem essa informação não é possível calcular descritores moleculares de misturas.

A estrutura necessária é:

```text
amostra_id | ingrediente | percentual
```

Exemplo:

| amostra | ingrediente | % |
|----------|-------------|---:|
| A1 | Andiroba | 20 |
| A1 | Copaíba | 35 |
| A1 | Decyl glucoside | 45 |

Quando essa informação estiver disponível será utilizada a função:

```python
descritor_da_mistura()
```

presente em:

```text
dados.py
```

---

# Problemas de Dados Conhecidos

## SMILES incorreto do Decyl Glucoside

O SMILES presente na aba **`5_SMILES_Ingredientes`** corresponde a um sal de nióbio obtido de uma página **Substance** do PubChem, e não ao composto **Decyl glucoside**.

Como o RDKit interpreta esse SMILES sem gerar erro, a inconsistência não é detectada automaticamente.

A correção deve ser realizada manualmente.

---

## Diferenças nas unidades experimentais

Os experimentos utilizam diferentes escalas de concentração celular:

- células/mL;
- ×10⁵ células/mL.

Essa padronização já é realizada automaticamente durante a etapa de normalização implementada em **`dados.py`**.

---

# Organização do Código

```text
architecture.md
│
└── qsar_ecotox/
    ├── main.py
    ├── descritores.py
    ├── dados.py
    ├── ecotox.py
    ├── modelo.py
    ├── validacao.py
    ├── README.md
    ├── requirements.txt
    └── Dados_QSAR_Saile_PJC2026.xlsx
```

## Responsabilidades dos módulos

| Arquivo | Responsabilidade |
|----------|------------------|
| `main.py` | Ponto de entrada da aplicação (`dados` ou `publico`). |
| `descritores.py` | Cálculo dos descritores moleculares utilizando RDKit. |
| `dados.py` | Processamento da planilha experimental, normalização dos dados e implementação inicial do descritor de mistura. |
| `ecotox.py` | Preparação da base pública ECOTOX, filtragem dos dados e obtenção de SMILES via PubChem. |
| `modelo.py` | Ajuste de curvas dose–resposta, treinamento do modelo QSAR, validação cruzada e validação externa. |
| `validacao.py` | Avaliação da qualidade dos dados, Y-scrambling, domínio de aplicabilidade (Leverage) e matriz de confusão das categorias GHS. |

Cada módulo possui uma responsabilidade específica, seguindo o princípio de **responsabilidade única**. A organização é feita por função (descritores, dados, modelagem e validação), e não pelas fases do projeto, permitindo o reaproveitamento de componentes como `descritores.py` em diferentes etapas do pipeline.

---

# Interface

O projeto utiliza exclusivamente **interface de linha de comando (CLI)**.

Os modos disponíveis são:

```bash
python3 main.py dados
```

Processa apenas os dados experimentais.

```bash
python3 main.py publico
```

Executa o treinamento utilizando a base pública ECOTOX, realiza as rotinas de validação do modelo e aplica o modelo treinado aos ativos amazônicos.

Não há interface gráfica, pois o objetivo do projeto é disponibilizar um pipeline científico modular voltado à pesquisa e experimentação computacional.