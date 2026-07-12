# Arquitetura — QSAR de Ecotoxicidade para Formulações Cosméticas Amazônicas

## Objetivo

Este projeto tem como objetivo prever a **ecotoxicidade** (inibição do crescimento de algas) de ativos amazônicos e, futuramente, de formulações cosméticas, utilizando modelos **QSAR (Quantitative Structure–Activity Relationship)** baseados em descritores moleculares.

---

# Visão Geral da Arquitetura

```text
                        ┌─────────────────────────────┐
                        │   Base Pública ECOTOX (EPA) │
                        └──────────────┬──────────────┘
                                       │
                                       ▼
                              CAS → PubChem → SMILES
                                       │
                                       ▼
                           Cálculo de descritores (RDKit)
                                       │
                                       ▼
                       Treinamento do modelo QSAR (Fase 1)
                                       │
                         Random Forest + Validação k-fold
                                       │
                                       ▼
                           Modelo treinado (congelado)
                                       │
                 ┌─────────────────────┴─────────────────────┐
                 ▼                                           ▼
      Validação com ativos amazônicos            Futuras formulações
              (Fase 2)                    (Fase 3 — descritores de mistura)
```

---

# Motivação

Os dados disponíveis respondem à mesma pergunta científica em escalas diferentes:

- **Base pública:** uma molécula → um efeito biológico;
- **Base experimental:** um ingrediente → um efeito biológico.

Ainda não existe uma tabela contendo a composição das formulações cosméticas. Por esse motivo, não é possível construir um modelo para misturas.

A estratégia adotada segue a prática recomendada para modelos QSAR:

1. treinar utilizando uma base pública ampla;
2. congelar o modelo treinado;
3. realizar validação externa utilizando os dados experimentais próprios.

Essa abordagem reduz o risco de **overfitting** em uma base pequena (6–17 amostras) e está alinhada aos princípios de validação propostos pela OECD.

---

# Fase 1 — Treinamento com Base Pública

## Objetivo

Construir um modelo QSAR utilizando a base **ECOTOX Knowledgebase (US EPA)**.

## Componentes

| Etapa | Ferramenta / Fonte |
|--------|--------------------|
| Dados de toxicidade em algas | ECOTOX Knowledgebase (US EPA) |
| Conversão CAS → SMILES | PubChem PUG REST |
| Descritores moleculares | RDKit |
| Variável-resposta | pEC50 = −log10(EC50 em mol/L) |
| Modelo | Random Forest |
| Validação | Cross Validation (k = 5) |

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

---

## Observação

O download da base **ECOTOX** deve ser realizado localmente, conforme descrito em **`README.md`**, pois o ambiente de execução não possui acesso direto aos arquivos disponibilizados pela EPA.

---

# Fase 2 — Validação Externa

Após o treinamento, o modelo permanece **congelado**.

Não é realizado retreinamento utilizando os dados próprios, evitando que uma base pequena substitua o conhecimento adquirido na base pública.

## Dados utilizados

| Etapa | Fonte |
|--------|-------|
| Ingredientes | `Dados_QSAR_Saile_PJC2026.xlsx` |
| Aba utilizada | `5_SMILES_Ingredientes` |
| Descritores | RDKit (mesma função da Fase 1) |
| Comparação | Predição do modelo × % de inibição experimental |

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

# Fase 3 — Predição de Formulações (Pendente)

Atualmente não existe uma tabela relacionando:

- formulação;
- ingredientes;
- concentração de cada ingrediente.

Sem essas informações não é possível construir descritores moleculares de misturas.

A estrutura esperada é:

```text
amostra_id | ingrediente | percentual
```

Exemplo:

| amostra | ingrediente | % |
|----------|-------------|---:|
| A1 | Andiroba | 20 |
| A1 | Copaíba | 35 |
| A1 | Decyl glucoside | 45 |

Quando essa tabela estiver disponível será possível utilizar a função:

```python
descritor_da_mistura()
```

implementada em:

```text
dados.py
```

para gerar descritores ponderados da formulação completa.

---

# Problemas de Dados Conhecidos

## SMILES incorreto do Decyl Glucoside

Na aba **`5_SMILES_Ingredientes`**, o composto **Decyl glucoside** possui um SMILES incorreto.

O registro corresponde a um **sal de nióbio**, proveniente de uma página **Substance** do PubChem em vez da página **Compound**.

Como o RDKit interpreta esse SMILES sem gerar erro, essa inconsistência não é detectada automaticamente e deve ser corrigida manualmente.

---

## Diferenças de unidades experimentais

Os experimentos utilizam diferentes escalas de concentração celular:

- células/mL;
- ×10⁵ células/mL.

Essa padronização já é realizada automaticamente durante a normalização implementada em **`dados.py`**.

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
    ├── README.md
    ├── requirements.txt
    └── Dados_QSAR_Saile_PJC2026.xlsx
```

## Responsabilidades dos módulos

| Arquivo | Responsabilidade |
|----------|------------------|
| `main.py` | Ponto de entrada da aplicação (`dados` ou `publico`). |
| `descritores.py` | Cálculo dos descritores moleculares via RDKit. Módulo compartilhado pelos demais componentes. |
| `dados.py` | Processamento da planilha experimental, normalização dos dados e implementação inicial do descritor de mistura. |
| `ecotox.py` | Download, filtragem e preparação da base pública ECOTOX, incluindo consultas ao PubChem. |
| `modelo.py` | Ajuste de curvas dose–resposta (EC50), treinamento do modelo QSAR, validação cruzada e validação externa. |

A separação em módulos segue o princípio de **responsabilidade única**, organizando o código por função (descritores, dados, base pública e modelagem), e não pelas fases do projeto. Dessa forma, componentes reutilizados, como `descritores.py`, permanecem centralizados e evitam duplicação de código.

---

# Interface

O projeto possui interface exclusivamente em **linha de comando (CLI)**.

Os modos disponíveis são:

```bash
python3 main.py dados
```

Processa apenas os dados experimentais.

```bash
python3 main.py publico
```

Executa o treinamento utilizando a base pública ECOTOX e realiza a validação externa nos ativos amazônicos.

Não há interface gráfica, pois o objetivo do projeto é disponibilizar um pipeline científico modular voltado à pesquisa e experimentação computacional.
