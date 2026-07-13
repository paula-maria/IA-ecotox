# QSAR Ecotoxicidade — Ativos Amazônicos

Pipeline em **Python** para prever a **ecotoxicidade** (inibição do crescimento de algas) de ativos cosméticos amazônicos utilizando modelos **QSAR (Quantitative Structure–Activity Relationship)** baseados em descritores moleculares.

O projeto está organizado em módulos independentes para facilitar a manutenção, reutilização do código e evolução para futuras etapas de modelagem de misturas. Para compreender a arquitetura completa e o fluxo das três fases do projeto, consulte **`architecture.md`**.

---

# Instalação

Requer **Python 3.10** ou superior.

## Instalação direta

```bash
pip install -r requirements.txt --break-system-packages
```

## Ambiente virtual (recomendado)

```bash
python3 -m venv venv

# Linux/macOS
source venv/bin/activate

# Windows
venv\Scripts\activate

pip install -r requirements.txt
```

---

## RDKit

O projeto utiliza o **RDKit** para:

- interpretar estruturas moleculares em formato SMILES;
- calcular descritores moleculares;
- gerar as variáveis de entrada dos modelos QSAR.

Na maioria dos sistemas basta instalar via:

```bash
pip install rdkit
```

Caso a instalação não funcione:

```bash
conda install -c conda-forge rdkit
```

---

# Estrutura do Projeto

```text
.
├── architecture.md
├── README.md
├── requirements.txt
├── Dados_QSAR_Saile_PJC2026.xlsx
├── descritores.py
├── dados.py
├── ecotox.py
├── modelo.py
└── main.py
```

---

# Arquivos

| Arquivo | Função |
|----------|--------|
| `descritores.py` | Calcula descritores moleculares utilizando o RDKit a partir de um SMILES. É utilizado pelos demais módulos do projeto. |
| `dados.py` | Carrega e normaliza os dados experimentais da planilha `Dados_QSAR_Saile_PJC2026.xlsx`, além de conter a implementação inicial do descritor de mistura (Fase 3). |
| `ecotox.py` | Carrega e filtra a base pública ECOTOX, obtém os SMILES via PubChem e monta a matriz de treinamento para o modelo QSAR. |
| `modelo.py` | Implementa o ajuste das curvas dose–resposta (EC50), o modelo exploratório com dados próprios, o treinamento utilizando a base ECOTOX e a validação externa. |
| `main.py` | Ponto de entrada da aplicação. Executa os módulos conforme o modo escolhido pelo usuário. |

---

# Como Executar

Coloque a planilha

```text
Dados_QSAR_Saile_PJC2026.xlsx
```

na mesma pasta dos scripts.

---

## 1. Utilizando apenas os dados experimentais

Esse modo **não depende da base ECOTOX**.

Execute:

```bash
python3 main.py dados
```

O pipeline realiza:

- carregamento da planilha experimental;
- normalização dos dados;
- cálculo dos descritores moleculares dos ingredientes com SMILES válido;
- impressão dos dados processados no terminal.

Esse modo é útil para verificar os dados experimentais antes da etapa de modelagem.

---

## 2. Treinamento na base pública + validação externa

Antes da execução, é necessário obter a base pública ECOTOX.

### Passo 1 — Download

Acesse:

> https://cfpub.epa.gov/ecotox/

Selecione:

> **Download the entire database as ASCII files**

---

### Passo 2 — Extração

Extraia o conteúdo do arquivo ZIP em um diretório local, por exemplo:

```text
./ecotox_ascii/
```

A estrutura esperada é semelhante a:

```text
ecotox_ascii/
├── tests.txt
├── results.txt
└── species.txt
```

Os nomes dos arquivos podem variar ligeiramente conforme a versão disponibilizada pela EPA.

---

### Passo 3 — Configuração

No arquivo `ecotox.py`, ajuste a constante:

```python
ECOTOX_DIR = "./ecotox_ascii"
```

caso os arquivos tenham sido extraídos para outro diretório.

---

### Passo 4 — Execução

```bash
python3 main.py publico
```

---

# Fluxo Executado

Durante a execução do modo **publico**, o pipeline realiza automaticamente as seguintes etapas:

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
Descritores Moleculares
      │
      ▼
Treinamento do modelo QSAR
      │
      ▼
Validação cruzada
      │
      ▼
Modelo treinado
      │
      ▼
Predição para os ativos amazônicos
```

Ao final da execução são apresentados:

- coeficiente de determinação (**R²**);
- erro quadrático médio (**RMSE**);
- predição de **pEC50** para cada ingrediente presente na planilha experimental.

---

# Treinamento e Configuração do Modelo QSAR

O modelo preditivo de ecotoxicidade é treinado e configurado utilizando as seguintes definições e etapas:

1. **Algoritmo de Aprendizado:**
   * Utiliza o **Random Forest Regressor** (`RandomForestRegressor` da biblioteca `scikit-learn`).
   * Configurado com `n_estimators=300` (300 árvores de decisão) para obter estimativas estáveis de predição.
   * Executa em paralelo utilizando todas as CPUs disponíveis (`n_jobs=-1`) para otimizar o tempo de processamento.
   * Semente de aleatoriedade fixa (`random_state=42`) para garantir a reprodutibilidade dos resultados.

2. **Conjunto de Entrada (Features):**
   * Descritores moleculares gerados a partir do SMILES via **RDKit**: `MolWt` (Peso Molecular), `LogP` (lipofilicidade), `TPSA` (área polar), `NumHDonors`/`NumHAcceptors` (ligações de hidrogênio), `NumRotatableBonds` (flexibilidade molecular), `RingCount` e `AromaticRings` (quantidade de anéis).

3. **Variável Alvo (Target):**
   * Escala de pEC50, definida pela transformação:
     $$pEC50 = -\log_{10}(EC50 \text{ molar em mol/L})$$
     Isso padroniza os efeitos de toxicidade de compostos químicos independentemente de suas massas molares.

4. **Validação Cruzada (Cross-Validation):**
   * Validação cruzada do tipo **5-Fold Standard Cross-Validation** (divisão dos dados em 5 partes, onde o modelo é treinado em 4 e testado na parte restante, rotacionando até que todos os dados tenham sido testados).
   * As métricas de avaliação calculadas são o Coeficiente de Determinação ($R^2$) e o Erro Quadrático Médio (RMSE) de validação cruzada.

5. **Modelo Final para Produção:**
   * Após a avaliação na validação cruzada, o modelo de regressão final é treinado com **100% dos dados filtrados da base ECOTOX**, sendo então congelado para servir de preditor externo dos novos compostos amazônicos.

## Por que usar Random Forest?

A escolha do **Random Forest** como modelo principal de regressão para este pipeline QSAR baseia-se em fatores técnicos e científicos essenciais:

*   **Robustez a ruídos e outliers:** Dados biológicos de ecotoxicidade pública (ECOTOX) são provenientes de múltiplos estudos científicos independentes, contendo inevitavelmente ruído experimental e outliers. O Random Forest é intrinsecamente resiliente a esses fatores por realizar a média das predições de um conjunto (*ensemble*) de árvores de decisão.
*   **Captura de não linearidades e interações:** A resposta ecotoxicológica (pEC50) não possui relação puramente linear com os descritores químicos. O algoritmo consegue capturar relações não lineares de alta ordem e interações complexas entre múltiplos descritores de forma nativa, sem necessidade de transformação prévia de variáveis.
*   **Resistência à multicolinearidade:** Os descritores químicos calculados (como peso molecular, LogP e área polar) tendem a ser muito correlacionados. O Random Forest seleciona aleatoriamente subconjuntos de atributos em cada nó da árvore de decisão, mitigando os efeitos negativos da correlação entre as variáveis.
*   **Prevenção contra overfitting:** A técnica de *bootstrap aggregating* (bagging) aliada à seleção aleatória de atributos diminui a variância geral do modelo. Isso evita a memorização do conjunto de treinamento e propicia alta capacidade de generalização para os novos compostos (validação externa).
*   **Independência de escalas:** Diferente de algoritmos como regressão linear regularizada, redes neurais ou SVM, as árvores de decisão não são afetadas pela escala dos dados. Isso dispensa etapas complexas de padronização/normalização de descritores com grandezas muito diferentes (ex. Peso Molecular versus LogP).
*   **Importância das variáveis:** O modelo permite extrair a importância relativa de cada descritor molecular para as previsões de ecotoxicidade, viabilizando uma análise interpretável e validação química mecanística do modelo.

---

# Consultas ao PubChem

Durante a construção da base de treinamento, o módulo `ecotox.py` realiza consultas automáticas à API **PubChem PUG REST** para converter números CAS em estruturas SMILES.

Cada composto é consultado apenas uma vez. Como a base filtrada para organismos do tipo alga normalmente contém poucas centenas de compostos, a execução costuma ser rápida.

Mesmo assim, recomenda-se evitar múltiplas execuções consecutivas para não gerar requisições desnecessárias ao serviço.

---

# Problemas de Dados Conhecidos

## SMILES incorreto do Decyl Glucoside

Na aba **`5_SMILES_Ingredientes`** da planilha experimental, o composto **Decyl glucoside** possui um SMILES incorreto.

O registro atual corresponde a um **sal de nióbio**, obtido a partir de uma página *Substance* do PubChem em vez da página *Compound*.

Como o RDKit interpreta esse SMILES sem gerar erro, essa inconsistência não é detectada automaticamente.

É necessário substituir manualmente o SMILES antes da execução do pipeline.

---

## Ausência da tabela de composição das formulações

Ainda não existe uma tabela relacionando cada formulação cosmética aos ingredientes que a compõem e às respectivas concentrações.

Sem essa informação não é possível calcular descritores moleculares de misturas nem validar o modelo para formulações completas.

Atualmente, a validação externa é realizada apenas para **ingredientes isolados**.

A estrutura esperada para a próxima etapa do projeto é semelhante a:

| amostra_id | ingrediente | percentual |
|------------|-------------|-----------:|
| A1 | Andiroba | 20 |
| A1 | Copaíba | 35 |
| A1 | Decyl glucoside | 45 |

Essa informação permitirá implementar a **Fase 3**, descrita em **`architecture.md`**, baseada em descritores moleculares ponderados pela composição das formulações.