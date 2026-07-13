# QSAR Ecotoxicidade — Ativos Amazônicos

Pipeline em **Python** para prever a **ecotoxicidade** (inibição do crescimento de algas) de ativos cosméticos amazônicos utilizando modelos **QSAR (Quantitative Structure–Activity Relationship)** baseados em descritores moleculares.

O projeto está organizado em módulos independentes para facilitar a manutenção, reutilização do código e evolução para futuras etapas de modelagem de misturas. Para compreender a arquitetura completa e o fluxo das três fases do projeto, consulte **`architecture.md`**.

---

# Instalação

Requer **Python 3.10** ou superior.

## Ambiente virtual (recomendado)

```bash
python3 -m venv venv

# Linux/macOS
source venv/bin/activate

# Windows
venv\Scripts\activate

pip install -r requirements.txt
```

## Instalação direta

```bash
pip install -r requirements.txt --break-system-packages
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
├── validacao.py
└── main.py
```

---

# Arquivos

| Arquivo | Função |
|----------|--------|
| `descritores.py` | Calcula descritores moleculares utilizando o RDKit a partir de um SMILES. É utilizado pelos demais módulos do projeto. |
| `dados.py` | Carrega e normaliza os dados experimentais da planilha `Dados_QSAR_Saile_PJC2026.xlsx`. |
| `ecotox.py` | Carrega e filtra a base pública ECOTOX, obtém os SMILES via PubChem e monta a matriz de treinamento para o modelo QSAR. |
| `modelo.py` | Implementa o treinamento do modelo Random Forest com a base ECOTOX e a validação externa com os dados amazônicos. |
| `validacao.py` | Avaliação da qualidade dos dados, Y-scrambling, domínio de aplicabilidade (leverage) e matriz de confusão das categorias GHS. |
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

---

## 2. Treinamento na base pública + validação externa

Antes da execução, é necessário obter a base pública ECOTOX (https://cfpub.epa.gov/ecotox/), extrair o conteúdo (arquivos ASCII) em `./ecotox_ascii/` e configurar o caminho em `ecotox.py`.

### Execução

```bash
python3 main.py publico
```

---

## 3. Validação completa do modelo

```bash
python3 main.py validar
```

Executa:

- verificação de qualidade dos dados (duplicatas, outliers, descritores sem variação);
- **Y-scrambling** — testa se o modelo aprendeu sinal químico real ou apenas ruído;
- **domínio de aplicabilidade** (leverage) — verifica se os ativos amazônicos pertencem ao espaço químico conhecido pelo modelo;
- **matriz de confusão GHS** — avalia o desempenho por categoria de toxicidade aguda aquática.

---

# Fluxo Executado

```text
ECOTOX
      │
      ▼
Filtragem dos testes (Chlorella vulgaris)
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
Treinamento do modelo QSAR (Random Forest)
      │
      ▼
Validação cruzada (5-fold)
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

# Modelo QSAR

## Algoritmo

O modelo preditivo de ecotoxicidade é treinado com as seguintes definições:

| Configuração | Valor |
|---|---|
| Algoritmo | Random Forest Regressor (`scikit-learn`) |
| Número de árvores | 300 (`n_estimators=300`) |
| Paralelismo | todas as CPUs disponíveis (`n_jobs=-1`) |
| Reprodutibilidade | semente fixa (`random_state=42`) |
| Validação | 5-Fold Cross-Validation |

## Variáveis de entrada (descritores moleculares)

Calculados via **RDKit** a partir do SMILES de cada composto:

| Descritor | Significado |
|---|---|
| `MolWt` | Peso molecular |
| `LogP` | Lipofilicidade |
| `TPSA` | Área polar total |
| `NumHDonors` | Doadores de ligação de hidrogênio |
| `NumHAcceptors` | Aceptores de ligação de hidrogênio |
| `NumRotatableBonds` | Flexibilidade molecular |
| `RingCount` | Número de anéis |
| `AromaticRings` | Número de anéis aromáticos |

## Variável-resposta

$$pEC50 = -\log_{10}(EC50 \text{ em mol/L})$$

Padroniza os efeitos de toxicidade independentemente da massa molar dos compostos.

## Por que usar Random Forest?

A escolha do **Random Forest** como modelo principal de regressão para este pipeline QSAR baseia-se em fatores técnicos e científicos:

- **Robustez a ruídos e outliers:** dados biológicos de ecotoxicidade pública (ECOTOX) são provenientes de múltiplos estudos independentes, contendo inevitavelmente ruído experimental. O Random Forest é resiliente a esses fatores por realizar a média de um conjunto (*ensemble*) de árvores de decisão.
- **Captura de não linearidades:** a resposta ecotoxicológica (pEC50) não possui relação linear com os descritores químicos. O algoritmo captura relações complexas de forma nativa.
- **Resistência à multicolinearidade:** descritores químicos tendem a ser correlacionados. O Random Forest seleciona aleatoriamente subconjuntos de atributos em cada nó, mitigando esse efeito.
- **Prevenção contra overfitting:** a técnica de *bagging* aliada à seleção aleatória de atributos reduz a variância geral do modelo.
- **Independência de escalas:** árvores de decisão não são afetadas pela escala dos dados — dispensando normalização prévia dos descritores.
- **Importância das variáveis:** permite extrair a importância relativa de cada descritor, viabilizando interpretação química dos resultados.

### Respaldo científico

A escolha metodológica é apoiada por **Ferro (2025)** — *Aprendizado de Máquina para predição de crescimento de Chlorella vulgaris* (UFAL, 2025) —, que comparou Random Forest, SVR e redes neurais para a mesma espécie (*C. vulgaris*) e identificou o Random Forest como o algoritmo de melhor desempenho, com validação cruzada.

A diferença fundamental que preserva o ineditismo deste projeto:

| Aspecto | Ferro (UFAL, 2025) | Este projeto (PJC 2026) |
|---|---|---|
| **Objetivo** | Predizer crescimento/produção da alga | Predizer toxicidade de ingredientes sobre a alga |
| **Variáveis de entrada** | Nutrientes, intensidade luminosa, tempo | Descritores moleculares QSAR |
| **Contexto** | Cultivo biotecnológico | Avaliação ecotoxicológica |
| **Aplicação** | Biorefinaria | Segurança cosmética ambiental |

São abordagens complementares: a de Ferro valida metodologicamente a escolha do algoritmo; a deste projeto aplica essa metodologia a uma pergunta científica inédita.

---

# Consultas ao PubChem

Durante a construção da base de treinamento, o módulo `ecotox.py` realiza consultas automáticas à API **PubChem PUG REST** para converter números CAS em estruturas SMILES.

As respostas são armazenadas em cache local (`pubchem_cache.json`). Cada composto é consultado apenas uma vez, acelerando execuções subsequentes.

---

# Problemas de Dados Conhecidos

## SMILES incorreto do Decyl Glucoside

Na aba **`5_SMILES_Ingredientes`** da planilha experimental, o composto **Decyl glucoside** possui um SMILES incorreto. O registro atual corresponde a um **sal de nióbio**. É necessário substituir manualmente o SMILES antes da execução do pipeline.

## Ausência da tabela de composição das formulações

Ainda não existe uma tabela relacionando cada formulação cosmética aos ingredientes que a compõem e às respectivas concentrações. Atualmente, a validação externa é realizada apenas para **ingredientes isolados**.

A estrutura esperada para a próxima etapa do projeto é semelhante a:

| amostra_id | ingrediente | percentual |
|------------|-------------|-----------:|
| A1 | Andiroba | 20 |
| A1 | Copaíba | 35 |
| A1 | Decyl glucoside | 45 |

Essa informação permitirá implementar a **Fase 3**, descrita em **`architecture.md`**, baseada em descritores moleculares ponderados pela composição das formulações.