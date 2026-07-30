# Módulos do Sistema e Funcionamento

Este documento descreve detalhadamente o papel de cada módulo Python presente no repositório, suas principais funções e como eles interagem ao longo do pipeline QSAR de ecotoxicidade.

---

## 1. Módulo Principal (`main.py`)
O [main.py](file:///home/paula/Downloads/IA-ecotox/main.py) serve como o **ponto de entrada e orquestrador** de toda a aplicação. Ele expõe uma interface de linha de comando para invocar os fluxos de processamento baseados no argumento passado.

*   **Principais Funções:**
    *   `rodar_dados()`: Carrega e exibe os dados experimentais dos ativos amazônicos e seus descritores calculados.
    *   `rodar_publico()`: Executa o pipeline com a base **ECOTOX** (Fase 1 — treino, Fase 2 — predição). Comando: `python3 main.py publico`.
    *   `rodar_publico_combinado()`: Executa o mesmo pipeline, mas treinando com **ECOTOX + EnviroTox combinados** (~2.140 compostos). Comando: `python3 main.py publico-combinado`.
    *   `rodar_validacao()`: Executa a auditoria completa de validação estatística do modelo público (avaliação de qualidade de dados, teste de Y-scrambling, cálculo do Domínio de Aplicabilidade e geração da Matriz de Confusão por categoria GHS).

---

## 2. Módulo de Descritores Moleculares (`descritores.py`)
O [descritores.py](file:///home/paula/Downloads/IA-ecotox/descritores.py) é a interface com a biblioteca **RDKit**. Ele recebe uma estrutura molecular textual em formato **SMILES** e computa propriedades físico-químicas e estruturais da molécula.

*   **Principais Funções:**
    *   `calcular_descritores(smiles: str)`: Tenta converter o SMILES em um objeto molecular do RDKit (`Chem.MolFromSmiles`). Se bem-sucedido, calcula um conjunto de descritores padrão:
        *   `MolWt` (Peso Molecular)
        *   `LogP` (Coeficiente de partição Octanol-Água)
        *   `TPSA` (Área de superfície polar topológica)
        *   `NumHDonors` / `NumHAcceptors` (Doadores/Aceitadores de Hidrogênio)
        *   `NumRotatableBonds` (Ligações Rotacionáveis)
        *   `RingCount` (Quantidade total de anéis)
        *   `AromaticRings` (Quantidade de anéis aromáticos)
        *   `MACCS Keys` (166 features binárias com perfis estruturais de grupamentos funcionais)

---

## 3. Módulo de Dados Experimentais (`dados.py`)
O [dados.py](file:///home/paula/Downloads/IA-ecotox/dados.py) gerencia a leitura e normalização dos dados obtidos em ensaios de laboratório contidos na planilha Excel `Dados_QSAR_Saile_PJC2026.xlsx`.

*   **Principais Funções:**
    *   `carregar_dados_inibicao()`: Lê e limpa a tabela de ensaios experimentais (ex.: normalizando strings, tratando porcentagens de inibição).
    *   `carregar_descritores_ingredientes()`: Extrai os SMILES dos ingredientes da planilha e aciona o módulo `descritores.py` para gerar as variáveis independentes necessárias para predição.
    *   `descritor_mistura_ponderado()`: Implementação inicial para a futura Fase 3 (descritores moleculares ponderados baseados na composição percentual de formulações cosméticas completas).

---

## 4. Módulo da Base Pública ECOTOX (`ecotox.py`)
O [ecotox.py](file:///home/paula/Downloads/IA-ecotox/ecotox.py) constrói a base de treinamento a partir dos dados do US EPA ECOTOX. Ele filtra apenas ensaios aplicáveis e obtém estruturas SMILES.

*   **Principais Funções:**
    *   `carregar_ecotox_algas()`: Lê as tabelas brutas `species.txt`, `tests.txt` e `results.txt` sob demanda (usando filtros de colunas otimizados para economia de memória) e filtra ensaios de crescimento realizados com a espécie alvo (*Chlorella vulgaris*) ou com os gêneros de alga correlacionados (fallback).
    *   `cas_para_smiles(cas: str, cache: dict)`: Consulta o número CAS no serviço REST do PubChem para retornar o SMILES canônico correspondente, utilizando controle de taxa de requisições.
    *   `anexar_smiles(df: pd.DataFrame)`: Traduz os CAS da tabela para SMILES, utilizando um sistema de **cache local em disco** (`pubchem_cache.json`) que elimina requisições repetidas na web e exibe o progresso em tempo real.
    *   `montar_matriz_treino()`: Converte as concentrações de efeito (ex.: EC50 em mg/L) para escala molar baseada no peso molecular do RDKit e gera a matriz de treino contendo descritores moleculares estruturados e a variável-alvo `pEC50` ($-\log_{10}(\text{EC50}_{\text{molar}})$).
    *   `_conc_para_mol_por_l()` / `_UNIT_TO_G_PER_L` / `_UNIT_TO_MOL_PER_L`: Lógica de conversão de unidades de concentração (mg/L, μg/L, mmol/L, etc.) reutilizada também pelo módulo `envirotox.py`.

---

## 5. Módulo da Base Pública EnviroTox (`envirotox.py`)
O [envirotox.py](file:///home/paula/Downloads/IA-ecotox/envirotox.py) carrega e pré-processa os dados do banco **EnviroTox** (envirotoxdatabase.org), espelhando a estrutura de `ecotox.py`. Diferencia-se por usar os arquivos Excel (com sheets `test` e `substance`) e por aproveitar os SMILES já presentes no banco, sem chamar o PubChem desnecessariamente.

*   **Conjuntos carregados:**
    *   `algae/` — filtro amplo (Trophic Level = ALGAE): 11.102 linhas, 186 espécies.
    *   `chlorella/` — filtro específico (*Chlorella vulgaris*): 542 linhas.
    *   Os dois conjuntos são **combinados**, dando preferência à espécie mais específica quando um CAS aparece nos dois.

*   **Principais Funções:**
    *   `carregar_envirotox()`: Lê as sheets `test` e `substance` dos Excels, filtra endpoints em `{EC50, IC50, NOEC, LOEC}`, faz merge por CAS, usa o SMILES direto da coluna `Desalted Canonical SMILES` e chama PubChem apenas como fallback para CAS sem SMILES (~8% dos casos).
    *   `montar_matriz_envirotox(df)`: Constrói a matriz de treino com os mesmos descritores RDKit do ECOTOX, aplica a conversão de unidades (`_conc_para_mol_por_l` de `ecotox.py`), deduplica por CAS (mesmo critério std pEC50 ≤ 1.0) e **mantém a coluna `latin_name`** no resultado para uso como feature categórica.

---

## 6. Módulo de Combinação de Fontes (`fontes_externas.py`)
O [fontes_externas.py](file:///home/paula/Downloads/IA-ecotox/fontes_externas.py) une as matrizes do ECOTOX e do EnviroTox em uma única base de treino ampliada, com controle de qualidade entre fontes.

*   **Principais Funções:**
    *   `combinar_fontes(matriz_ecotox, matriz_envirotox)`: Alinha as colunas de descritores (interseção das numéricas), identifica CAS compartilhados e exclusivos de cada fonte. Para CAS presentes nas **duas** fontes, calcula o desvio padrão do pEC50 entre elas: se `std > 1.0` (mesma regra de `ecotox.montar_matriz_treino`), descarta e reporta; se concordam, tira a média. Gera colunas one-hot `especie_*` a partir de `latin_name` para uso como features no modelo. Reporta detalhadamente: CAS exclusivos, compartilhados, descartados por divergência.

*   **Resultado da combinação atual:**

    | Métrica | Valor |
    |---|---|
    | Compostos ECOTOX | 1.720 CAS |
    | Compostos EnviroTox | 1.220 CAS |
    | **CAS exclusivos do EnviroTox (novos)** | **447** |
    | CAS compartilhados (média aplicada) | 746 |
    | CAS descartados (divergência inter-fonte) | 27 |
    | **Total combinado** | **2.140 compostos** |

---

## 7. Módulo do Modelo QSAR (`modelo.py`)
O [modelo.py](file:///home/paula/Downloads/IA-ecotox/modelo.py) define o algoritmo de aprendizado de máquina para modelar a relação estrutura-atividade quantitativa.

*   **Principais Funções:**
    *   `treinar_modelo_publico()`: Executa um pipeline com `StandardScaler` e aplica o **GridSearchCV** para encontrar o melhor algoritmo entre o **Random Forest Regressor** e o **SVR**, utilizando `n_jobs=2` para não travar o PC do usuário. Reporta o desempenho estimado por validação cruzada 5-fold (R² e RMSE) e elege o melhor modelo.
    *   `importancia_descritores()`: Extrai a importância relativa de cada variável molecular calculada pelo RDKit usando a técnica **Permutation Importance**, que é agnóstica ao modelo vencedor.
    *   `validar_externamente()`: Executa o modelo QSAR previamente congelado sobre os ingredientes amazônicos extraídos pelo módulo `dados.py`, gerando predições de pEC50 e convertendo de volta para concentração equivalente em miligramas por litro (mg/L).

---

## 8. Módulo de Validação e Qualidade (`validacao.py`)
O [validacao.py](file:///home/paula/Downloads/IA-ecotox/validacao.py) audita o pipeline de acordo com as melhores práticas científicas recomendadas pelos Princípios OECD para validação de modelos QSAR.

*   **Principais Funções:**
    *   `checar_qualidade_dados()`: Inspeciona a base de treinamento em busca de duplicatas, valores nulos, outliers estatísticos na variável dependente e descritores sem variabilidade estrutural.
    *   `teste_y_scrambling()`: Embaralha aleatoriamente o vetor de respostas $y$ (`pEC50`) 30 vezes e retreina a Random Forest (`n_jobs=2`). Se o desempenho com o vetor real não for significativamente superior ao das repetições com alvos embaralhados, sinaliza risco de falso positivo estatístico (overfitting de ruído).
    *   `calcular_leverage(X_treino, X_novo)`: Determina se os novos compostos (ativos amazônicos) pertencem ao espaço químico do modelo (Domínio de Aplicabilidade) calculando a distância de leverage ($h$) e comparando-a com o limite crítico usual ($h^* = 3(p+1)/n$).
    *   `matriz_confusao_toxicidade()`: Transforma os valores de toxicidade contínuos em categorias discretas regulamentadas pelo **GHS (Sistema Globalmente Harmonizado)** da ONU para predições e valores observados.
    *   `plotar_matriz_confusao()`: Gera e salva um Heatmap gráfico (`matriz_confusao.png`) demonstrando os acertos e erros do modelo sob validação cruzada.
