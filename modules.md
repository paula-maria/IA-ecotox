# Módulos do Sistema e Funcionamento

Este documento descreve detalhadamente o papel de cada módulo Python presente no repositório, suas principais funções e como eles interagem ao longo do pipeline QSAR de ecotoxicidade.

---

## 1. Módulo Principal (`main.py`)
O [main.py](file:///home/paula/Downloads/IA-ecotox/main.py) serve como o **ponto de entrada e orquestrador** de toda a aplicação. Ele expõe uma interface de linha de comando para invocar os fluxos de processamento baseados no argumento passado.

*   **Principais Funções:**
    *   `rodar_dados()`: Carrega e exibe os dados experimentais dos ativos amazônicos e seus descritores calculados.
    *   `rodar_publico()`: Executa a Fase 1 (Treino do modelo usando a base pública ECOTOX) e a Fase 2 (Predição externa para os ativos amazônicos).
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

---

## 5. Módulo do Modelo QSAR (`modelo.py`)
O [modelo.py](file:///home/paula/Downloads/IA-ecotox/modelo.py) define o algoritmo de aprendizado de máquina para modelar a relação estrutura-atividade quantitativa.

*   **Principais Funções:**
    *   `treinar_modelo_publico()`: Treina o algoritmo **Random Forest Regressor** (configurado com 300 árvores de decisão e otimizado em múltiplos núcleos com `n_jobs=-1`) utilizando toda a matriz de dados da base pública ECOTOX. Reporta o desempenho estimado por validação cruzada 5-fold (R² e RMSE).
    *   `importancia_descritores()`: Extrai a importância relativa de cada variável molecular calculada pelo RDKit para a tomada de decisão do modelo.
    *   `validar_externamente()`: Executa o modelo QSAR previamente congelado sobre os ingredientes amazônicos extraídos pelo módulo `dados.py`, gerando predições de pEC50 e convertendo de volta para concentração equivalente em miligramas por litro (mg/L).

---

## 6. Módulo de Validação e Qualidade (`validacao.py`)
O [validacao.py](file:///home/paula/Downloads/IA-ecotox/validacao.py) audita o pipeline de acordo com as melhores práticas científicas recomendadas pelos Princípios OECD para validação de modelos QSAR.

*   **Principais Funções:**
    *   `checar_qualidade_dados()`: Inspeciona a base de treinamento em busca de duplicatas, valores nulos, outliers estatísticos na variável dependente e descritores sem variabilidade estrutural.
    *   `teste_y_scrambling()`: Embaralha aleatoriamente o vetor de respostas $y$ (`pEC50`) 30 vezes e retreina a Random Forest em paralelo (`n_jobs=-1`). Se o desempenho com o vetor real não for significativamente superior ao das repetições com alvos embaralhados, sinaliza risco de falso positivo estatístico (overfitting de ruído).
    *   `calcular_leverage(X_treino, X_novo)`: Determina se os novos compostos (ativos amazônicos) pertencem ao espaço químico do modelo (Domínio de Aplicabilidade) calculando a distância de leverage ($h$) e comparando-a com o limite crítico usual ($h^* = 3(p+1)/n$).
    *   `matriz_confusao_toxicidade()`: Transforma os valores de toxicidade contínuos em categorias discretas regulamentadas pelo **GHS (Sistema Globalmente Harmonizado)** da ONU para predições e valores observados.
    *   `plotar_matriz_confusao()`: Gera e salva um Heatmap gráfico (`matriz_confusao.png`) demonstrando os acertos e erros do modelo sob validação cruzada.
