import os
import streamlit as st
import pandas as pd

import dados
import ecotox
import modelo
import validacao

st.set_page_config(page_title="QSAR Ecotoxicidade", layout="wide")

# ---- Injeção de CSS Customizado ----
st.markdown("""
<style>
    /* Arredondar botões e adicionar hover */
    .stButton>button {
        border-radius: 8px !important;
        transition: all 0.3s ease;
        border: none !important;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(46, 125, 50, 0.2);
    }
    
    /* Estilizar cartões de DataFrame */
    [data-testid="stDataFrame"] {
        border-radius: 10px;
        overflow: hidden;
        box-shadow: 0 4px 15px rgba(0,0,0,0.03);
    }
    
    /* Remover borda do file uploader e deixar mais sutil */
    [data-testid="stFileUploader"] {
        background-color: #ffffff;
        border-radius: 10px;
        padding: 10px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.02);
    }
</style>
""", unsafe_allow_html=True)

st.title("Pipeline QSAR Ecotoxicidade")
st.markdown("Interface para o pipeline de aprendizado de máquina para previsão de toxicidade de ativos amazônicos.")

# ---- Cache do Modelo ----
@st.cache_resource
def load_public_data():
    """Carrega e prepara a base pública apenas uma vez."""
    arquivo_matriz = "matriz_treino_ecotox.csv"
    if os.path.exists(arquivo_matriz):
        # Se a matriz já foi processada (ideal para deploy), carrega direto
        return pd.read_csv(arquivo_matriz)
        
    # Senão, processa a partir da base bruta (ecotox_ascii)
    dados_ecotox = ecotox.carregar_ecotox_algas()
    dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
    matriz = ecotox.montar_matriz_treino(dados_ecotox)
    return matriz

# ---- Menu Lateral ----
modo = st.sidebar.radio(
    "Escolha a etapa de execução:",
    [
        "0. Como Funciona / Tutorial",
        "1. Dados Experimentais", 
        "2. Treinamento Público + Previsão", 
        "3. Validação do Modelo"
    ]
)

caminho_planilha = dados.ARQUIVO_PADRAO
st.sidebar.markdown(f"**Planilha padrão:** `{caminho_planilha}`")
upload = st.sidebar.file_uploader("Ou envie outra planilha (Opcional):", type=["xlsx"])

if upload is not None:
    caminho_planilha = "temp_upload.xlsx"
    with open(caminho_planilha, "wb") as f:
        f.write(upload.getbuffer())

# ==========================================
# MODO 0: COMO FUNCIONA / TUTORIAL
# ==========================================
if modo == "0. Como Funciona / Tutorial":
    st.header("Sobre o Projeto e Tutorial")
    st.info("Abaixo você encontra a documentação oficial do pipeline, explicando o algoritmo, a base de dados e a base científica das previsões.")
    
    st.markdown("""
### Como usar este aplicativo
Utilize o **Menu Lateral** para navegar pelas etapas do pipeline:

1. **Dados Experimentais**: 
   - *O que faz:* Normaliza sua planilha e calcula descritores dos seus ativos (sem depender da base ECOTOX).
   - *Como usar:* Caso tenha modificado a planilha de dados, faça o upload no menu lateral e clique em **Executar Etapa 1**. Ele mostrará a tabela de amostras e gerará os gráficos de curvas de crescimento.

2. **Treinamento Público + Previsão**: 
   - *O que faz:* Treina o algoritmo Random Forest na base pública e gera as previsões (pEC50) finais.
   - *Como usar:* Clique em **Executar Etapa 2**. Se você fez upload de uma planilha, ele usará os SMILES dela. Ao final, a tabela com o pEC50 previsto de todos os ingredientes aparecerá na tela e você poderá baixar como `.csv`.

3. **Validação do Modelo**: 
   - *O que faz:* Prova a robustez e confiabilidade das previsões feitas.
   - *Como usar:* Clique em **Executar Etapa 3**. Ele verificará se os dados dos ativos estão dentro do Domínio de Aplicabilidade (leverage), fará o teste de ruído (Y-Scrambling) e gerará a Matriz de Confusão.

---
### Documentação Técnica
""")
    st.page_link("https://github.com/paula-maria/IA-ecotox", label="Repositório Oficial no GitHub", icon="🔗")
    
    try:
        with open("readme.md", "r", encoding="utf-8") as f:
            readme_content = f.read()
            st.markdown(readme_content)
    except FileNotFoundError:
        st.warning("O arquivo readme.md não foi encontrado nesta pasta.")

# ==========================================
# MODO 1: DADOS
# ==========================================
if modo == "1. Dados Experimentais":
    st.header("1. Processamento dos Dados Experimentais")
    st.markdown("Apenas normaliza a planilha própria e calcula os descritores, sem usar a base ECOTOX.")
    
    if st.button("Executar Etapa 1", type="primary"):
        with st.spinner("Processando dados e descritores..."):
            try:
                df = dados.carregar_dados_inibicao(caminho_planilha)
                st.subheader("Dados de Inibição (Amostra)")
                st.dataframe(df.head(), use_container_width=True)
                
                st.info("Gerando relatórios e gráficos das curvas de crescimento...")
                dados.gerar_relatorios_e_graficos(df)
                st.success("Arquivos Excel e PNGs de curvas gerados na pasta do projeto.")
                
                # Exibir os gráficos gerados no Streamlit
                if os.path.exists("assets/curvas_crescimento_tcc_original.png"):
                    st.image("assets/curvas_crescimento_tcc_original.png", caption="Curvas de Crescimento - TCC", use_container_width=True)
                if os.path.exists("assets/curvas_crescimento_novo_exp.png"):
                    st.image("assets/curvas_crescimento_novo_exp.png", caption="Curvas de Crescimento - Novo Experimento", use_container_width=True)
                    
                # Exibir botão de download para o arquivo Excel
                if os.path.exists("relatorio_contagem_celular.xlsx"):
                    with open("relatorio_contagem_celular.xlsx", "rb") as file:
                        st.download_button(
                            label="Baixar Dados",
                            data=file,
                            file_name="relatorio_contagem_celular.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )

                
                st.subheader("Descritores Moleculares Calculados")
                descritores_df = dados.carregar_descritores_ingredientes(caminho_planilha)
                st.dataframe(descritores_df, use_container_width=True)
            except Exception as e:
                st.error(f"Erro ao processar dados: {e}")

# ==========================================
# MODO 2: PÚBLICO
# ==========================================
elif modo == "2. Treinamento Público + Previsão":
    st.header("2. Treinamento (ECOTOX) e Previsão")
    st.markdown("Treina o modelo na base pública e realiza a validação externa (previsão) para os ativos da planilha.")
    
    if st.button("Executar Etapa 2", type="primary"):
        with st.spinner("Carregando base pública ECOTOX e treinando o modelo (Isso pode demorar na 1ª vez)..."):
            try:
                matriz = load_public_data()
                
                # Para suprimir os prints originais do `treinar_modelo_publico` e colocar no Streamlit,
                # chamamos os métodos, mas o `treinar_modelo_publico` faz prints internos de RMSE.
                modelo_treinado, colunas_x, scaler = modelo.treinar_modelo_publico(matriz)
                
                st.subheader("Importância dos Descritores")
                ranking = modelo.importancia_descritores(modelo_treinado, colunas_x, scaler, matriz)
                st.dataframe(ranking, use_container_width=True)
                
                st.subheader("Previsão para os ativos amazônicos (pEC50)")
                previsao = modelo.validar_externamente(modelo_treinado, colunas_x, caminho_planilha, scaler)
                st.dataframe(previsao, use_container_width=True)
                
                # Botão de download
                csv = previsao.to_csv(index=False).encode('utf-8')
                st.download_button("Baixar Previsões (CSV)", csv, "previsoes_ativos.csv", "text/csv")
                
            except Exception as e:
                st.error(f"Erro na execução pública: {e}")

# ==========================================
# MODO 3: VALIDAÇÃO
# ==========================================
elif modo == "3. Validação do Modelo":
    st.header("3. Validação Completa")
    st.markdown("Avalia a qualidade dos dados, domínio de aplicabilidade, y-scrambling e matriz de confusão GHS.")
    
    if st.button("Executar Etapa 3", type="primary"):
        with st.spinner("Realizando testes de validação no modelo (Pode demorar alguns minutos)..."):
            try:
                matriz = load_public_data()
                
                st.subheader("Qualidade dos Dados")
                relatorio_qualidade = validacao.checar_qualidade_dados(matriz)
                # O método original imprime, então vamos usar o st.text pro output se não for dicionário puro
                # Mas como `relatorio_qualidade` é dict, podemos renderizar:
                st.json(relatorio_qualidade)
                
                colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
                X, y = matriz[colunas_x], matriz["pEC50"]
                
                st.subheader("Teste de Y-Scrambling")
                resultado_scrambling = validacao.teste_y_scrambling(X, y)
                st.json(resultado_scrambling)
                if not resultado_scrambling["passou_no_teste"]:
                    st.warning("[ATENÇÃO] R² real não se distanciou o suficiente do embaralhado.")
                else:
                    st.success("O modelo aprendeu sinal químico real (passou no Y-scrambling).")
                
                st.subheader("Domínio de Aplicabilidade (Ativos Amazônicos)")
                ingredientes = dados.carregar_descritores_ingredientes(caminho_planilha)
                colunas_comuns = [c for c in colunas_x if c in ingredientes.columns]
                leverage_df = validacao.calcular_leverage(X[colunas_comuns], ingredientes[colunas_comuns])
                leverage_df["Ingrediente"] = ingredientes["Ingrediente"].values
                
                st.write(f"**h* (limite de corte) = {leverage_df.attrs['h_estrela']:.3f}**")
                st.dataframe(leverage_df[["Ingrediente", "leverage", "dentro_do_dominio"]], use_container_width=True)
                
                st.subheader("Matriz de Confusão (Categorias GHS)")
                from sklearn.ensemble import RandomForestRegressor
                from sklearn.model_selection import KFold, cross_val_predict
                
                modelo_cv = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=2)
                kf = KFold(n_splits=5, shuffle=True, random_state=42)
                pred_cv = cross_val_predict(modelo_cv, X, y, cv=kf)
                
                df_matriz = validacao.matriz_confusao_toxicidade(y, pd.Series(pred_cv), matriz["MolWt"])
                st.dataframe(df_matriz, use_container_width=True)
                
                caminho_png = validacao.plotar_matriz_confusao(df_matriz)
                st.image(caminho_png, caption="Matriz de Confusão GHS")
                
            except Exception as e:
                st.error(f"Erro na validação: {e}")
