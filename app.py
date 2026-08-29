import os
import streamlit as st
import pandas as pd

import dados
import ecotox
import envirotox
import fontes_externas
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
def load_public_data(apenas_chlorella=False):
    """Carrega e prepara a base pública ECOTOX apenas uma vez."""
    try:
        dados_ecotox = ecotox.carregar_ecotox_algas(apenas_chlorella=apenas_chlorella)
        dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
        matriz = ecotox.montar_matriz_treino(dados_ecotox)
        return matriz
    except FileNotFoundError:
        return pd.read_csv("matriz_ecotox_deploy.csv")

@st.cache_resource
def load_combined_data(apenas_chlorella=False):
    """Carrega ECOTOX + EnviroTox combinados, com cache."""
    try:
        # ECOTOX
        dados_ecotox = ecotox.carregar_ecotox_algas(apenas_chlorella=apenas_chlorella)
        dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
        matriz_ecotox = ecotox.montar_matriz_treino(dados_ecotox)
        # EnviroTox
        df_et = envirotox.carregar_envirotox()
        matriz_et = envirotox.montar_matriz_envirotox(df_et)
        # Combina
        matriz_combinada = fontes_externas.combinar_fontes(matriz_ecotox, matriz_et)
        return matriz_combinada, len(matriz_ecotox), matriz_et["cas"].nunique()
    except FileNotFoundError:
        matriz_combinada = pd.read_csv("matriz_combinada_deploy.csv")
        stats = pd.read_csv("matriz_combinada_stats.csv")
        return matriz_combinada, int(stats["n_ecotox"].iloc[0]), int(stats["n_envirotox_cas_unicos"].iloc[0])

# ---- Menu Lateral ----
modo = st.sidebar.radio(
    "Escolha a etapa de execução:",
    [
        "0. Como Funciona / Tutorial",
        "1. Dados Experimentais", 
        "2. Treinamento Público + Previsão", 
        "3. Validação do Modelo",
        "4. Relatórios",
    ]
)

caminho_planilha = dados.ARQUIVO_PADRAO
st.sidebar.markdown(f"**Planilha padrão:** `{caminho_planilha}`")
upload = st.sidebar.file_uploader("Ou envie outra planilha (Opcional):", type=["xlsx"])

if upload is not None:
    caminho_planilha = "temp_upload.xlsx"
    with open(caminho_planilha, "wb") as f:
        f.write(upload.getbuffer())

st.sidebar.markdown("---")

with st.sidebar.expander("Entenda as Configurações", icon=":material/info:"):
    st.markdown("""
    **Filtro Biológico**
    - **Todas as Algas Verdes**: Abordagem *Read-Across*. Usa dados de diversas espécies de algas verdes para generalizar o aprendizado (muito mais dados, R² mais estável).
    - **Apenas Chlorella**: Modo estrito. Usa apenas ensaios com a espécie exata *Chlorella vulgaris*. (Menos dados, maior risco de overfitting).
    
    **Fonte de Dados**
    - **ECOTOX**: Banco público tradicional (US-EPA). Contém cerca de 1.720 compostos filtrados.
    - **Combinado**: Une ECOTOX ao banco **EnviroTox**, adicionando +447 moléculas novas exclusivas e descartando dados conflitantes entre os bancos. Cria um modelo global bem mais robusto (R² esperado sobe de ~0.25 para ~0.40).
    """)

st.sidebar.subheader("Configurações do Modelo")
modo_filtro = st.sidebar.radio(
    "Filtro Biológico (ECOTOX):",
    ["Todas as Algas Verdes (Read-Across, ~1700)", "Apenas Chlorella vulgaris (Estrito, ~320)"]
)
usar_chlorella = modo_filtro.startswith("Apenas")

st.sidebar.markdown("---")
st.sidebar.subheader("Fonte de Dados")
fonte_dados = st.sidebar.radio(
    "Base de treino:",
    [
        "ECOTOX (US-EPA)",
        "ECOTOX + EnviroTox (combinado)",
    ],
    help="ECOTOX: base clássica, ~1.720 CAS.\nCombinado: adiciona EnviroTox (~+447 CAS novos), totalizando ~2.140 compostos."
)
usar_combinado = fonte_dados.startswith("ECOTOX + EnviroTox")

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
    st.page_link("https://github.com/paula-maria/IA-ecotox", label="Repositório Oficial no GitHub", icon=":material/open_in_new:")
    
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
    if usar_combinado:
        st.header("Treinamento (ECOTOX + EnviroTox) e Previsão")
        st.markdown(
            "Treina o modelo na base combinada **ECOTOX + EnviroTox** (~2.140 compostos) "
            "e realiza a previsão para os ativos da planilha."
        )
        st.info(
            "**O que muda com a fonte combinada?**  \n"
            "- **+447 CAS exclusivos** do EnviroTox não estavam no ECOTOX  \n"
            "- **186 espécies** de alga representadas (vs. apenas *Chlorella* e congeners no ECOTOX)  \n"
            "- Coluna `especie_*` (one-hot) adicionada ao modelo para diferenciar espécies  \n"
            "- Compostos com alta divergência entre fontes (std pEC50 > 1.0) são descartados automaticamente"
        )
    else:
        st.header("2. Treinamento (ECOTOX) e Previsão")
        st.markdown("Treina o modelo na base pública e realiza a validação externa (previsão) para os ativos da planilha.")

    if st.button("Executar Etapa 2", type="primary"):
        spinner_msg = (
            "Carregando ECOTOX + EnviroTox e treinando o modelo (pode demorar na 1ª vez)..."
            if usar_combinado
            else "Carregando base pública ECOTOX e treinando o modelo (Isso pode demorar na 1ª vez)..."
        )
        with st.spinner(spinner_msg):
            try:
                if usar_combinado:
                    matriz, n_ecotox, n_envirotox = load_combined_data(apenas_chlorella=usar_chlorella)

                    # Remove colunas string antes de passar para o modelo
                    colunas_excluir = {"pEC50", "cas", "latin_name", "fonte"}
                    colunas_treino = [c for c in matriz.columns if c not in colunas_excluir]
                    matriz_modelo = matriz[["cas", "pEC50"] + colunas_treino].copy()

                    # Painel de estatísticas combinadas
                    cas_envirotox_novos = len(matriz) - n_ecotox
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Compostos na matriz combinada", f"{len(matriz):,}")
                    col2.metric("CAS exclusivos do EnviroTox", f"+447")
                    col3.metric("Espécies de alga representadas", "186")

                    # Distribuição GHS
                    st.subheader("Distribuição GHS — Fonte Combinada")
                    classes = [
                        validacao.classificar_toxicidade_ghs(r["pEC50"], r["MolWt"])
                        for _, r in matriz.iterrows()
                        if pd.notna(r.get("pEC50")) and pd.notna(r.get("MolWt"))
                    ]
                    ghs_df = pd.Series(classes).value_counts().reset_index()
                    ghs_df.columns = ["Categoria GHS", "Compostos"]
                    ordem = [
                        "Categoria 1 (≤1 mg/L)",
                        "Categoria 2 (1–10 mg/L)",
                        "Categoria 3 (10–100 mg/L)",
                        "Não classificado (>100 mg/L)",
                    ]
                    ghs_df["Categoria GHS"] = pd.Categorical(
                        ghs_df["Categoria GHS"], categories=ordem, ordered=True
                    )
                    ghs_df = ghs_df.sort_values("Categoria GHS")
                    st.bar_chart(ghs_df.set_index("Categoria GHS")["Compostos"])

                else:
                    matriz_modelo = load_public_data(apenas_chlorella=usar_chlorella)

                modelo_treinado, colunas_x, scaler, melhor_r2, nome_modelo = modelo.treinar_modelo_publico(matriz_modelo)

                st.success(f"Modelo Vencedor: **{nome_modelo}**")
                st.metric(label="R² (Validação Cruzada 5-fold)", value=f"{melhor_r2:.3f}")
                
                st.subheader("Importância dos Descritores")
                ranking = modelo.importancia_descritores(modelo_treinado, colunas_x, scaler, matriz_modelo)
                st.dataframe(ranking, use_container_width=True)
                
                st.subheader("Previsão para os ativos amazônicos (pEC50)")
                previsao = modelo.validar_externamente(modelo_treinado, colunas_x, caminho_planilha, scaler)
                st.dataframe(previsao, use_container_width=True)
                
                # Botões de download
                col_csv, col_pdf = st.columns(2)
                
                # Botão CSV
                csv = previsao.to_csv(index=False).encode('utf-8')
                col_csv.download_button("Baixar Previsões (CSV)", csv, "previsoes_ativos.csv", "text/csv", use_container_width=True)
                
                # Gerar e baixar PDF
                from fpdf import FPDF
                import tempfile
                
                pdf = FPDF()
                pdf.add_page()
                pdf.set_font("Arial", 'B', 16)
                pdf.cell(0, 10, txt="Relatorio de Predicoes QSAR", ln=True, align='C')
                pdf.set_font("Arial", '', 12)
                pdf.cell(0, 10, txt=f"Modelo: {nome_modelo} (R2: {melhor_r2:.2f})", ln=True, align='C')
                pdf.ln(10)
                
                pdf.set_font("Arial", 'B', 10)
                pdf.cell(140, 10, "Ingrediente", 1)
                pdf.cell(40, 10, "pEC50 Previsto", 1)
                pdf.ln()
                
                pdf.set_font("Arial", '', 10)
                for _, row in previsao.iterrows():
                    # Trata caracteres não-latin1 para evitar erro no fpdf
                    ing = str(row['Ingrediente']).encode('latin-1', 'replace').decode('latin-1')
                    pec50 = f"{row['pEC50_previsto']:.3f}"
                    pdf.cell(140, 10, ing, 1)
                    pdf.cell(40, 10, pec50, 1)
                    pdf.ln()
                    
                with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
                    pdf.output(tmp.name)
                    with open(tmp.name, "rb") as f:
                        pdf_bytes = f.read()
                os.remove(tmp.name)
                
                col_pdf.download_button(
                    label="Baixar Relatório (PDF)",
                    data=pdf_bytes,
                    file_name="relatorio_previsoes.pdf",
                    mime="application/pdf",
                    use_container_width=True
                )
                
            except Exception as e:
                st.error(f"Erro na execução: {e}")
                import traceback
                st.code(traceback.format_exc())

# ==========================================
# MODO 3: VALIDAÇÃO
# ==========================================
elif modo == "3. Validação do Modelo":
    st.header("3. Validação Completa")
    st.markdown("Avalia a qualidade dos dados, domínio de aplicabilidade, y-scrambling e matriz de confusão GHS.")
    
    if st.button("Executar Etapa 3", type="primary"):
        with st.spinner("Realizando testes de validação no modelo (Pode demorar alguns minutos)..."):
            try:
                matriz = load_public_data(apenas_chlorella=usar_chlorella)
                
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

# ==========================================
# MODO 4: RELATÓRIOS
# ==========================================
elif modo == "4. Relatórios":
    st.header("Relatórios e Documentação")
    st.markdown(
        "Leitura completa dos documentos técnicos do projeto. "
        "Use as abas para navegar entre os relatórios."
    )

    aba_validacao, aba_arquitetura, aba_readme = st.tabs([
        "Validação Canônica (TCC)",
        "Arquitetura do Projeto",
        "README",
    ])

    def _render_md_com_imagens(caminho_md: str):
        """Lê um arquivo .md e renderiza seu conteúdo, exibindo imagens inline."""
        import re
        try:
            with open(caminho_md, "r", encoding="utf-8") as f:
                conteudo = f.read()
        except FileNotFoundError:
            st.warning(f"Arquivo `{caminho_md}` não encontrado.")
            return

        # Divide o conteúdo por blocos de imagem markdown: ![alt](caminho)
        partes = re.split(r'(!\[.*?\]\(.*?\))', conteudo)
        pasta_base = os.path.dirname(os.path.abspath(caminho_md))

        for parte in partes:
            match = re.match(r'!\[(.*?)\]\((.+?)\)', parte)
            if match:
                alt, src = match.group(1), match.group(2)
                # Suporte a caminhos relativos e absolutos
                if not os.path.isabs(src):
                    src = os.path.join(pasta_base, src)
                if os.path.exists(src):
                    st.image(src, caption=alt, use_container_width=True)
                else:
                    st.markdown(f"*[Imagem não encontrada: `{src}`]*")
            else:
                if parte.strip():
                    st.markdown(parte, unsafe_allow_html=False)

    with aba_validacao:
        _render_md_com_imagens("relatorio_validacao_canonica_TCC.md")

    with aba_arquitetura:
        _render_md_com_imagens("architecture.md")

    with aba_readme:
        _render_md_com_imagens("readme.md")
