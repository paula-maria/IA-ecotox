"""Ponto de entrada do pipeline. Uso:

    python3 main.py dados       # só normaliza a planilha própria e mostra os descritores
    python3 main.py publico     # Fase 1 (treino no ECOTOX) + Fase 2 (validação externa)

Ver README.md para instalação e pré-requisitos (em especial, o download
manual do ECOTOX antes de rodar "publico").
"""

import sys

import pandas as pd

import dados
import ecotox
import modelo
import validacao


def rodar_dados(caminho_planilha: str = dados.ARQUIVO_PADRAO):
    df = dados.carregar_dados_inibicao(caminho_planilha)
    print(df.head())

    # Gera automaticamente os relatórios em Excel e gráficos das curvas de crescimento
    dados.gerar_relatorios_e_graficos(df)

    descritores_df = dados.carregar_descritores_ingredientes(caminho_planilha)
    print(descritores_df)


def rodar_publico(caminho_planilha: str = dados.ARQUIVO_PADRAO):
    dados_ecotox = ecotox.carregar_ecotox_algas()
    dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
    matriz = ecotox.montar_matriz_treino(dados_ecotox)

    modelo_treinado, colunas_x, scaler, melhor_r2, nome_modelo = modelo.treinar_modelo_publico(matriz)

    ranking = modelo.importancia_descritores(modelo_treinado, colunas_x, scaler, matriz)
    print("\nDescritores que mais influenciam a toxicidade prevista (Permutation Importance):")
    print(ranking)

    previsao = modelo.validar_externamente(modelo_treinado, colunas_x, caminho_planilha, scaler)
    print("\nPrevisão para os ativos amazônicos:")
    print(previsao)


def rodar_validacao(caminho_planilha: str = dados.ARQUIVO_PADRAO):
    dados_ecotox = ecotox.carregar_ecotox_algas()
    dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
    matriz = ecotox.montar_matriz_treino(dados_ecotox)

    relatorio_qualidade = validacao.checar_qualidade_dados(matriz)
    validacao.imprimir_relatorio_qualidade(relatorio_qualidade)

    colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
    X, y = matriz[colunas_x], matriz["pEC50"]

    print("\n=== Y-scrambling (o modelo aprendeu sinal real?) ===")
    resultado_scrambling = validacao.teste_y_scrambling(X, y)
    print(resultado_scrambling)
    if not resultado_scrambling["passou_no_teste"]:
        print("[ATENÇÃO] R² real não se distanciou o suficiente do embaralhado "
              "— revisar o modelo antes de reportar como preditivo.")


    print("\n=== Domínio de aplicabilidade (ativos amazônicos) ===")
    ingredientes = dados.carregar_descritores_ingredientes(caminho_planilha)
    colunas_comuns = [c for c in colunas_x if c in ingredientes.columns]
    leverage_df = validacao.calcular_leverage(X[colunas_comuns], ingredientes[colunas_comuns])
    leverage_df["Ingrediente"] = ingredientes["Ingrediente"].values
    print(f"h* (limite de corte) = {leverage_df.attrs['h_estrela']:.3f}")
    print(leverage_df[["Ingrediente", "leverage", "dentro_do_dominio"]])

    print("\n=== Distribuição de Classes GHS (Treino) ===")
    from validacao import classificar_toxicidade_ghs

    classes_obs = [classificar_toxicidade_ghs(row["pEC50"], row["MolWt"]) for _, row in matriz.iterrows()]
    contagem = pd.Series(classes_obs).value_counts()
    print(contagem)

    print("\n=== Matriz de confusão e Análise de Resíduos (Regressão pEC50 contínuo) ===")
    # Usa validação cruzada no pipeline para extrair y_pred do melhor modelo sem vazar dados
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import KFold, cross_val_predict
    
    # Para simplificar a validação e gerar gráficos rapidamente, vamos usar 
    # o modelo treinado de forma mais direta, mas idealmente re-rodariamos o pipeline.
    # Usaremos o RandomForest otimizado como representativo para os gráficos, 
    # mas o grid_search no modelo publico é onde a mágica acontece.
    
    # Vamos gerar os gráficos:
    # 1. Matriz de confusão clássica
    modelo_cv, col_x, scaler_cv = modelo.treinar_modelo_publico(matriz)
    
    # Para o CV limpo, refazemos:
    X_scaled = scaler_cv.transform(X)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    pred_cv = cross_val_predict(modelo_cv, X_scaled, y, cv=kf)

    df_matriz_reg = validacao.matriz_confusao_toxicidade(y, pd.Series(pred_cv), matriz["MolWt"])
    print(df_matriz_reg)
    caminho_png_reg = validacao.plotar_matriz_confusao(df_matriz_reg, "matriz_confusao_regressao.png")
    
    caminho_png_residuo = validacao.plotar_analise_residuos(y, pd.Series(pred_cv), "analise_residuos.png")
    print(f"[INFO] Gráfico de Análise de Resíduos salvo em: {caminho_png_residuo}")

    print("\n=== Matriz de confusão (Classificador com class_weight='balanced') ===")
    modelo_clf, pred_cv_clf, y_class = modelo.treinar_classificador_ghs(matriz)
    from sklearn.metrics import confusion_matrix
    from validacao import LABELS_GHS
    matriz_clf = confusion_matrix(y_class, pred_cv_clf, labels=LABELS_GHS)
    df_matriz_clf = pd.DataFrame(matriz_clf, index=LABELS_GHS, columns=LABELS_GHS)
    df_matriz_clf.index.name = "Observado"
    df_matriz_clf.columns.name = "Previsto"
    print(df_matriz_clf)
    caminho_png_clf = validacao.plotar_matriz_confusao(df_matriz_clf, "matriz_confusao_classificador.png")
    print(f"\nMatrizes de confusão salvas em: {caminho_png_reg} e {caminho_png_clf}")


if __name__ == "__main__":
    comando = sys.argv[1] if len(sys.argv) > 1 else "dados"

    if comando == "dados":
        rodar_dados()
    elif comando == "publico":
        rodar_publico()
    elif comando == "validar":
        rodar_validacao()
    else:
        print(f"Comando desconhecido: {comando!r}. Use 'dados', 'publico' ou 'validar'.")