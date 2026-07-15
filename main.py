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

    modelo_treinado, colunas_x = modelo.treinar_modelo_publico(matriz)

    ranking = modelo.importancia_descritores(modelo_treinado, colunas_x)
    print("\nDescritores que mais influenciam a toxicidade prevista:")
    print(ranking)

    previsao = modelo.validar_externamente(modelo_treinado, colunas_x, caminho_planilha)
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

    print("\n=== Matriz de confusão (categorias GHS) ===")
    # usa o mesmo esquema de validação cruzada 5-fold do treino público como
    # 'previsto', comparado ao 'observado' real da base — dá uma view geral
    # de acerto por categoria antes mesmo de aplicar aos dados amazônicos
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold, cross_val_predict
    # n_jobs=-1: usa todas as CPUs disponíveis para paralelizar a validação cruzada
    modelo_cv = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    pred_cv = cross_val_predict(modelo_cv, X, y, cv=kf)

    df_matriz = validacao.matriz_confusao_toxicidade(y, pd.Series(pred_cv), matriz["MolWt"])
    print(df_matriz)
    caminho_png = validacao.plotar_matriz_confusao(df_matriz)
    print(f"Matriz de confusão salva em: {caminho_png}")


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