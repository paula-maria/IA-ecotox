"""Ponto de entrada do pipeline. Uso:

    python3 main.py dados       # só normaliza a planilha própria e mostra os descritores
    python3 main.py publico     # Fase 1 (treino no ECOTOX) + Fase 2 (validação externa)

Ver README.md para instalação e pré-requisitos (em especial, o download
manual do ECOTOX antes de rodar "publico").
"""

import sys

import dados
import ecotox
import modelo


def rodar_dados(caminho_planilha: str = dados.ARQUIVO_PADRAO):
    df = dados.carregar_dados_inibicao(caminho_planilha)
    print(df.head())

    descritores_df = dados.carregar_descritores_ingredientes(caminho_planilha)
    print(descritores_df)


def rodar_publico(caminho_planilha: str = dados.ARQUIVO_PADRAO):
    dados_ecotox = ecotox.carregar_ecotox_algas()
    dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
    matriz = ecotox.montar_matriz_treino(dados_ecotox)

    modelo_treinado, colunas_x = modelo.treinar_modelo_publico(matriz)

    previsao = modelo.validar_externamente(modelo_treinado, colunas_x, caminho_planilha)
    print(previsao)


if __name__ == "__main__":
    comando = sys.argv[1] if len(sys.argv) > 1 else "dados"

    if comando == "dados":
        rodar_dados()
    elif comando == "publico":
        rodar_publico()
    else:
        print(f"Comando desconhecido: {comando!r}. Use 'dados' ou 'publico'.")