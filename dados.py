"""Carregamento e normalização dos dados próprios (planilha do TCC/PJC) e
montagem do descritor de mistura (Fase 3, pendente).
"""

import numpy as np
import pandas as pd

from descritores import calcular_descritores

ARQUIVO_PADRAO = "Dados_QSAR_Saile_PJC2026.xlsx"


def carregar_dados_inibicao(caminho: str = ARQUIVO_PADRAO) -> pd.DataFrame:
    """Lê a aba '4_Formato_Python_QSAR' e normaliza as unidades: o TCC
    original está em células/mL, o novo experimento em ×10⁵ células/mL."""
    df = pd.read_excel(caminho, sheet_name="4_Formato_Python_QSAR", header=1)
    df = df.dropna(subset=["amostra_id"]).copy()

    fator = df["unidade"].map({"celulas/mL": 1, "x10^5 celulas/mL": 1e5})
    df["contagem_normalizada"] = df["contagem_media"] * fator
    df["controle_normalizado"] = df["controle_media"] * fator

    return df


def carregar_descritores_ingredientes(caminho: str = ARQUIVO_PADRAO) -> pd.DataFrame:
    """Lê a aba '5_SMILES_Ingredientes' e calcula os descritores de cada um.
    Ingredientes com SMILES inválido são reportados e descartados."""
    df = pd.read_excel(caminho, sheet_name="5_SMILES_Ingredientes", header=1)
    df = df.dropna(subset=["Ingrediente"]).copy()
    df["SMILES_canonical"] = df["SMILES_canonical"].str.strip()

    registros = []
    for _, row in df.iterrows():
        desc = calcular_descritores(row["SMILES_canonical"])
        if desc is None:
            print(f"[ATENÇÃO] SMILES inválido/suspeito para '{row['Ingrediente']}': "
                  f"{row['SMILES_canonical']!r} — revisar no PubChem.")
            continue
        registros.append({"Ingrediente": row["Ingrediente"], "CAS": row["CAS"], **desc})

    return pd.DataFrame(registros)


# ---------------------------------------------------------------------------
# FASE 3 (PENDENTE) — descritor de mistura, depende da tabela de composição
# ---------------------------------------------------------------------------
# Formato esperado, a ser preenchido pela sua irmã a partir das fichas
# técnicas das formulações:
#
#   amostra_id   | Ingrediente              | percentual
#   -------------|--------------------------|-----------
#   Amostra 9    | Terpinen-4-ol            | 2.5
#   Amostra 9    | Glicerina (Glycerol)     | 5.0
#
# composicao = pd.read_excel(ARQUIVO_PADRAO, sheet_name="6_Composicao_Formulas")

def descritor_da_mistura(tabela_composicao: pd.DataFrame,
                          descritores_ingredientes: pd.DataFrame) -> pd.DataFrame:
    """Calcula o descritor de uma FORMULAÇÃO como média ponderada pelo % de
    cada ingrediente. Só funciona quando a tabela de composição acima
    existir — ver 'Fora do escopo atual' em architecture.md."""
    merged = tabela_composicao.merge(descritores_ingredientes, on="Ingrediente")
    cols_desc = [c for c in descritores_ingredientes.columns if c not in ("Ingrediente", "CAS")]

    def media_ponderada(g):
        pesos = g["percentual"] / g["percentual"].sum()
        return pd.Series({c: np.average(g[c], weights=pesos) for c in cols_desc})

    return merged.groupby("amostra_id").apply(media_ponderada).reset_index()