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


def gerar_relatorios_e_graficos(df: pd.DataFrame, pasta_assets: str = "assets", arquivo_saida: str = "relatorio_contagem_celular.xlsx"):
    """Gera relatórios em Excel e gráficos de crescimento celular automaticamente."""
    import os
    import matplotlib.pyplot as plt

    # 1. Salva tabela limpa e ordenada em Excel
    df_ordenado = df.sort_values(["experimento", "amostra_id", "diluicao_label", "dia"])
    df_ordenado.to_excel(arquivo_saida, index=False, sheet_name="Contagens_Celulares")
    print(f"[INFO] Relatório completo em Excel salvo em: {arquivo_saida}")

    # 2. Cria pasta de assets se não existir
    os.makedirs(pasta_assets, exist_ok=True)

    # 3. Gráfico para TCC_original (6 amostras)
    df_tcc = df[df["experimento"] == "TCC_original"].copy()
    if not df_tcc.empty:
        amostras = sorted(df_tcc["amostra_id"].unique(), key=lambda x: int(x.split()[-1]) if x.split()[-1].isdigit() else 0)
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for i, amostra in enumerate(amostras):
            if i >= len(axes):
                break
            ax = axes[i]
            sub_df = df_tcc[df_tcc["amostra_id"] == amostra].sort_values("dia")

            # Controle
            ctrl_df = sub_df[["dia", "controle_normalizado"]].drop_duplicates().sort_values("dia")
            ax.plot(ctrl_df["dia"], ctrl_df["controle_normalizado"], "k--", label="Controle", marker="o", linewidth=2)

            # Diluições
            for diluicao in sub_df["diluicao_label"].dropna().unique():
                dil_df = sub_df[sub_df["diluicao_label"] == diluicao].sort_values("dia")
                ax.plot(dil_df["dia"], dil_df["contagem_normalizada"], label=diluicao, marker="s", alpha=0.8)

            ax.set_title(f"{amostra} (TCC Original)")
            ax.set_ylabel("Células/mL")
            ax.set_xlabel("Dia")
            ax.grid(True, linestyle=":", alpha=0.6)
            ax.legend(fontsize=8)

        # Ocultar subplots sobressalentes se houver
        for j in range(len(amostras), len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        caminho_tcc = os.path.join(pasta_assets, "curvas_crescimento_tcc_original.png")
        plt.savefig(caminho_tcc, dpi=300)
        plt.close()
        print(f"[INFO] Gráfico TCC Original salvo em: {caminho_tcc}")

    # 4. Gráfico para Novo_exp_diluicoes (4 amostras)
    df_novo = df[df["experimento"] == "Novo_exp_diluicoes"].copy()
    if not df_novo.empty:
        amostras = sorted(df_novo["amostra_id"].unique(), key=lambda x: int(x.split()[-1]) if x.split()[-1].isdigit() else 0)
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        for i, amostra in enumerate(amostras):
            if i >= len(axes):
                break
            ax = axes[i]
            sub_df = df_novo[df_novo["amostra_id"] == amostra].sort_values("dia")

            # Controle
            ctrl_df = sub_df[["dia", "controle_normalizado"]].drop_duplicates().sort_values("dia")
            ax.plot(ctrl_df["dia"], ctrl_df["controle_normalizado"], "k--", label="Controle", marker="o", linewidth=2)

            # Diluições
            for diluicao in sub_df["diluicao_label"].dropna().unique():
                dil_df = sub_df[sub_df["diluicao_label"] == diluicao].sort_values("dia")
                ax.plot(dil_df["dia"], dil_df["contagem_normalizada"], label=diluicao, marker="s", alpha=0.8)

            ax.set_title(f"{amostra} (Novo Exp. Diluições)")
            ax.set_ylabel("Células/mL")
            ax.set_xlabel("Dia")
            ax.grid(True, linestyle=":", alpha=0.6)
            ax.legend(fontsize=8)

        # Ocultar subplots sobressalentes se houver
        for j in range(len(amostras), len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        caminho_novo = os.path.join(pasta_assets, "curvas_crescimento_novo_exp.png")
        plt.savefig(caminho_novo, dpi=300)
        plt.close()
        print(f"[INFO] Gráfico Novo Experimento salvo em: {caminho_novo}")