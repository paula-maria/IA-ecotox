"""Combinação das matrizes de treino provenientes de múltiplas fontes
públicas (ECOTOX + EnviroTox), com controle de qualidade inter-fontes.

Uso principal:
    from fontes_externas import combinar_fontes
    matriz = combinar_fontes(matriz_ecotox, matriz_envirotox)
    modelo.treinar_modelo_publico(matriz)
"""

import numpy as np
import pandas as pd


# Limiar de desvio padrão de pEC50 entre fontes para um mesmo CAS.
# Mesmo critério usado em ecotox.montar_matriz_treino para réplicas.
LIMITE_STD_INTER_FONTE = 1.0


def _one_hot_latin_name(df: pd.DataFrame) -> pd.DataFrame:
    """Adiciona colunas binárias (0/1) para cada espécie presente na coluna
    'latin_name', com prefixo 'especie_'.

    Exemplo de colunas geradas:
        especie_Chlorella vulgaris  →  1 se a linha é dessa espécie, 0 caso contrário
        especie_Scenedesmus obliquus → 1/0

    A coluna 'latin_name' original é mantida para rastreabilidade.
    """
    if "latin_name" not in df.columns:
        return df

    df = df.copy()
    # Normaliza para evitar variações de capitalização/espaço
    df["latin_name"] = df["latin_name"].fillna("").str.strip()

    # Gera dummies; dtype int para compatibilidade com StandardScaler/RF
    dummies = pd.get_dummies(df["latin_name"], prefix="especie", dtype=int)
    # Garante que nomes de coluna sejam seguros (sem caracteres problemáticos)
    dummies.columns = [
        c.replace(" ", "_").replace(".", "_").replace("/", "_")
        for c in dummies.columns
    ]
    df = pd.concat([df, dummies], axis=1)
    return df


def combinar_fontes(
    matriz_ecotox: pd.DataFrame,
    matriz_envirotox: pd.DataFrame,
) -> pd.DataFrame:
    """Combina as matrizes de treino do ECOTOX e do EnviroTox em uma única
    matriz, com controle de qualidade inter-fontes.

    Passos:
    1. Adiciona coluna 'fonte' para rastreabilidade ('ecotox' / 'envirotox').
    2. Adiciona coluna 'latin_name' no ECOTOX (não existe lá, preenche com
       'Chlorella vulgaris' pois é o filtro de espécie usado em ecotox.py).
    3. Alinha colunas de descritores (interseção das colunas numéricas).
    4. Para CAS presentes nas DUAS fontes:
       - Calcula std de pEC50 entre as fontes.
       - std > LIMITE_STD_INTER_FONTE → reporta e descarta (não tira média cega).
       - std ≤ LIMITE_STD_INTER_FONTE → tira a média (mesma lógica de deduplicação).
    5. Para CAS exclusivos de uma fonte → inclui diretamente.
    6. Gera colunas one-hot para 'latin_name'.
    7. Retorna a matriz final pronta para modelo.treinar_modelo_publico.

    Parâmetros
    ----------
    matriz_ecotox : DataFrame com colunas de descritores + 'pEC50' + 'cas'
    matriz_envirotox : DataFrame com colunas de descritores + 'pEC50' + 'cas'
                       + 'latin_name'

    Retorna
    -------
    DataFrame combinado com descritores + 'pEC50' + 'cas' + 'latin_name'
    + colunas 'especie_*' (one-hot) + coluna 'fonte' (rastreabilidade).
    """
    # ------------------------------------------------------------------ #
    # 1. Adiciona marcadores de fonte e latin_name ao ECOTOX              #
    # ------------------------------------------------------------------ #
    ecotox = matriz_ecotox.copy()
    envirotox = matriz_envirotox.copy()

    ecotox["fonte"] = "ecotox"
    envirotox["fonte"] = "envirotox"

    # ECOTOX não tem latin_name — preenche com a espécie filtrada (Chlorella vulgaris)
    if "latin_name" not in ecotox.columns:
        ecotox["latin_name"] = "Chlorella vulgaris"

    # ------------------------------------------------------------------ #
    # 2. Alinha colunas de descritores (interseção das numéricas)         #
    # ------------------------------------------------------------------ #
    colunas_nao_desc = {"pEC50", "cas", "latin_name", "fonte"}

    cols_ecotox_desc = set(ecotox.select_dtypes(include="number").columns) - colunas_nao_desc
    cols_envirotox_desc = set(envirotox.select_dtypes(include="number").columns) - colunas_nao_desc
    descritores_comuns = sorted(cols_ecotox_desc & cols_envirotox_desc)

    colunas_desc_faltantes_ecotox = cols_envirotox_desc - cols_ecotox_desc
    colunas_desc_faltantes_envirotox = cols_ecotox_desc - cols_envirotox_desc
    if colunas_desc_faltantes_ecotox:
        print(
            f"[AVISO] {len(colunas_desc_faltantes_ecotox)} descritores do EnviroTox "
            f"ausentes no ECOTOX — ignorados na combinação."
        )
    if colunas_desc_faltantes_envirotox:
        print(
            f"[AVISO] {len(colunas_desc_faltantes_envirotox)} descritores do ECOTOX "
            f"ausentes no EnviroTox — ignorados na combinação."
        )

    colunas_finais = descritores_comuns + ["pEC50", "cas", "latin_name", "fonte"]
    ecotox = ecotox[[c for c in colunas_finais if c in ecotox.columns]]
    envirotox = envirotox[[c for c in colunas_finais if c in envirotox.columns]]

    # ------------------------------------------------------------------ #
    # 3. Identifica CAS compartilhados vs. exclusivos                     #
    # ------------------------------------------------------------------ #
    cas_ecotox = set(ecotox["cas"].dropna().astype(str))
    cas_envirotox = set(envirotox["cas"].dropna().astype(str))

    cas_compartilhados = cas_ecotox & cas_envirotox
    cas_exclusivos_ecotox = cas_ecotox - cas_envirotox
    cas_exclusivos_envirotox = cas_envirotox - cas_ecotox

    print(f"\n[INFO] CAS no ECOTOX:     {len(cas_ecotox)}")
    print(f"[INFO] CAS no EnviroTox:  {len(cas_envirotox)}")
    print(f"[INFO] CAS compartilhados (ambas as fontes): {len(cas_compartilhados)}")
    print(f"[INFO] CAS exclusivos do ECOTOX:             {len(cas_exclusivos_ecotox)}")
    print(f"[INFO] CAS exclusivos do EnviroTox (novos):  {len(cas_exclusivos_envirotox)}")

    # ------------------------------------------------------------------ #
    # 4. Trata CAS compartilhados com controle de variância inter-fonte   #
    # ------------------------------------------------------------------ #
    # Empilha as linhas dos CAS compartilhados de ambas as fontes
    df_todas = pd.concat([ecotox, envirotox], ignore_index=True)
    df_compartilhados = df_todas[
        df_todas["cas"].astype(str).isin(cas_compartilhados)
    ].copy()

    # Calcula desvio padrão do pEC50 por CAS entre as fontes
    std_inter_fonte = df_compartilhados.groupby("cas")["pEC50"].std()
    cas_alta_divergencia = std_inter_fonte[std_inter_fonte > LIMITE_STD_INTER_FONTE].index
    cas_concordantes = std_inter_fonte[
        std_inter_fonte.isna() | (std_inter_fonte <= LIMITE_STD_INTER_FONTE)
    ].index

    if len(cas_alta_divergencia) > 0:
        print(
            f"\n[AVISO] {len(cas_alta_divergencia)} CAS descartados por alta divergência "
            f"inter-fonte (std pEC50 > {LIMITE_STD_INTER_FONTE}):"
        )
        divergencia_detalhada = std_inter_fonte[cas_alta_divergencia].sort_values(ascending=False)
        for cas, std_val in divergencia_detalhada.items():
            vals = df_compartilhados[df_compartilhados["cas"] == cas]["pEC50"].values
            print(f"  CAS {cas}: std={std_val:.2f} | valores={np.round(vals, 2)}")

    # Média dos CAS concordantes
    df_concordantes = df_compartilhados[
        df_compartilhados["cas"].isin(cas_concordantes)
    ].copy()

    colunas_num = [c for c in descritores_comuns + ["pEC50"] if c in df_concordantes.columns]
    df_concordantes_media = (
        df_concordantes.groupby("cas", sort=False)[colunas_num]
        .mean()
        .reset_index()
    )
    # latin_name: espécie mais específica — dá preferência a Chlorella vulgaris
    latin_concordantes = (
        df_concordantes.groupby("cas", sort=False)["latin_name"]
        .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else "")
        .reset_index()
    )
    df_concordantes_final = df_concordantes_media.merge(latin_concordantes, on="cas")
    df_concordantes_final["fonte"] = "ecotox+envirotox"

    # ------------------------------------------------------------------ #
    # 5. CAS exclusivos — inclui diretamente                              #
    # ------------------------------------------------------------------ #
    df_excl_ecotox = ecotox[ecotox["cas"].astype(str).isin(cas_exclusivos_ecotox)].copy()
    df_excl_envirotox = envirotox[
        envirotox["cas"].astype(str).isin(cas_exclusivos_envirotox)
    ].copy()

    # ------------------------------------------------------------------ #
    # 6. Concatena tudo                                                   #
    # ------------------------------------------------------------------ #
    df_final = pd.concat(
        [df_excl_ecotox, df_excl_envirotox, df_concordantes_final],
        ignore_index=True,
    )

    # ------------------------------------------------------------------ #
    # 7. One-hot encoding de latin_name                                   #
    # ------------------------------------------------------------------ #
    df_final = _one_hot_latin_name(df_final)

    # ------------------------------------------------------------------ #
    # 8. Resumo final                                                     #
    # ------------------------------------------------------------------ #
    n_final = len(df_final)
    print(f"\n[INFO] ===== MATRIZ COMBINADA =====")
    print(f"[INFO] Linhas finais:         {n_final} compostos únicos (CAS)")
    print(f"[INFO] CAS exclusivos ECOTOX: {len(cas_exclusivos_ecotox)}")
    print(f"[INFO] CAS exclusivos EnviroTox (novos): {len(cas_exclusivos_envirotox)}")
    print(f"[INFO] CAS de ambas as fontes (média usada): {len(cas_concordantes)}")
    print(
        f"[INFO] CAS descartados por divergência inter-fonte: {len(cas_alta_divergencia)}"
    )

    # Distribuição GHS
    try:
        from validacao import classificar_toxicidade_ghs
        classes = [
            classificar_toxicidade_ghs(row["pEC50"], row["MolWt"])
            for _, row in df_final.iterrows()
            if pd.notna(row["pEC50"]) and pd.notna(row.get("MolWt", np.nan))
        ]
        print("\n[INFO] Distribuição GHS na matriz combinada:")
        print(pd.Series(classes).value_counts().to_string())
    except Exception as e:
        print(f"[AVISO] Não foi possível calcular distribuição GHS: {e}")

    return df_final
