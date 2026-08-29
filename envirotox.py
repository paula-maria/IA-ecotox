"""Carregamento e pré-processamento da base pública EnviroTox
(envirotoxdatabase.org), espelhando a estrutura de ecotox.py.

Diferenças-chave em relação ao ECOTOX:
  - Os arquivos fonte são Excel (não CSV com pipes) com 3 sheets: test,
    substance, taxonomy.
  - O SMILES já vem pronto na coluna "Desalted Canonical SMILES" da sheet
    substance — PubChem só é consultado como fallback para CAS sem SMILES.
  - A unidade de concentração no EnviroTox (nos conjuntos baixados) é sempre
    mg/L — a lógica de conversão reaproveita _conc_para_mol_por_l de ecotox.py.
  - A coluna latin_name é mantida na matriz final para uso como descritor
    categórico (one-hot encoding) em fontes_externas.py.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from descritores import calcular_descritores

# Reusa as constantes e a lógica de conversão de unidades de ecotox.py,
# evitando duplicação.
from ecotox import (
    ENDPOINTS_CRESCIMENTO,
    _conc_para_mol_por_l,
    cas_para_smiles,
)

# Caminhos padrão — espelham a estrutura de pastas do projeto.
DIR_ALGAE = Path(".")
DIR_CHLORELLA = Path(".")

# Arquivo canônico de cada conjunto (os 3 arquivos de cada pasta são idênticos).
ALGAE_TEST_FILE = DIR_ALGAE / "algae.xlsx"
ALGAE_SUBSTANCE_FILE = DIR_ALGAE / "algae.xlsx"
CHLORELLA_TEST_FILE = DIR_CHLORELLA / "chlorella.xlsx"
CHLORELLA_SUBSTANCE_FILE = DIR_CHLORELLA / "chlorella.xlsx"

# Sheet names dentro de cada Excel
SHEET_TEST = "test"
SHEET_SUBSTANCE = "substance"


def _normalizar_cas(serie: pd.Series) -> pd.Series:
    """Normaliza a coluna CAS para string sem espaços, removendo sufixos do tipo
    'Metalgrp.Se' (CAS especiais de metais agrupados no EnviroTox).
    Retorna strings limpas para permitir merge consistente entre test e substance."""
    return serie.astype(str).str.strip()


def _carregar_sheet_test(caminho: Path) -> pd.DataFrame:
    """Lê a sheet 'test' de um Excel EnviroTox e retorna apenas as colunas
    relevantes, com o filtro de endpoints aplicado."""
    colunas = ["CAS", "Latin name", "Test statistic", "Effect value", "Unit"]
    df = pd.read_excel(caminho, sheet_name=SHEET_TEST, usecols=colunas)
    df["CAS"] = _normalizar_cas(df["CAS"])
    df["Effect value"] = pd.to_numeric(df["Effect value"], errors="coerce")
    df = df.dropna(subset=["CAS", "Effect value"])
    df = df[df["Effect value"] > 0]
    df = df[df["Test statistic"].isin(ENDPOINTS_CRESCIMENTO)]
    return df


def _carregar_sheet_substance(caminho: Path) -> pd.DataFrame:
    """Lê a sheet 'substance' de um Excel EnviroTox e retorna apenas as colunas
    relevantes para o cálculo de pEC50 e descritores."""
    colunas = ["CAS", "Desalted Canonical SMILES", "MW (g/mol)"]
    df = pd.read_excel(caminho, sheet_name=SHEET_SUBSTANCE, usecols=colunas)
    df["CAS"] = _normalizar_cas(df["CAS"])
    df = df.rename(columns={
        "Desalted Canonical SMILES": "smiles_envirotox",
        "MW (g/mol)": "mol_wt_envirotox",
    })
    # Deduplicar por CAS — substance tem uma linha por composto
    df = df.drop_duplicates(subset=["CAS"])
    return df


def _resolver_smiles_faltantes(df: pd.DataFrame) -> pd.DataFrame:
    """Para linhas sem SMILES na planilha EnviroTox, tenta resolver via PubChem
    (reusa o cache pubchem_cache.json de ecotox.py)."""
    import json

    cache_file = Path("pubchem_cache.json")
    cache: dict = {}
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                cache = json.load(f)
        except Exception as e:
            print(f"[AVISO] Erro ao ler cache do PubChem: {e}")

    if cache and all(v is None for v in cache.values()):
        print("[AVISO] Cache do PubChem está corrompido (todos None) — descartando.")
        cache = {}

    sem_smiles = df[df["smiles_envirotox"].isna()]["CAS"].unique()
    if len(sem_smiles) == 0:
        return df

    new_cas = [c for c in sem_smiles if c not in cache]
    if new_cas:
        print(
            f"[INFO] EnviroTox: {len(sem_smiles)} CAS sem SMILES — "
            f"consultando PubChem para {len(new_cas)} novos..."
        )
        total = len(new_cas)
        for i, c in enumerate(new_cas, 1):
            cas_para_smiles(c, cache)
            if i % 50 == 0 or i == total:
                print(f"  -> Progresso PubChem: {i}/{total}")
                try:
                    with open(cache_file, "w") as f:
                        json.dump(cache, f, indent=4)
                except Exception:
                    pass
    else:
        print(f"[INFO] EnviroTox: {len(sem_smiles)} CAS sem SMILES já estão no cache.")

    # Preenche smiles_envirotox apenas onde estava nulo
    mascara_nulo = df["smiles_envirotox"].isna()
    df = df.copy()
    df.loc[mascara_nulo, "smiles_envirotox"] = df.loc[mascara_nulo, "CAS"].map(cache)
    return df


def carregar_envirotox(
    caminho_test_algae: Path = ALGAE_TEST_FILE,
    caminho_substance_algae: Path = ALGAE_SUBSTANCE_FILE,
    caminho_test_chlorella: Path = CHLORELLA_TEST_FILE,
    caminho_substance_chlorella: Path = CHLORELLA_SUBSTANCE_FILE,
) -> pd.DataFrame:
    """Carrega e combina os datasets EnviroTox de algas e Chlorella.

    Estratégia:
      - Lê os dois conjuntos (algae + chlorella) e os combina, desduplicando
        por CAS para que a espécie mais específica (*Chlorella vulgaris*)
        seja preservada quando o CAS aparecer nos dois.
      - Filtra Test statistic em ENDPOINTS_CRESCIMENTO (EC50, IC50, NOEC, LOEC).
      - Faz merge com a sheet substance para obter SMILES e MW.
      - Para CAS sem SMILES na planilha, consulta PubChem como fallback.
      - Mantém a coluna 'latin_name' para uso posterior como feature categórica.

    Retorna DataFrame com colunas:
        CAS, latin_name, Test statistic, Effect value, Unit,
        smiles_envirotox, mol_wt_envirotox
    """
    print("[INFO] Carregando EnviroTox — conjunto amplo (algae)...")
    test_algae = _carregar_sheet_test(caminho_test_algae)
    substance_algae = _carregar_sheet_substance(caminho_substance_algae)
    print(f"       {len(test_algae)} linhas após filtro de endpoints.")

    print("[INFO] Carregando EnviroTox — conjunto específico (Chlorella)...")
    test_chlorella = _carregar_sheet_test(caminho_test_chlorella)
    substance_chlorella = _carregar_sheet_substance(caminho_substance_chlorella)
    print(f"       {len(test_chlorella)} linhas após filtro de endpoints.")

    # Combina as sheets substance dos dois conjuntos — preferindo chlorella
    # quando o mesmo CAS aparece nos dois (mais específico).
    substance_combinada = pd.concat(
        [substance_algae, substance_chlorella], ignore_index=True
    )
    # Mantém a primeira ocorrência por CAS; como chlorella vem por último,
    # invertemos e fazemos drop_duplicates(keep='last') para dar preferência.
    substance_combinada = substance_combinada.drop_duplicates(subset=["CAS"], keep="last")

    # Combina os datasets de teste, depois desduplicamos por CAS preservando
    # a espécie mais específica (chlorella tem precedência).
    test_combinado = pd.concat(
        [test_algae, test_chlorella], ignore_index=True
    )

    print(
        f"[INFO] Total bruto combinado: {len(test_combinado)} linhas | "
        f"{test_combinado['CAS'].nunique()} CAS únicos."
    )

    # Merge com substance para obter SMILES e MW
    df = test_combinado.merge(substance_combinada, on="CAS", how="left")
    df = df.rename(columns={"Latin name": "latin_name"})

    # Resolve SMILES ausentes via PubChem
    df = _resolver_smiles_faltantes(df)

    # Descarta linhas sem SMILES mesmo após consulta ao PubChem
    n_antes = len(df)
    df = df.dropna(subset=["smiles_envirotox"])
    n_descartados = n_antes - len(df)
    if n_descartados > 0:
        print(f"[AVISO] {n_descartados} linhas descartadas por CAS sem SMILES resolvível.")

    df = df.drop_duplicates()
    print(
        f"[INFO] EnviroTox carregado: {len(df)} linhas | "
        f"{df['CAS'].nunique()} CAS | "
        f"{df['latin_name'].nunique()} espécies."
    )
    return df


def montar_matriz_envirotox(df: pd.DataFrame) -> pd.DataFrame:
    """Constrói a matriz de treino a partir do DataFrame bruto do EnviroTox.

    Espelha ecotox.montar_matriz_treino, com as seguintes diferenças:
      - Usa 'smiles_envirotox' em vez de 'smiles'.
      - Usa 'Effect value' em vez de 'conc1_mean' e 'Unit' em vez de 'conc1_unit'.
      - Mantém 'latin_name' na matriz final (como coluna string; codificação
        one-hot é feita em fontes_externas.combinar_fontes).
      - Usa mol_wt_envirotox como peso molar auxiliar quando disponível
        (RDKit pode diferir levemente; usa-se o RDKit como fonte primária).

    Deduplicação por CAS: mesmo critério de ecotox.py (std pEC50 ≤ 1.0).
    Para latin_name: preserva a espécie mais frequente por CAS após o agrupamento.
    """
    registros = []
    skipped_unit = 0
    skipped_invalid = 0

    for _, row in df.iterrows():
        smiles = row["smiles_envirotox"]
        desc = calcular_descritores(str(smiles).strip())
        if desc is None:
            skipped_invalid += 1
            continue

        unit_raw = row.get("Unit", "mg/L")  # fallback — EnviroTox quase sempre é mg/L
        conc_molar = _conc_para_mol_por_l(
            float(row["Effect value"]), unit_raw, desc["MolWt"]
        )
        if conc_molar is None:
            skipped_unit += 1
            continue
        if conc_molar <= 0:
            skipped_invalid += 1
            continue

        registros.append({
            **desc,
            "pEC50": -np.log10(conc_molar),
            "cas": row["CAS"],
            "latin_name": str(row.get("latin_name", "")).strip(),
        })

    if skipped_unit > 0:
        print(
            f"[AVISO] EnviroTox: {skipped_unit} linhas descartadas por unidade de "
            f"concentração desconhecida (adicionar em _UNIT_TO_G_PER_L/_UNIT_TO_MOL_PER_L se necessário)."
        )
    if skipped_invalid > 0:
        print(
            f"[INFO] EnviroTox: {skipped_invalid} linhas descartadas por SMILES inválido "
            f"ou concentração ≤ 0."
        )

    matriz = pd.DataFrame(registros)
    if matriz.empty or "pEC50" not in matriz.columns:
        raise ValueError(
            "A matriz EnviroTox ficou vazia. Verifique:\n"
            "  1. Se os arquivos Excel estão nos caminhos esperados (algae/ e chlorella/).\n"
            "  2. Se os endpoints filtrados estão presentes (EC50, IC50, NOEC, LOEC).\n"
            "  3. Se os SMILES foram resolvidos corretamente."
        )

    n_antes = len(matriz)
    n_cas_antes = matriz["cas"].nunique()

    # --- Deduplicação: uma linha por composto (CAS) ----------------------------
    limite_std_pec50 = 1.0
    std_pec50 = matriz.groupby("cas", sort=False)["pEC50"].std()
    cas_validos = std_pec50[std_pec50.isna() | (std_pec50 <= limite_std_pec50)].index
    n_cas_descartados = len(std_pec50) - len(cas_validos)

    matriz_filtrada = matriz[matriz["cas"].isin(cas_validos)].copy()
    linhas_descartadas = n_antes - len(matriz_filtrada)

    # Para colunas numéricas: média por CAS
    colunas_num = matriz_filtrada.select_dtypes(include="number").columns.tolist()
    matriz_num = (
        matriz_filtrada.groupby("cas", sort=False)[colunas_num]
        .mean()
        .reset_index()
    )

    # Para latin_name: espécie mais frequente por CAS (moda)
    latin_por_cas = (
        matriz_filtrada.groupby("cas", sort=False)["latin_name"]
        .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else "")
        .reset_index()
    )

    matriz = matriz_num.merge(latin_por_cas, on="cas", how="left")

    n_depois = len(matriz)
    linhas_condensadas = len(matriz_filtrada) - n_depois

    print(f"[INFO] Matriz EnviroTox: {n_antes} amostras brutas originais.")
    if n_cas_descartados > 0:
        print(
            f"[INFO] {n_cas_descartados} compostos (CAS) descartados por alta divergência "
            f"(std pEC50 > {limite_std_pec50}). ({linhas_descartadas} linhas afetadas)."
        )
    print(
        f"[INFO] Resultado final: {n_depois} compostos únicos (CAS) após deduplicação "
        f"({linhas_condensadas} réplicas condensadas pela média)."
    )
    print(f"[INFO] pEC50 NaN após deduplicação: {matriz['pEC50'].isna().sum()}.")

    # Espécies representadas
    contagem_especies = matriz["latin_name"].value_counts()
    print("[INFO] Espécies na matriz EnviroTox (top 10):")
    print(contagem_especies.head(10).to_string())

    return matriz
