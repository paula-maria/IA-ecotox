"""Carregamento e filtro da base pública ECOTOX (US EPA) e conversão
CAS → SMILES via PubChem. Ver README.md para instruções de download do
ECOTOX (não é feito automaticamente por este módulo).
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from descritores import calcular_descritores

ECOTOX_DIR = Path("./ecotox_ascii")

# Espécie exata usada no ensaio de laboratório (OECD TG 201) — sensibilidade
# a um mesmo composto varia entre espécies, então o ideal é treinar só com
# a espécie que ela realmente usa, não o gênero inteiro.
ESPECIE_ALVO = "Chlorella vulgaris"

# Fallback: se sobrarem poucos registros (ver LIMIAR_MINIMO_REGISTROS abaixo),
# amplia para o gênero — documentar essa decisão explicitamente no texto do
# TCC se isso acontecer, não deixar implícito.
GENEROS_ALGA = ["Chlorella", "Raphidocelis", "Pseudokirchneriella",
                "Scenedesmus", "Selenastrum", "Desmodesmus"]
LIMIAR_MINIMO_REGISTROS = 20

ENDPOINTS_CRESCIMENTO = ["EC50", "IC50", "NOEC", "LOEC"]


def carregar_ecotox_algas(diretorio: Path = ECOTOX_DIR,
                           especie: str = ESPECIE_ALVO) -> pd.DataFrame:
    tests = pd.read_csv(diretorio / "tests.txt", sep="|", low_memory=False,
                         encoding="latin1", on_bad_lines="skip")
    results = pd.read_csv(diretorio / "results.txt", sep="|", low_memory=False,
                           encoding="latin1", on_bad_lines="skip")
    species = pd.read_csv(diretorio / "species.txt", sep="|", low_memory=False,
                           encoding="latin1", on_bad_lines="skip")

    for df in (tests, results, species):
        df.columns = [c.strip().lower() for c in df.columns]

    species_alga = species[species["latin_name"].str.strip().str.lower()
                            == especie.lower()]
    origem_filtro = f"espécie exata ({especie})"

    if len(species_alga) == 0 or len(tests.merge(
            species_alga, on="species_number", how="inner")) < LIMIAR_MINIMO_REGISTROS:
        print(f"[AVISO] Poucos ou nenhum registro para '{especie}' isoladamente — "
              f"ampliando para o gênero. Documentar essa decisão no TCC.")
        padrao_genero = "|".join(GENEROS_ALGA)
        species_alga = species[species["latin_name"].str.contains(
            padrao_genero, case=False, na=False)]
        origem_filtro = "gênero (fallback)"

    print(f"[INFO] Filtro de espécie usado: {origem_filtro} — "
          f"{len(species_alga)} espécie(s) encontrada(s).")

    df = tests.merge(species_alga, on="species_number", how="inner")
    df = df.merge(results, on="test_id", how="inner")
    df = df[df["endpoint"].isin(ENDPOINTS_CRESCIMENTO)]

    colunas = ["test_id", "result_id", "test_cas", "latin_name",
               "endpoint", "effect", "conc1_mean", "conc1_unit"]
    df = df[[c for c in colunas if c in df.columns]].dropna(
        subset=["test_cas", "conc1_mean"])

    df["test_cas"] = df["test_cas"].astype(str).str.replace(r"\D", "", regex=True)

    return df.drop_duplicates()


def cas_para_smiles(cas: str, cache: dict) -> str | None:
    if cas in cache:
        return cache[cas]

    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/"
           f"RegistryID/{cas}/property/CanonicalSMILES/JSON")
    try:
        resp = requests.get(url, timeout=10)
        time.sleep(0.25)  # PubChem pede no máx. ~5 requisições/segundo
        if resp.status_code == 200:
            smiles = resp.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            cache[cas] = smiles
            return smiles
    except Exception:
        pass

    cache[cas] = None
    return None


def anexar_smiles(df: pd.DataFrame) -> pd.DataFrame:
    cache: dict[str, str | None] = {}
    df = df.copy()
    df["smiles"] = df["test_cas"].apply(lambda c: cas_para_smiles(c, cache))
    return df.dropna(subset=["smiles"])


def montar_matriz_treino(df: pd.DataFrame) -> pd.DataFrame:
    """Junta descritores + variável-resposta em escala molar (pEC50)."""
    registros = []
    for _, row in df.iterrows():
        desc = calcular_descritores(row["smiles"])
        if desc is None:
            continue
        # conc1_mean normalmente vem em mg/L; convertendo para mol/L com o
        # peso molecular do próprio composto.
        conc_molar = (row["conc1_mean"] / 1000) / desc["MolWt"]
        if conc_molar <= 0:
            continue
        registros.append({**desc, "pEC50": -np.log10(conc_molar), "cas": row["test_cas"]})
    return pd.DataFrame(registros)