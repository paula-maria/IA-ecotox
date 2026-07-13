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
    path_tests = diretorio / "tests.txt"
    path_results = diretorio / "results.txt"
    path_species = diretorio / "species.txt"

    if not path_tests.exists():
        candidatos = list(diretorio.glob("**/tests.txt"))
        if candidatos:
            path_tests = candidatos[0]
            base_dir = path_tests.parent
            path_results = base_dir / "results.txt"
            path_species = base_dir / "species.txt"
            if not path_species.exists():
                candidatos_species = list(base_dir.glob("**/species.txt"))
                if candidatos_species:
                    path_species = candidatos_species[0]
                else:
                    candidatos_species = list(diretorio.glob("**/species.txt"))
                    if candidatos_species:
                        path_species = candidatos_species[0]
        else:
            raise FileNotFoundError(f"Não foi possível encontrar tests.txt dentro de {diretorio}")

    if not path_results.exists():
        raise FileNotFoundError(f"Não foi possível encontrar results.txt (procurado em {path_results})")
    if not path_species.exists():
        raise FileNotFoundError(f"Não foi possível encontrar species.txt (procurado em {path_species})")

    # Carrega species.txt filtrando apenas colunas necessárias
    species = pd.read_csv(path_species, sep="|",
                          usecols=lambda c: c.strip().lower() in ["species_number", "latin_name"],
                          low_memory=False, encoding="latin1", on_bad_lines="skip")
    species.columns = [c.strip().lower() for c in species.columns]

    species_alga = species[species["latin_name"].str.strip().str.lower()
                            == especie.lower()]
    origem_filtro = f"espécie exata ({especie})"

    # Carrega tests.txt filtrando apenas colunas necessárias
    tests = pd.read_csv(path_tests, sep="|",
                        usecols=lambda c: c.strip().lower() in ["test_id", "test_cas", "species_number"],
                        low_memory=False, encoding="latin1", on_bad_lines="skip")
    tests.columns = [c.strip().lower() for c in tests.columns]

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

    # Carrega results.txt filtrando apenas colunas necessárias
    results = pd.read_csv(path_results, sep="|",
                          usecols=lambda c: c.strip().lower() in ["result_id", "test_id", "endpoint", "effect", "conc1_mean", "conc1_unit"],
                          low_memory=False, encoding="latin1", on_bad_lines="skip")
    results.columns = [c.strip().lower() for c in results.columns]

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
    import json
    cache_file = Path("pubchem_cache.json")
    cache = {}
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                cache = json.load(f)
        except Exception as e:
            print(f"[AVISO] Erro ao ler cache do PubChem: {e}")

    df = df.copy()
    unique_cas = df["test_cas"].unique()
    new_cas = [c for c in unique_cas if c not in cache]

    if new_cas:
        print(f"[INFO] Buscando SMILES para {len(new_cas)} novos CAS de {len(unique_cas)} no PubChem...")
        resolved_count = 0
        total_new = len(new_cas)
        for c in new_cas:
            cas_para_smiles(c, cache)
            resolved_count += 1
            if resolved_count % 50 == 0 or resolved_count == total_new:
                print(f"  -> Progresso: {resolved_count}/{total_new} novos CAS consultados...")
                try:
                    with open(cache_file, "w") as f:
                        json.dump(cache, f, indent=4)
                except Exception:
                    pass
    else:
        print(f"[INFO] Todos os {len(unique_cas)} CAS já estão no cache local ({cache_file.name}).")

    df["smiles"] = df["test_cas"].map(cache)
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