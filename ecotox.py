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

    # conc1_mean vem como string no ECOTOX (pode conter 'NR', '~10', '<5', etc.)
    # — converter para numérico, descartando valores não parseáveis.
    df["conc1_mean"] = pd.to_numeric(df["conc1_mean"], errors="coerce")
    df = df.dropna(subset=["conc1_mean"])
    df = df[df["conc1_mean"] > 0]

    return df.drop_duplicates()


def _formatar_cas(cas_digits: str) -> str:
    """Converte CAS sem hífens (como armazenado no ECOTOX) para o formato
    padrão com hífens: RRRRRR-DD-C (ex: '108952' → '108-95-2').
    Se já tiver hífens, retorna sem modificar."""
    s = str(cas_digits).strip()
    if "-" in s:
        return s
    if len(s) >= 3:
        return f"{s[:-3]}-{s[-3:-1]}-{s[-1]}"
    return s


def cas_para_smiles(cas: str, cache: dict) -> str | None:
    if cas in cache:
        return cache[cas]

    cas_fmt = _formatar_cas(cas)
    # /name/ endpoint resolve CAS formatado de forma unívoca e retorna
    # CanonicalSMILES (ao contrário de /xref/RegistryID/ que pode retornar
    # múltiplos CIDs ambíguos com ConnectivitySMILES).
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{cas_fmt}/property/CanonicalSMILES/JSON")
    try:
        resp = requests.get(url, timeout=10)
        time.sleep(0.25)  # PubChem pede no máx. ~5 requisições/segundo
        if resp.status_code == 200:
            props = resp.json()["PropertyTable"]["Properties"]
            smiles = (props[0].get("CanonicalSMILES")
                      or props[0].get("ConnectivitySMILES"))
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

    # Se o cache existir mas todos os valores forem None (corrompido),
    # descartar e re-consultar tudo.
    if cache and all(v is None for v in cache.values()):
        print("[AVISO] Cache do PubChem está corrompido (todos None) — descartando e re-consultando.")
        cache = {}

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
    n_sem_smiles = df["smiles"].isna().sum()
    if n_sem_smiles > 0:
        print(f"[INFO] {n_sem_smiles} linhas sem SMILES resolvido (CAS não encontrado no PubChem).")
    return df.dropna(subset=["smiles"])


# Fatores de conversão de conc1_unit para mol/L.
# ECOTOX usa muitas unidades; as mais comuns para ensaios com algas são:
#   mg/L  → divide por 1000 → g/L → divide por MolWt → mol/L
#   ug/L  → divide por 1e6  → g/L → divide por MolWt → mol/L
#   ng/L  → divide por 1e9  → g/L → divide por MolWt → mol/L
#   ppb   ≈ ug/L em água (densidade ≈ 1 g/mL)
#   ppm   ≈ mg/L em água
#   mol/L → já está em mol/L (não precisa de MolWt)
#   mmol/L, umol/L, nmol/L → factor direto
_UNIT_TO_G_PER_L: dict[str, float] = {
    # mg/L e variantes
    "mg/l":     1e-3,
    "mg/L":     1e-3,
    "mg/dm3":   1e-3,   # dm³ = L
    "mg/dm³":   1e-3,
    "ppm":      1e-3,   # mg/L ≈ ppm em água
    # AI mg/L e ae mg/L: "active ingredient" / "acid equivalent" — mesma unidade de massa
    "ai mg/l":  1e-3,
    "ai mg/L":  1e-3,
    "ae mg/l":  1e-3,
    "ae mg/L":  1e-3,
    # ug/L e variantes
    "ug/l":     1e-6,
    "ug/L":     1e-6,
    "ug/ml":    1e-3,   # ug/mL = mg/L
    "ug/mL":    1e-3,
    "ppb":      1e-6,   # ug/L ≈ ppb em água
    "ai ug/l":  1e-6,
    "ai ug/L":  1e-6,
    # ng/L
    "ng/l":     1e-9,
    "ng/L":     1e-9,
    # g/L
    "g/l":      1.0,
    "g/L":      1.0,
}
_UNIT_TO_MOL_PER_L: dict[str, float] = {
    "mol/l":    1.0,
    "mol/L":    1.0,
    # mmol
    "mmol/l":   1e-3,
    "mmol/L":   1e-3,
    "mm":       1e-3,   # mM (millimolar)
    "mm3/l":    1e-3,   # provavelmente mM mal formatado no ECOTOX
    "mm3/L":    1e-3,
    "mM":       1e-3,
    # umol / uM
    "umol/l":   1e-6,
    "umol/L":   1e-6,
    "umol/dm3": 1e-6,
    "umol/dm³": 1e-6,
    "um":       1e-6,
    "uM":       1e-6,
    "um/l":     1e-6,
    "uM/L":     1e-6,
    # nmol / nM
    "nmol/l":   1e-9,
    "nmol/L":   1e-9,
    "nm":       1e-9,
    "nM":       1e-9,
}


def _conc_para_mol_por_l(conc_mean: float, unit_raw: str, mol_wt: float) -> float | None:
    """Converte conc1_mean (em unidade unit_raw) para mol/L.
    Retorna None se a unidade for desconhecida ou se o resultado for inválido."""
    unit = str(unit_raw).strip()
    if unit in _UNIT_TO_MOL_PER_L:
        return conc_mean * _UNIT_TO_MOL_PER_L[unit]
    if unit in _UNIT_TO_G_PER_L:
        if mol_wt is None or mol_wt <= 0:
            return None
        return (conc_mean * _UNIT_TO_G_PER_L[unit]) / mol_wt
    return None  # unidade desconhecida — descartar


def montar_matriz_treino(df: pd.DataFrame) -> pd.DataFrame:
    """Junta descritores + variável-resposta em escala molar (pEC50)."""
    registros = []
    skipped_unit = 0
    skipped_invalid = 0
    for _, row in df.iterrows():
        desc = calcular_descritores(row["smiles"])
        if desc is None:
            skipped_invalid += 1
            continue

        unit_raw = row.get("conc1_unit", "mg/L")  # fallback se coluna ausente
        conc_molar = _conc_para_mol_por_l(
            float(row["conc1_mean"]), unit_raw, desc["MolWt"]
        )
        if conc_molar is None:
            skipped_unit += 1
            continue
        if conc_molar <= 0:
            skipped_invalid += 1
            continue

        registros.append({**desc, "pEC50": -np.log10(conc_molar), "cas": row["test_cas"]})

    if skipped_unit > 0:
        print(f"[AVISO] {skipped_unit} linhas descartadas por unidade de concentração "
              f"desconhecida (adicionar em _UNIT_TO_G_PER_L/_UNIT_TO_MOL_PER_L se necessário).")
    if skipped_invalid > 0:
        print(f"[INFO] {skipped_invalid} linhas descartadas por SMILES inválido ou "
              f"concentração ≤ 0.")

    matriz = pd.DataFrame(registros)
    if matriz.empty or "pEC50" not in matriz.columns:
        raise ValueError(
            "A matriz de treino ficou vazia. Verifique:\n"
            "  1. Se o arquivo ECOTOX está correto e contém linhas para "
            f"'{ESPECIE_ALVO}'.\n"
            "  2. Se as unidades de concentração em conc1_unit estão mapeadas "
            "em _UNIT_TO_G_PER_L ou _UNIT_TO_MOL_PER_L.\n"
            "  3. Se o PubChem retornou SMILES válidos (verifique pubchem_cache.json)."
        )
    print(f"[INFO] Matriz de treino montada: {len(matriz)} amostras, "
          f"{matriz['pEC50'].isna().sum()} pEC50 NaN.")
    return matriz