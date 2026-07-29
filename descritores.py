"""Cálculo de descritores moleculares a partir de SMILES (RDKit).

Usado tanto pela Fase 1 (treino com o ECOTOX, em ecotox.py) quanto pela
Fase 2 (validação externa nos ativos amazônicos, em dados.py/modelo.py) —
é o mesmo conjunto de features nas duas fases, por isso mora em um único
lugar.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def calcular_descritores(smiles: str) -> dict | None:
    """Retorna um dicionário de descritores QSAR clássicos e MACCS Keys para um SMILES.
    Retorna None se o SMILES for inválido — nesse caso, revisar manualmente
    (RDKit aceita SMILES quimicamente válidos mesmo que sejam o composto
    errado; veja o aviso sobre o 'Decyl glucoside' no README)."""
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    
    descritores = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),          # lipofilicidade — forte preditor em ecotox aquática
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }

    # Adiciona MACCS Keys (166 features binárias, 1 a 166)
    from rdkit.Chem import MACCSkeys
    maccs = MACCSkeys.GenMACCSKeys(mol)
    for i in range(1, 167):
        descritores[f"MACCS_{i}"] = int(maccs[i])
        
    return descritores