"""Ajuste de EC50, treino do modelo com a base pública, modelo exploratório
com os dados próprios, e aplicação do modelo (validação externa) aos
ativos amazônicos.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from descritores import calcular_descritores


# ---------------------------------------------------------------------------
# EC50 POR DOSE-RESPOSTA (dados próprios, quando houver >=3 concentrações)
# ---------------------------------------------------------------------------
def _modelo_logistico(conc, top, bottom, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / conc) ** hill)


def ajustar_ec50(concentracoes, inibicoes) -> dict | None:
    """Ajusta uma curva log-logística simples. Precisa de >=3 pontos de
    concentração distintos e >0 (ex.: 100 %, 1 %, 0,1 % da formulação)."""
    conc = np.asarray(concentracoes, dtype=float)
    inib = np.asarray(inibicoes, dtype=float)
    try:
        popt, _ = curve_fit(
            _modelo_logistico, conc, inib,
            p0=[100, 0, np.median(conc), 1],
            bounds=([0, 0, 1e-6, 0.1], [100, 50, 1e6, 10]),
            maxfev=5000,
        )
        return {"top": popt[0], "bottom": popt[1], "EC50": popt[2], "hill": popt[3]}
    except RuntimeError:
        return None


# ---------------------------------------------------------------------------
# MODELO EXPLORATÓRIO (n pequeno — dados próprios, sem a base pública)
# ---------------------------------------------------------------------------
def treinar_modelo_exploratorio(X: pd.DataFrame, y: pd.Series) -> dict:
    """Com n ~ 6 formulações, usar regressão linear regularizada (Ridge) com
    leave-one-out em vez de random forest / rede neural, que decoraria os
    dados em vez de generalizar."""
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler

    loo = LeaveOneOut()
    erros = []
    for train_idx, test_idx in loo.split(X):
        scaler = StandardScaler().fit(X.iloc[train_idx])
        X_train = scaler.transform(X.iloc[train_idx])
        X_test = scaler.transform(X.iloc[test_idx])

        modelo = Ridge(alpha=10.0).fit(X_train, y.iloc[train_idx])
        pred = modelo.predict(X_test)
        erros.append(abs(pred[0] - y.iloc[test_idx].values[0]))

    return {"MAE_leave_one_out": np.mean(erros), "n_amostras": len(X)}


# ---------------------------------------------------------------------------
# TREINO NA BASE PÚBLICA (Fase 1) + VALIDAÇÃO EXTERNA (Fase 2)
# ---------------------------------------------------------------------------
def treinar_modelo_publico(matriz: pd.DataFrame):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.metrics import r2_score, mean_squared_error

    colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
    X, y = matriz[colunas_x], matriz["pEC50"]

    modelo = RandomForestRegressor(n_estimators=300, random_state=42)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    pred_cv = cross_val_predict(modelo, X, y, cv=kf)

    print(f"R² (validação cruzada 5-fold): {r2_score(y, pred_cv):.3f}")
    print(f"RMSE (validação cruzada 5-fold): {mean_squared_error(y, pred_cv, squared=False):.3f}")

    modelo.fit(X, y)  # modelo final treinado com toda a base pública
    return modelo, colunas_x


def validar_externamente(modelo, colunas_x, caminho_planilha: str) -> pd.DataFrame:
    """Aplica o modelo (sem retreino) aos ingredientes isolados da planilha
    própria — validação externa da Fase 2."""
    ingredientes = pd.read_excel(caminho_planilha,
                                  sheet_name="5_SMILES_Ingredientes", header=1)
    ingredientes = ingredientes.dropna(subset=["Ingrediente"])

    linhas = []
    for _, row in ingredientes.iterrows():
        desc = calcular_descritores(row["SMILES_canonical"])
        if desc is None:
            continue
        X_novo = pd.DataFrame([desc])[colunas_x]
        pred_pEC50 = modelo.predict(X_novo)[0]
        linhas.append({"Ingrediente": row["Ingrediente"], "pEC50_previsto": pred_pEC50})

    return pd.DataFrame(linhas)