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
    from sklearn.svm import SVR
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import KFold, GridSearchCV

    colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
    X, y = matriz[colunas_x], matriz["pEC50"]

    print("\n[INFO] Padronizando as variáveis (StandardScaler)...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    # Definir modelos e grades de parâmetros para otimização
    modelos = {
        "RandomForest": {
            "estimator": RandomForestRegressor(random_state=42, n_jobs=2),
            "params": {
                "n_estimators": [100, 300],
                "max_depth": [None, 10, 20],
                "min_samples_leaf": [1, 2]
            }
        },
        "SVR": {
            "estimator": SVR(),
            "params": {
                "C": [0.1, 1, 10],
                "gamma": ["scale", "auto", 0.01],
                "epsilon": [0.01, 0.1, 0.2]
            }
        }
    }

    melhor_modelo = None
    melhor_nome = ""
    melhor_score = -float('inf')

    print("[INFO] Iniciando Grid Search (Busca em Grade)... Isso pode demorar um pouco.")
    for nome, config in modelos.items():
        print(f"  -> Otimizando {nome}...")
        grid = GridSearchCV(config["estimator"], config["params"], cv=kf, scoring="r2", n_jobs=2)
        grid.fit(X_scaled, y)
        
        score_cv = grid.best_score_
        print(f"     Melhor R² CV: {score_cv:.3f} | Params: {grid.best_params_}")
        
        if score_cv > melhor_score:
            melhor_score = score_cv
            melhor_modelo = grid.best_estimator_
            melhor_nome = nome

    print(f"\n[INFO] Modelo vencedor: {melhor_nome} (R² CV = {melhor_score:.3f})")

    # Treina o modelo final com toda a base pública (o GridSearchCV já faz isso com best_estimator_, mas fica explícito)
    melhor_modelo.fit(X_scaled, y)
    return melhor_modelo, colunas_x, scaler


def treinar_classificador_ghs(matriz: pd.DataFrame):
    """Treina um classificador RandomForest especificamente para as categorias
    GHS, utilizando class_weight='balanced' para lidar com o desbalanceamento
    das classes, em especial a classe 'Não classificado (>100 mg/L)'."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import KFold, cross_val_predict
    from validacao import classificar_toxicidade_ghs, LABELS_GHS

    colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
    X = matriz[colunas_x]
    
    # Criar a coluna de classes GHS alvo
    y_class = [classificar_toxicidade_ghs(row["pEC50"], row["MolWt"]) for _, row in matriz.iterrows()]
    y_class = pd.Series(y_class)

    # n_jobs=2 para paralelizar de leve, class_weight="balanced" para compensar o desbalanceamento
    modelo = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2, class_weight="balanced")
    
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    pred_cv = cross_val_predict(modelo, X, y_class, cv=kf)
    
    modelo.fit(X, y_class)
    return modelo, pred_cv, y_class


def importancia_descritores(modelo, colunas_x: list[str], scaler=None, matriz=None) -> pd.DataFrame:
    """Ranking de quais descritores mais pesam na predição.
    Usa Permutation Importance (agnóstico ao modelo) para suportar tanto RF quanto SVR."""
    from sklearn.inspection import permutation_importance

    if matriz is not None and scaler is not None:
        X = matriz[colunas_x]
        y = matriz["pEC50"]
        X_scaled = scaler.transform(X)
        resultado = permutation_importance(modelo, X_scaled, y, n_repeats=10, random_state=42, n_jobs=2)
        importancias = resultado.importances_mean
    elif hasattr(modelo, "feature_importances_"):
        # Fallback para classificador ou RF legado
        importancias = modelo.feature_importances_
    else:
        importancias = [0] * len(colunas_x)

    ranking = pd.DataFrame({
        "descritor": colunas_x,
        "importancia": importancias,
    }).sort_values("importancia", ascending=False).reset_index(drop=True)
    return ranking


def validar_externamente(modelo, colunas_x, caminho_planilha: str, scaler=None) -> pd.DataFrame:
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
        if scaler is not None:
            X_novo_scaled = scaler.transform(X_novo)
        else:
            X_novo_scaled = X_novo
        pred_pEC50 = modelo.predict(X_novo_scaled)[0]
        linhas.append({"Ingrediente": row["Ingrediente"], "pEC50_previsto": pred_pEC50})

    return pd.DataFrame(linhas)