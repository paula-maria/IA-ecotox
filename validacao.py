"""Validação em duas frentes:

  A. Os DADOS estão corretos? — duplicatas, valores ausentes, outliers
     estatísticos, SMILES que não bateram com o CAS esperado.
  B. O MODELO está prevendo de verdade? — y-scrambling (embaralhamento do
     alvo) e domínio de aplicabilidade (leverage), os dois testes centrais
     dos Princípios OECD de validação de QSAR.

Ver o checklist completo (5 itens) no chat — este módulo implementa os
itens que dá pra automatizar (1, 2 e 3); os itens 4 (comparar com
ECOSAR/VEGA) e 5 (plausibilidade biológica) são análise manual, mas a
função `importancia_descritores` do modelo.py já dá o ranking pra isso.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# A. QUALIDADE DOS DADOS
# ---------------------------------------------------------------------------
def checar_qualidade_dados(matriz: pd.DataFrame, coluna_y: str = "pEC50") -> dict:
    """Roda antes de treinar. Retorna um relatório com o que encontrar —
    não corrige nada sozinho, porque decidir o que fazer com um outlier
    (remover? manter? investigar?) é uma decisão que precisa de critério
    científico, não deve ser automática."""
    relatorio = {}

    relatorio["n_linhas"] = len(matriz)
    relatorio["duplicatas"] = int(matriz.duplicated().sum())
    relatorio["valores_ausentes"] = matriz.isna().sum().to_dict()

    # outliers estatísticos no alvo (pEC50), método IQR — só sinaliza, não remove
    y = matriz[coluna_y].dropna()
    q1, q3 = y.quantile([0.25, 0.75])
    iqr = q3 - q1
    limite_inf, limite_sup = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    outliers = matriz[(matriz[coluna_y] < limite_inf) | (matriz[coluna_y] > limite_sup)]
    relatorio["outliers_pEC50"] = {
        "n": len(outliers),
        "cas": outliers["cas"].tolist() if "cas" in outliers.columns else None,
        "limites": (round(limite_inf, 2), round(limite_sup, 2)),
    }

    # descritores com variância zero não ajudam o modelo e podem indicar
    # erro de cálculo (ex.: todos os SMILES caindo no mesmo valor)
    colunas_numericas = matriz.select_dtypes("number").columns
    variancia_zero = [c for c in colunas_numericas if matriz[c].nunique() <= 1]
    relatorio["descritores_sem_variacao"] = variancia_zero

    return relatorio


def imprimir_relatorio_qualidade(relatorio: dict):
    print("=== Qualidade dos dados ===")
    print(f"Linhas: {relatorio['n_linhas']}  |  Duplicatas: {relatorio['duplicatas']}")
    if any(relatorio["valores_ausentes"].values()):
        print(f"Valores ausentes: {relatorio['valores_ausentes']}")
    n_out = relatorio["outliers_pEC50"]["n"]
    if n_out:
        print(f"[ATENÇÃO] {n_out} outlier(s) em pEC50 fora de "
              f"{relatorio['outliers_pEC50']['limites']} — investigar antes de treinar.")
    if relatorio["descritores_sem_variacao"]:
        print(f"[ATENÇÃO] Descritores sem variação (suspeitos): "
              f"{relatorio['descritores_sem_variacao']}")


# ---------------------------------------------------------------------------
# B.1 — Y-SCRAMBLING (o modelo aprendeu sinal real ou está decorando ruído?)
# ---------------------------------------------------------------------------
def teste_y_scrambling(X: pd.DataFrame, y: pd.Series, n_repeticoes: int = 30,
                        random_state: int = 42) -> dict:
    """Embaralha o alvo (pEC50) aleatoriamente e retreina o modelo várias
    vezes. Se o R² embaralhado ficar perto do R² real, o modelo não está
    aprendendo relação estrutura-atividade nenhuma — é falso positivo
    estatístico. Espera-se que R²_real seja MUITO maior que a média dos
    R²_embaralhados."""
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.metrics import r2_score

    rng = np.random.RandomState(random_state)
    kf = KFold(n_splits=5, shuffle=True, random_state=random_state)

    # n_jobs=-1: usa todos os núcleos de CPU disponíveis para paralelizar e acelerar o processamento
    modelo_real = RandomForestRegressor(n_estimators=300, random_state=random_state, n_jobs=-1)
    pred_real = cross_val_predict(modelo_real, X, y, cv=kf)
    r2_real = r2_score(y, pred_real)

    r2_embaralhados = []
    for _ in range(n_repeticoes):
        y_embaralhado = y.sample(frac=1.0, random_state=rng.randint(0, 1_000_000)).reset_index(drop=True)
        # n_jobs=-1: paraleliza o treino de cada repetição do embaralhado
        modelo_emb = RandomForestRegressor(n_estimators=300, random_state=random_state, n_jobs=-1)
        pred_emb = cross_val_predict(modelo_emb, X, y_embaralhado, cv=kf)
        r2_embaralhados.append(r2_score(y_embaralhado, pred_emb))

    media_embaralhado = float(np.mean(r2_embaralhados))
    return {
        "R2_real": round(r2_real, 3),
        "R2_embaralhado_media": round(media_embaralhado, 3),
        "R2_embaralhado_max": round(max(r2_embaralhados), 3),
        "passou_no_teste": r2_real > (media_embaralhado + 2 * np.std(r2_embaralhados)),
    }


# ---------------------------------------------------------------------------
# B.2 — DOMÍNIO DE APLICABILIDADE (leverage) para os ativos amazônicos
# ---------------------------------------------------------------------------
def calcular_leverage(X_treino: pd.DataFrame, X_novo: pd.DataFrame) -> pd.DataFrame:
    """Leverage (h) de cada composto novo em relação ao espaço de treino,
    pela matriz-chapéu clássica de QSAR: h_i = x_i (X'X)^-1 x_i'.
    Limite de corte usual: h* = 3(p+1)/n, onde p = nº de descritores e
    n = nº de compostos de treino. Acima de h*, a predição está fora do
    domínio de aplicabilidade — reportar como pouco confiável, não
    descartar silenciosamente."""
    X_t = X_treino.to_numpy()
    XtX_inv = np.linalg.pinv(X_t.T @ X_t)  # pinv: robusto a colinearidade

    p = X_treino.shape[1]
    n = X_treino.shape[0]
    h_estrela = 3 * (p + 1) / n

    leverages = []
    for _, linha in X_novo.iterrows():
        x = linha.to_numpy()
        h = float(x @ XtX_inv @ x.T)
        leverages.append(h)

    resultado = X_novo.copy()
    resultado["leverage"] = leverages
    resultado["dentro_do_dominio"] = resultado["leverage"] <= h_estrela
    resultado.attrs["h_estrela"] = h_estrela
    return resultado


# ---------------------------------------------------------------------------
# C. MATRIZ DE CONFUSÃO — pEC50 (contínuo) → categoria GHS → comparação visual
# ---------------------------------------------------------------------------
LABELS_GHS = ["Categoria 1 (≤1 mg/L)", "Categoria 2 (1–10 mg/L)",
              "Categoria 3 (10–100 mg/L)", "Não classificado (>100 mg/L)"]


def classificar_toxicidade_ghs(pEC50: float, mol_wt: float) -> str:
    """Converte pEC50 (escala molar) de volta para mg/L usando o peso
    molecular, e classifica pelas faixas de toxicidade aguda aquática do
    GHS (Sistema Globalmente Harmonizado de classificação de produtos
    químicos, ONU) — referencial reconhecido, não uma categorização ad hoc."""
    ec50_molar = 10 ** (-pEC50)
    ec50_mg_l = ec50_molar * mol_wt * 1000

    if ec50_mg_l <= 1:
        return LABELS_GHS[0]
    elif ec50_mg_l <= 10:
        return LABELS_GHS[1]
    elif ec50_mg_l <= 100:
        return LABELS_GHS[2]
    return LABELS_GHS[3]


def matriz_confusao_toxicidade(pEC50_observado: pd.Series, pEC50_previsto: pd.Series,
                                mol_wt: pd.Series):
    """Monta a matriz de confusão (observado × previsto) nas categorias
    GHS. Retorna a matriz como DataFrame (fácil de imprimir/exportar) e os
    rótulos usados. Requer scikit-learn."""
    from sklearn.metrics import confusion_matrix

    classes_obs = [classificar_toxicidade_ghs(p, m) for p, m in zip(pEC50_observado, mol_wt)]
    classes_prev = [classificar_toxicidade_ghs(p, m) for p, m in zip(pEC50_previsto, mol_wt)]

    matriz = confusion_matrix(classes_obs, classes_prev, labels=LABELS_GHS)
    df_matriz = pd.DataFrame(matriz, index=LABELS_GHS, columns=LABELS_GHS)
    df_matriz.index.name = "Observado"
    df_matriz.columns.name = "Previsto"
    return df_matriz


def plotar_matriz_confusao(df_matriz: pd.DataFrame, caminho_saida: str = "matriz_confusao.png"):
    """Gera um heatmap da matriz de confusão e salva como PNG."""
   

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(df_matriz.values, cmap="Blues")

    ax.set_xticks(range(len(df_matriz.columns)))
    ax.set_yticks(range(len(df_matriz.index)))
    ax.set_xticklabels(df_matriz.columns, rotation=30, ha="right")
    ax.set_yticklabels(df_matriz.index)
    ax.set_xlabel("Previsto pelo modelo")
    ax.set_ylabel("Observado (experimental)")
    ax.set_title("Matriz de confusão — categorias de toxicidade (GHS)")

    for i in range(df_matriz.shape[0]):
        for j in range(df_matriz.shape[1]):
            valor = df_matriz.values[i, j]
            cor_texto = "white" if valor > df_matriz.values.max() / 2 else "black"
            ax.text(j, i, str(valor), ha="center", va="center", color=cor_texto)

    fig.colorbar(im, ax=ax, label="Nº de compostos")
    fig.tight_layout()
    fig.savefig(caminho_saida, dpi=150)
    plt.close(fig)
    return caminho_saida