"""
Script de validação completo — roda uma vez e gera todos os resultados.
Evita o duplo GridSearch do main.py/validar para economizar tempo.
"""
import sys, time
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import ecotox, dados, modelo, validacao
from validacao import classificar_toxicidade_ghs, LABELS_GHS

start = time.time()

# ─────────────────────────────────────────────────────────────────
# 1. Carrega ECOTOX e monta matriz (deduplicada por CAS)
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 1 — Carregando ECOTOX e montando matriz de treino")
print("="*60)
dados_ecotox = ecotox.carregar_ecotox_algas()
dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
matriz = ecotox.montar_matriz_treino(dados_ecotox)

colunas_x = [c for c in matriz.columns if c not in ("pEC50", "cas")]
X, y = matriz[colunas_x], matriz["pEC50"]
print(f"\n  → Matriz final: {len(matriz)} compostos únicos, {len(colunas_x)} descritores")

# ─────────────────────────────────────────────────────────────────
# 2. Qualidade dos dados
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 2 — Qualidade dos dados")
print("="*60)
rel = validacao.checar_qualidade_dados(matriz)
validacao.imprimir_relatorio_qualidade(rel)

print(f"\n  pEC50 — min: {y.min():.2f} | max: {y.max():.2f} | média: {y.mean():.2f} | DP: {y.std():.2f}")
print(f"  Descritores sem variação: {len(rel['descritores_sem_variacao'])} → {rel['descritores_sem_variacao'][:5]}")

# Distribuição GHS
classes_obs = [classificar_toxicidade_ghs(row["pEC50"], row["MolWt"]) for _, row in matriz.iterrows()]
contagem_ghs = pd.Series(classes_obs).value_counts()
print("\n  Distribuição de Classes GHS (dados de treino):")
for label, n in contagem_ghs.items():
    print(f"    {label}: {n} compostos ({100*n/len(matriz):.1f}%)")

# ─────────────────────────────────────────────────────────────────
# 3. Treina o melhor modelo (GridSearchCV) UMA única vez
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 3 — Treinamento (GridSearchCV — Random Forest + SVR)")
print("="*60)
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.metrics import r2_score, mean_squared_error

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
kf = KFold(n_splits=5, shuffle=True, random_state=42)

modelos_grid = {
    "RandomForest": {
        "estimator": RandomForestRegressor(random_state=42, n_jobs=2),
        "params": {"n_estimators": [100, 300], "max_depth": [None, 10, 20], "min_samples_leaf": [1, 2]}
    },
    "SVR": {
        "estimator": SVR(),
        "params": {"C": [0.1, 1, 10], "gamma": ["scale", "auto", 0.01], "epsilon": [0.01, 0.1, 0.2]}
    }
}

melhor_modelo = None
melhor_nome = ""
melhor_r2_cv = -float("inf")
todos_r2 = {}

for nome, cfg in modelos_grid.items():
    print(f"  → Otimizando {nome}...")
    grid = GridSearchCV(cfg["estimator"], cfg["params"], cv=kf, scoring="r2", n_jobs=2)
    grid.fit(X_scaled, y)
    todos_r2[nome] = round(grid.best_score_, 3)
    print(f"     Melhor R² CV: {grid.best_score_:.3f} | Params: {grid.best_params_}")
    if grid.best_score_ > melhor_r2_cv:
        melhor_r2_cv = grid.best_score_
        melhor_modelo = grid.best_estimator_
        melhor_nome = nome

print(f"\n  ★ Modelo vencedor: {melhor_nome} (R² CV = {melhor_r2_cv:.3f})")
melhor_modelo.fit(X_scaled, y)

# Predições por CV cruzado para gráficos de resíduos
pred_cv = cross_val_predict(melhor_modelo, X_scaled, y, cv=kf)
rmse_cv = float(np.sqrt(mean_squared_error(y, pred_cv)))
r2_cv   = float(r2_score(y, pred_cv))
print(f"  R² (cross_val_predict): {r2_cv:.3f}")
print(f"  RMSE (cross_val_predict): {rmse_cv:.3f}")

# ─────────────────────────────────────────────────────────────────
# 4. Y-Scrambling
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 4 — Y-Scrambling (30 repetições)")
print("="*60)
resultado_scr = validacao.teste_y_scrambling(X, y, n_repeticoes=30)
print(f"  R² real (CV):              {resultado_scr['R2_real']:.3f}")
print(f"  R² embaralhado (média):    {resultado_scr['R2_embaralhado_media']:.3f}")
print(f"  R² embaralhado (máx):      {resultado_scr['R2_embaralhado_max']:.3f}")
delta = resultado_scr["R2_real"] - resultado_scr["R2_embaralhado_media"]
print(f"  Δ R² (real − embaralhado): {delta:.3f}")
print(f"  Passou no teste?           {'✅ SIM' if resultado_scr['passou_no_teste'] else '❌ NÃO'}")
if not resultado_scr["passou_no_teste"]:
    print("  [ATENÇÃO] Modelo não se distancia suficientemente do embaralhado!")

# ─────────────────────────────────────────────────────────────────
# 5. Domínio de aplicabilidade (leverage)
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 5 — Domínio de Aplicabilidade (Leverage Williams Plot)")
print("="*60)
ingredientes = dados.carregar_descritores_ingredientes()
colunas_comuns = [c for c in colunas_x if c in ingredientes.columns]
leverage_df = validacao.calcular_leverage(X[colunas_comuns], ingredientes[colunas_comuns])
leverage_df["Ingrediente"] = ingredientes["Ingrediente"].values
h_star = leverage_df.attrs["h_estrela"]

print(f"  h* (limite de corte) = {h_star:.4f}")
print(f"  nº descritores (p): {len(colunas_comuns)}  |  nº compostos treino (n): {len(X)}")
print()
dentro  = leverage_df["dentro_do_dominio"].sum()
fora    = (~leverage_df["dentro_do_dominio"]).sum()
print(f"  Dentro do domínio: {dentro}/{len(leverage_df)} ingredientes")
print(f"  Fora do domínio:   {fora}/{len(leverage_df)} ingredientes")
print()
print(f"  {'Ingrediente':<40} {'Leverage':>10} {'h ≤ h*':>12} {'Status':>12}")
print("  " + "-"*76)
for _, row in leverage_df.sort_values("leverage").iterrows():
    status = "✅ Dentro" if row["dentro_do_dominio"] else "⚠️  FORA"
    print(f"  {row['Ingrediente']:<40} {row['leverage']:>10.4f} {str(row['dentro_do_dominio']):>12} {status:>12}")

# ─────────────────────────────────────────────────────────────────
# 6. Importância dos descritores (top 15)
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 6 — Importância dos descritores (Permutation Importance)")
print("="*60)
ranking = modelo.importancia_descritores(melhor_modelo, colunas_x, scaler, matriz)
print(f"\n  Top 15 descritores mais preditivos:")
print(f"  {'#':<4} {'Descritor':<35} {'Importância':>12}")
print("  " + "-"*53)
for i, row in ranking.head(15).iterrows():
    print(f"  {i+1:<4} {row['descritor']:<35} {row['importancia']:>12.4f}")

# ─────────────────────────────────────────────────────────────────
# 7. Matriz de confusão GHS (regressão)
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 7 — Matriz de confusão GHS (regressão → CV)")
print("="*60)
df_conf = validacao.matriz_confusao_toxicidade(y, pd.Series(pred_cv), matriz["MolWt"])
print("\n  (linhas = Observado | colunas = Previsto pelo modelo)\n")
print(df_conf.to_string())

# Acurácia por classe
print("\n  Acertos por classe (diagonal):")
total = df_conf.values.sum()
for label in LABELS_GHS:
    acertos = df_conf.loc[label, label]
    total_obs = df_conf.loc[label].sum()
    pct = 100 * acertos / total_obs if total_obs > 0 else 0
    print(f"    {label}: {acertos}/{total_obs} ({pct:.0f}%)")
print(f"  Acurácia global: {100*np.diag(df_conf.values).sum()/total:.1f}%")

# Salva gráfico
caminho_conf = validacao.plotar_matriz_confusao(df_conf, "matriz_confusao_regressao.png")
caminho_res  = validacao.plotar_analise_residuos(y, pd.Series(pred_cv), "analise_residuos.png")
print(f"\n  Gráficos salvos: {caminho_conf} | {caminho_res}")

# ─────────────────────────────────────────────────────────────────
# 8. Previsão dos ativos amazônicos
# ─────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("ETAPA 8 — Previsão dos ativos amazônicos (validação externa)")
print("="*60)
previsao = modelo.validar_externamente(melhor_modelo, colunas_x, dados.ARQUIVO_PADRAO, scaler)

# Adiciona classificação GHS e leverage
previsao = previsao.merge(
    leverage_df[["Ingrediente", "leverage", "dentro_do_dominio"]], on="Ingrediente", how="left"
)

# Converte pEC50 → EC50 mg/L para facilitar interpretação
desc_ing = ingredientes.set_index("Ingrediente")
previsao["EC50_mgL"] = previsao.apply(
    lambda r: round(10**(-r["pEC50_previsto"]) * desc_ing.loc[r["Ingrediente"], "MolWt"] * 1000, 2)
    if r["Ingrediente"] in desc_ing.index else None,
    axis=1
)
previsao["Classe_GHS"] = previsao.apply(
    lambda r: classificar_toxicidade_ghs(r["pEC50_previsto"], desc_ing.loc[r["Ingrediente"], "MolWt"])
    if r["Ingrediente"] in desc_ing.index else None,
    axis=1
)

print(f"\n  {'Ingrediente':<40} {'pEC50':>7} {'EC50 mg/L':>12} {'Classe GHS':<35} {'Domínio':>10}")
print("  " + "-"*106)
for _, r in previsao.sort_values("pEC50_previsto", ascending=False).iterrows():
    dom = "✅" if r.get("dentro_do_dominio") else "⚠️ fora"
    print(f"  {r['Ingrediente']:<40} {r['pEC50_previsto']:>7.2f} {str(r['EC50_mgL']):>12} {str(r['Classe_GHS']):<35} {dom:>10}")

elapsed = time.time() - start
print(f"\n{'='*60}")
print(f"Pipeline concluído em {elapsed/60:.1f} minutos.")
print("="*60)
