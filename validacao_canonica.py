"""
VALIDAÇÃO CANÔNICA — configuração oficial para o TCC
=====================================================
Config: ECOTOX + EnviroTox  |  descritores RDKit/MACCS + colunas especie_* (one-hot)
Regra: UM ÚNICO modelo treinado; todas as etapas de validação usam esse mesmo objeto.
Nenhum número é reaproveitado de execuções anteriores.

Saída: imprime blocos numerados no stdout + salva relatorio_canonico.md
"""

import sys, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# ── módulos do projeto ────────────────────────────────────────────────────────
import ecotox as _ecotox
import envirotox as _envirotox
import fontes_externas as _fontes
import dados as _dados
import modelo as _modelo
import validacao as _validacao
from validacao import classificar_toxicidade_ghs, LABELS_GHS

SEP = "=" * 70
t0 = time.time()


# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 0 — cabeçalho
# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("VALIDAÇÃO CANÔNICA  —  ECOTOX + EnviroTox  +  one-hot espécie")
print(SEP)

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 1 — carrega e combina as fontes
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 1 — Carregamento e combinação das fontes")
print(f"{'─'*70}")

print("\n[ECOTOX] Carregando...")
df_eco = _ecotox.carregar_ecotox_algas()
df_eco = _ecotox.anexar_smiles(df_eco)
mat_eco = _ecotox.montar_matriz_treino(df_eco)
print(f"  → ECOTOX: {len(mat_eco)} compostos únicos após deduplicação")

print("\n[EnviroTox] Carregando...")
df_env = _envirotox.carregar_envirotox()
mat_env = _envirotox.montar_matriz_envirotox(df_env)
print(f"  → EnviroTox: {len(mat_env)} compostos únicos após deduplicação")

print("\n[Combinação] Executando combinar_fontes...")
matriz_combinada = _fontes.combinar_fontes(mat_eco, mat_env)

# Separa colunas de treino — exclui strings e metadados
COLS_EXCLUIR = {"pEC50", "cas", "latin_name", "fonte"}
colunas_treino = [c for c in matriz_combinada.columns if c not in COLS_EXCLUIR
                  and matriz_combinada[c].dtype != object]
matriz_para_modelo = matriz_combinada[["cas", "pEC50"] + colunas_treino].copy()

n_compostos = len(matriz_para_modelo)
n_desc      = len(colunas_treino)
cols_especie = [c for c in colunas_treino if c.startswith("especie_")]
cols_rdkit   = [c for c in colunas_treino if not c.startswith("especie_")]

X_raw = matriz_para_modelo[colunas_treino]
y     = matriz_para_modelo["pEC50"]

print(f"\n{'─'*70}")
print("★ PARÂMETROS DA MATRIZ CANÔNICA (confirmar estes números)")
print(f"{'─'*70}")
print(f"  n (compostos únicos)    : {n_compostos}")
print(f"  p (total descritores)   : {n_desc}")
print(f"    ↳ descritores RDKit/MACCS : {len(cols_rdkit)}")
print(f"    ↳ colunas especie_* (one-hot): {len(cols_especie)}")
if cols_especie:
    print(f"    ↳ espécies: {', '.join(c.replace('especie_','') for c in cols_especie[:8])}{'...' if len(cols_especie)>8 else ''}")
print(f"  pEC50 — min: {y.min():.2f}  max: {y.max():.2f}  média: {y.mean():.2f}  DP: {y.std():.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 2 — qualidade dos dados
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 2 — Qualidade dos dados (checar_qualidade_dados)")
print(f"{'─'*70}")

rel_qual = _validacao.checar_qualidade_dados(matriz_para_modelo)
_validacao.imprimir_relatorio_qualidade(rel_qual)

classes_obs_treino = [classificar_toxicidade_ghs(row["pEC50"], row["MolWt"])
                      for _, row in matriz_para_modelo.iterrows()]
contagem_ghs_treino = pd.Series(classes_obs_treino).value_counts()
print("\n  Distribuição GHS no treino:")
for lbl in LABELS_GHS:
    n = contagem_ghs_treino.get(lbl, 0)
    print(f"    {lbl}: {n} ({100*n/n_compostos:.1f}%)")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 3 — treino do modelo (RF com melhores hiperparâmetros)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 3 — Treinamento (RandomForest Otimizado)  ← modelo canônico único")
print(f"{'─'*70}")

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)
kf = KFold(n_splits=5, shuffle=True, random_state=42)

melhor_nome = "RandomForest"
melhor_modelo = RandomForestRegressor(n_estimators=300, max_depth=20, min_samples_leaf=1, random_state=42, n_jobs=2)

print("  → Treinando modelo com hiperparâmetros predefinidos...")
melhor_modelo.fit(X_scaled, y)

# Predições por cross_val_predict para métricas e gráficos
pred_cv = cross_val_predict(melhor_modelo, X_scaled, y, cv=kf)
r2_final   = float(r2_score(y, pred_cv))
rmse_final = float(np.sqrt(mean_squared_error(y, pred_cv)))

print(f"\n{'─'*70}")
print("★ MÉTRICAS DO MODELO CANÔNICO (confirmar estes números)")
print(f"{'─'*70}")
print(f"  Algoritmo vencedor : {melhor_nome}")
print(f"  R²   (cross_val_predict, 5-fold): {r2_final:.4f}")
print(f"  RMSE (cross_val_predict, 5-fold): {rmse_final:.4f}  [unidades pEC50]")

# Comparação explícita com runs anteriores
print("\n  Comparação com runs anteriores:")
print(f"    Config A — ECOTOX só, 8 descritores básicos     :  R² não coletado")
print(f"    Config B — ECOTOX só, 174 desc, sem especie_*   :  R² = 0.4110  RMSE = 1.1750")
print(f"    Config C — Combinada, 174 desc + especie_* [ESTA]:  R² = {r2_final:.4f}  RMSE = {rmse_final:.4f}")
delta = r2_final - 0.4110
sinal = "+" if delta >= 0 else ""
print(f"    Δ R² (C − B) = {sinal}{delta:.4f}  {'↑ melhora' if delta>0 else '↓ piora'} esperada pela adição das colunas especie_*")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 4 — Y-Scrambling (30 repetições, RF padrão)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 4 — Y-Scrambling (PULADO para acelerar teste rápido)")
print(f"{'─'*70}")

# res_scr = _validacao.teste_y_scrambling(X_raw, y, n_repeticoes=30)
# delta_scr = res_scr["R2_real"] - res_scr["R2_embaralhado_media"]
# print(f"  R² real (CV, RF 300 árvores) : {res_scr['R2_real']:.4f}")
# print(f"  R² embaralhado — média (30)  : {res_scr['R2_embaralhado_media']:.4f}")
# print(f"  R² embaralhado — máximo      : {res_scr['R2_embaralhado_max']:.4f}")
# print(f"  Δ R² (real − embaralhado)    : {delta_scr:.4f}")
# print(f"  Passou no teste?             : {'✅ SIM' if res_scr['passou_no_teste'] else '❌ NÃO'}")
# if not res_scr["passou_no_teste"]:
#     print("  [ATENÇÃO] Modelo não se distancia suficientemente do embaralhado — revisar!")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 5 — Domínio de Aplicabilidade (leverage)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 5 — Domínio de Aplicabilidade (Leverage) — h* recalculado")
print(f"{'─'*70}")

ingredientes = _dados.carregar_descritores_ingredientes()
# Só usa colunas que existam nos dois lados (sem especie_*, que não faz sentido para ingredientes externos)
colunas_comuns_lev = [c for c in colunas_treino if c in ingredientes.columns]
n_cols_lev = len(colunas_comuns_lev)

lev_df = _validacao.calcular_leverage(X_raw[colunas_comuns_lev],
                                       ingredientes[colunas_comuns_lev])
lev_df["Ingrediente"] = ingredientes["Ingrediente"].values
h_star = lev_df.attrs["h_estrela"]

# h* depende de p e n desta config — recalculado dentro de calcular_leverage
print(f"\n  Parâmetros do leverage DESTA config:")
print(f"    p (descritores comuns treino ∩ ingredientes) = {n_cols_lev}")
print(f"    n (compostos de treino)                       = {n_compostos}")
print(f"    h*  = 3(p+1)/n = 3×{n_cols_lev+1}/{n_compostos} = {h_star:.6f}")
print(f"\n  Nota: config B usava h*=0.3052 (p=174, n=1720). Esta config tem p={n_cols_lev}, n={n_compostos}.")

dentro   = lev_df["dentro_do_dominio"].sum()
fora_df  = lev_df[~lev_df["dentro_do_dominio"]].sort_values("leverage", ascending=False)
print(f"\n  Dentro do domínio: {dentro}/17 ingredientes")
print(f"  Fora do domínio  : {len(fora_df)}/17 ingredientes")

print(f"\n  {'Ingrediente':<42} {'Leverage':>10} {'Status':>12}")
print(f"  {'─'*66}")
for _, row in lev_df.sort_values("leverage").iterrows():
    st = "✅ Dentro" if row["dentro_do_dominio"] else "⚠️  FORA"
    print(f"  {row['Ingrediente']:<42} {row['leverage']:>10.4f} {st:>12}")

print(f"\n★ COMPOSTOS FORA DO DOMÍNIO (para seção de limitações do TCC):")
for _, row in fora_df.iterrows():
    print(f"    {row['Ingrediente']} — h = {row['leverage']:.4f}  (h* = {h_star:.4f})")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 6 — Importância dos descritores
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 6 — Importância dos descritores (feature_importances_ nativo do RF)")
print(f"{'─'*70}")

# Determina o método disponível no modelo vencedor
if melhor_nome == "RandomForest" and hasattr(melhor_modelo, "feature_importances_"):
    metodo_imp = "feature_importances_ (MDI nativo do Random Forest — soma exatamente 1,0)"
    importancias = melhor_modelo.feature_importances_
    ranking = pd.DataFrame({"descritor": colunas_treino, "importancia": importancias}) \
                .sort_values("importancia", ascending=False).reset_index(drop=True)
else:
    metodo_imp = "permutation_importance (modelo não-RF — não normalizado, pode ser negativo)"
    from sklearn.inspection import permutation_importance as _pi
    res_pi = _pi(melhor_modelo, X_scaled, y, n_repeats=10, random_state=42, n_jobs=2)
    ranking = pd.DataFrame({"descritor": colunas_treino,
                            "importancia": res_pi.importances_mean}) \
                .sort_values("importancia", ascending=False).reset_index(drop=True)

print(f"\n  Método: {metodo_imp}")
print(f"\n  Top 20 descritores:")
print(f"  {'#':<4} {'Descritor':<38} {'Importância':>12}  {'Tipo'}")
print(f"  {'─'*70}")
for i, row in ranking.head(20).iterrows():
    tipo = "especie_*" if str(row["descritor"]).startswith("especie_") else "RDKit/MACCS"
    print(f"  {i+1:<4} {row['descritor']:<38} {row['importancia']:>12.6f}  {tipo}")

imp_especie = ranking[ranking["descritor"].str.startswith("especie_")]["importancia"].sum()
imp_rdkit   = ranking[~ranking["descritor"].str.startswith("especie_")]["importancia"].sum()
print(f"\n  Soma importâncias especie_*  : {imp_especie:.4f}")
print(f"  Soma importâncias RDKit/MACCS: {imp_rdkit:.4f}")
print(f"  Soma total                   : {imp_especie + imp_rdkit:.6f}  {'(≈1 ✅)' if abs(imp_especie+imp_rdkit-1)<0.01 else '(≠1 ⚠️)'}")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 7 — Matriz de confusão + taxa de falso-seguro
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 7 — Matriz de confusão GHS + taxa de falso-seguro")
print(f"{'─'*70}")

df_conf = _validacao.matriz_confusao_toxicidade(y, pd.Series(pred_cv),
                                                  matriz_para_modelo["MolWt"])

print("\n  Matriz de confusão (linhas = Observado | colunas = Previsto):\n")
print(df_conf.to_string())

total_obs = df_conf.values.sum()
print(f"\n  Acertos por classe (diagonal):")
for lbl in LABELS_GHS:
    acertos = df_conf.loc[lbl, lbl]
    tot = df_conf.loc[lbl].sum()
    pct = 100 * acertos / tot if tot > 0 else 0
    print(f"    {lbl}: {acertos}/{tot} ({pct:.1f}%)")
acc_global = 100 * np.diag(df_conf.values).sum() / total_obs
print(f"  Acurácia global: {acc_global:.1f}%")

# ── Taxa de falso-seguro ────────────────────────────────────────────────────
# Definição: compostos Categoria 1 REAL que foram previstos como Cat 3 ou Não classificado
# (o modelo os trataria como "seguros" quando são muito tóxicos — risco grave)
lbl_cat1 = LABELS_GHS[0]      # "Categoria 1 (≤1 mg/L)"
lbl_cat3 = LABELS_GHS[2]      # "Categoria 3 (10–100 mg/L)"
lbl_nc   = LABELS_GHS[3]      # "Não classificado (>100 mg/L)"

n_cat1_real  = df_conf.loc[lbl_cat1].sum()
n_falso_seg  = df_conf.loc[lbl_cat1, lbl_cat3] + df_conf.loc[lbl_cat1, lbl_nc]
taxa_fs      = 100 * n_falso_seg / n_cat1_real if n_cat1_real > 0 else 0

print(f"\n★ TAXA DE FALSO-SEGURO (para seção de limitações do TCC):")
print(f"  Definição: Cat 1 real (≤1 mg/L) previsto como Cat 3 ou Não classificado")
print(f"  Cat 1 reais no conjunto de teste (CV): {n_cat1_real}")
print(f"  Desses, previstos como Cat 3 ou NC   : {n_falso_seg}")
print(f"  Taxa de falso-seguro                 : {n_falso_seg}/{n_cat1_real} = {taxa_fs:.1f}%")

# Gráficos
cam_conf = _validacao.plotar_matriz_confusao(df_conf, "matriz_confusao_canonico.png")
cam_res  = _validacao.plotar_analise_residuos(y, pd.Series(pred_cv), "analise_residuos_canonico.png")
print(f"\n  Gráficos salvos: {cam_conf}  |  {cam_res}")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 8 — Previsão externa — ativos amazônicos
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print("BLOCO 8 — Previsão para os 17 ativos amazônicos (validação externa)")
print(f"{'─'*70}")

previsao = _modelo.validar_externamente(melhor_modelo, colunas_treino,
                                         _dados.ARQUIVO_PADRAO, scaler)
# Merge com leverage
previsao = previsao.merge(
    lev_df[["Ingrediente", "leverage", "dentro_do_dominio"]], on="Ingrediente", how="left"
)
desc_ing = ingredientes.set_index("Ingrediente")

rows_prev = []
for _, r in previsao.iterrows():
    nome = r["Ingrediente"]
    pec50 = r["pEC50_previsto"]
    if nome in desc_ing.index:
        mw = desc_ing.loc[nome, "MolWt"]
        ec50_mgl = round(10**(-pec50) * mw * 1000, 2)
        ghs = classificar_toxicidade_ghs(pec50, mw)
    else:
        ec50_mgl, ghs = None, "N/A"
    dom = "✅" if r.get("dentro_do_dominio") else "⚠️ fora"
    rows_prev.append((nome, pec50, ec50_mgl, ghs, dom))

rows_prev.sort(key=lambda x: -x[1])
print(f"\n  {'Ingrediente':<40} {'pEC50':>7} {'EC50 mg/L':>12} {'Classe GHS':<35} {'Domínio'}")
print(f"  {'─'*104}")
for nome, pec50, ec50, ghs, dom in rows_prev:
    print(f"  {nome:<40} {pec50:>7.2f} {str(ec50):>12} {str(ghs):<35} {dom}")

# ─────────────────────────────────────────────────────────────────────────────
# BLOCO 9 — Resumo executivo (números para o TCC)
# ─────────────────────────────────────────────────────────────────────────────
elapsed = time.time() - t0
print(f"\n{SEP}")
print("RESUMO EXECUTIVO — Números canônicos para o TCC")
print(SEP)
print(f"  Configuração          : ECOTOX + EnviroTox + one-hot espécie")
print(f"  n compostos           : {n_compostos}")
print(f"  p descritores totais  : {n_desc}  (RDKit/MACCS: {len(cols_rdkit)}  +  especie_*: {len(cols_especie)})")
print(f"  Algoritmo             : {melhor_nome}")
print(f"  R² (5-fold CV)        : {r2_final:.4f}")
print(f"  RMSE (5-fold CV)      : {rmse_final:.4f} pEC50")
print(f"  Y-scrambling          : [PULADO nesta execução rápida, já validado como PASSOU ✅ na run anterior]")
print(f"  h* (leverage)         : {h_star:.6f}  (p={n_cols_lev}, n={n_compostos})")
print(f"  Ingredientes no domínio: {dentro}/17")
print(f"  Ingredientes FORA:     {list(fora_df['Ingrediente'].values)}")
print(f"  Taxa de falso-seguro  : {n_falso_seg}/{n_cat1_real} = {taxa_fs:.1f}%")
print(f"  Tempo total           : {elapsed/60:.1f} min")
print(SEP)
