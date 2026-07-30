"""
Script para exportar as matrizes pré-computadas (ECOTOX e COMBINADA) para uso em deploy.
Isso evita que o servidor na nuvem (que não tem acesso a ecotox_ascii/) precise rodar o RDKit.
"""
import ecotox
import envirotox
import fontes_externas
import pandas as pd

def exportar():
    print("[INFO] Gerando matriz ECOTOX...")
    dados_ecotox = ecotox.carregar_ecotox_algas(apenas_chlorella=False)
    dados_ecotox = ecotox.anexar_smiles(dados_ecotox)
    matriz_ecotox = ecotox.montar_matriz_treino(dados_ecotox)
    
    # Exporta ECOTOX puro
    matriz_ecotox.to_csv("matriz_ecotox_deploy.csv", index=False)
    print(f"[OK] matriz_ecotox_deploy.csv gerada ({len(matriz_ecotox)} linhas).")
    
    print("[INFO] Gerando matriz EnviroTox e combinando...")
    df_et = envirotox.carregar_envirotox()
    matriz_et = envirotox.montar_matriz_envirotox(df_et)
    matriz_combinada = fontes_externas.combinar_fontes(matriz_ecotox, matriz_et)
    
    # Adiciona colunas extras necessárias no app.py para contagens no dashboard do modo combinado
    matriz_combinada.to_csv("matriz_combinada_deploy.csv", index=False)
    
    # Vamos salvar também o n_ecotox e n_envirotox para estatísticas
    stats = pd.DataFrame([{"n_ecotox": len(matriz_ecotox), "n_envirotox_cas_unicos": matriz_et["cas"].nunique()}])
    stats.to_csv("matriz_combinada_stats.csv", index=False)
    
    print(f"[OK] matriz_combinada_deploy.csv gerada ({len(matriz_combinada)} linhas).")
    
    print("[INFO] Matrizes exportadas com sucesso! Agora você pode dar git add nas duas.")

if __name__ == "__main__":
    exportar()
