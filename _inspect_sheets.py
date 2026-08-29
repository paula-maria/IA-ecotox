import pandas as pd, sys

for arq in ['algae.xlsx', 'chlorella.xlsx']:
    try:
        xl = pd.ExcelFile(arq)
        print(f'\n=== {arq} === abas: {xl.sheet_names}')
        for aba in xl.sheet_names:
            df = xl.parse(aba, nrows=2)
            print(f'  [{aba}] colunas: {list(df.columns)[:12]}')
    except Exception as e:
        print(f'ERRO {arq}: {e}')
sys.stdout.flush()
