import pandas as pd
fm = pd.read_csv('results/feature_matrix.csv', low_memory=False)
for name in ['DP-944', 'DP-955', 'c*[PSLYF]']:
    row = fm[fm['Original_Name_in_Source_Literature'].str.contains(name, na=False, case=False)]
    if len(row):
        r = row.iloc[0]
        print(f"{name}")
        print(f"  ID: {r['ID']}")
        print(f"  SMILES: {r['SMILES']}")
        print(f"  PAMPA: {r['Permeability']}")
        print(f"  Monomer_Length: {r['Monomer_Length']}")
        print(f"  Source: {r['Source']}")
        print()
