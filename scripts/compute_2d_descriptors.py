"""
compute_2d_descriptors.py
-------------------------
Canonical 2D / lipophilicity descriptors for the reference compounds, to confirm
which descriptors are blind to stereochemistry (i.e. identical for R/S isomers).

Reads data/new_6mer_compounds.csv (name,smiles), writes results/2d_descriptors.csv.
"""
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
RDLogger.DisableLog("rdApp.*")

SRC = Path("data/new_6mer_compounds.csv")
OUT = Path("results/2d_descriptors.csv")

DESCRIPTORS = {
    "MolWt":            Descriptors.MolWt,
    "ExactMolWt":       Descriptors.ExactMolWt,
    "TPSA_2d":          Descriptors.TPSA,
    "MolLogP_Crippen":  Crippen.MolLogP,
    "MolMR_Crippen":    Crippen.MolMR,          # molar refractivity (lipophilicity-related)
    "NumHDonors":       Descriptors.NumHDonors,
    "NumHAcceptors":    Descriptors.NumHAcceptors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "FractionCSP3":     Descriptors.FractionCSP3,
    "NumAromaticRings": Descriptors.NumAromaticRings,
    "NumHeavyAtoms":    Descriptors.HeavyAtomCount,
    "LabuteASA":        Descriptors.LabuteASA,   # approx surface area
    "QED":              QED.qed,
}


def main():
    df = pd.read_csv(SRC)
    rows = []
    for _, r in df.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            continue
        row = {"name": r["name"]}
        for k, fn in DESCRIPTORS.items():
            try:
                row[k] = round(float(fn(mol)), 4)
            except Exception:
                row[k] = None
        rows.append(row)
    out = pd.DataFrame(rows).set_index("name")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT)

    # focus: the R/S isomer pair
    iso = [n for n in out.index if "3-12-8-12" in n]
    print("=== All reference 6-mers — 2D descriptors ===")
    print(out.T.to_string())
    if len(iso) == 2:
        a, b = iso
        print(f"\n=== R vs S isomer comparison ({a} vs {b}) ===")
        diff = (out.loc[a] != out.loc[b])
        identical = out.columns[~diff.values] if hasattr(diff, "values") else []
        differing = [c for c in out.columns if out.loc[a, c] != out.loc[b, c]]
        print(f"IDENTICAL 2D descriptors: {len(out.columns) - len(differing)}/{len(out.columns)}")
        if differing:
            print(f"DIFFERING: {differing}")
            print(out.loc[[a, b], differing].T.to_string())
        else:
            print("ALL 2D / lipophilicity descriptors are IDENTICAL between R and S.")
    print(f"\nSaved -> {OUT}")


if __name__ == "__main__":
    main()
