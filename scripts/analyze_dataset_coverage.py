"""
analyze_dataset_coverage.py
---------------------------
Inventory of data availability for the permeability ML model, by residue count:
  - CycPeptMPDB permeability labels (feature_matrix.csv)
  - CREMP CHCl3 conformer ensembles (cremp_deltapsa.csv)
  - vacuum RDKit dynamic descriptors (feature_matrix.csv aq_psa3d etc.)

Motivation: chameleonicity matters most above ~9 residues, but CREMP ensembles
cover only 6-7mers. This quantifies the label-vs-ensemble coverage gap.

Output: results/dataset_coverage.csv
See docs/experiments/2026-06-09_dataset_coverage_gap.md
"""
from pathlib import Path
import numpy as np
import pandas as pd


def nres_from_id(cid):
    try:
        return len(str(cid).split("."))
    except Exception:
        return np.nan


def main():
    fm = pd.read_csv("results/archive/feature_matrix.csv", low_memory=False)
    fm["Monomer_Length"] = pd.to_numeric(fm["Monomer_Length"], errors="coerce")
    fm["PAMPA"] = pd.to_numeric(fm["PAMPA"], errors="coerce")
    has_vac = "aq_psa3d" in fm.columns

    cremp = pd.read_csv("results/archive/cremp_deltapsa.csv", low_memory=False)
    cremp["n_res"] = cremp["compound_id"].apply(nres_from_id)
    cremp_ok = cremp[cremp["error"].isna()] if "error" in cremp.columns else cremp
    cremp_counts = cremp_ok["n_res"].value_counts().to_dict()

    rows = []
    for n in range(4, 16):
        sub = fm[fm["Monomer_Length"] == n]
        if len(sub) == 0 and n not in cremp_counts:
            continue
        rows.append({
            "n_res": n,
            "cycpeptmpdb_total": len(sub),
            "cycpeptmpdb_PAMPA": int(sub["PAMPA"].notna().sum()),
            "vacuum_descriptors": int(sub["aq_psa3d"].notna().sum()) if has_vac else 0,
            "cremp_chcl3_ensembles": int(cremp_counts.get(n, 0)),
        })
    cov = pd.DataFrame(rows)

    # regime rollups
    big9 = fm[fm["Monomer_Length"] >= 9]
    big11 = fm[fm["Monomer_Length"] >= 11]

    out = Path("results/dataset_coverage.csv")
    cov.to_csv(out, index=False)

    print(cov.to_string(index=False))
    print()
    print(f">=9-mer  (chameleonic regime): {len(big9):5d} labeled | "
          f"{cov[cov.n_res>=9]['cremp_chcl3_ensembles'].sum():4d} CREMP ensembles")
    print(f">=11-mer (CsA regime):         {len(big11):5d} labeled | "
          f"{cov[cov.n_res>=11]['cremp_chcl3_ensembles'].sum():4d} CREMP ensembles")
    print(f"\nSaved -> {out}")


if __name__ == "__main__":
    main()
