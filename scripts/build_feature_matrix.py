"""
02_build_feature_matrix.py
--------------------------
Merge conformer Δ descriptors with:
  - CycPeptMPDB database 3DPSA values (ΔΨ_db = H2O_3DPSA - CHCl3_3DPSA)
  - RDKit 2D descriptors (MolWt, MolLogP, TPSA, HBA, HBD) as baseline

Final feature matrix columns (all compounds in PAMPA subset):
  Group A — DB Δ features (available for ~88% of PAMPA subset):
    delta_3DPSA_db, H2O_3DPSA, CHCl3_3DPSA

  Group B — Tier-1 conformer Δ features (from conformer_engine.py):
    delta_psa3d, delta_hb, delta_Rg, delta_NPR1, delta_NPR2,
    delta_Asphericity, psa3d_spread, psa3d_std, hb_spread

  Group C — 2D baseline descriptors (from RDKit, all compounds):
    MolWt, MolLogP, TPSA, NumHAcceptors, NumHDonors,
    NumRotatableBonds, FractionCSP3, RingCount

  Target: PAMPA (log permeability)
  Binary target: permeable = (PAMPA >= -6.0)

Usage:
  python build_feature_matrix.py [--pampa data/pampa_curated.csv]
                                  [--conformers results/conformer_descriptors_raw.csv]
                                  [--outdir results]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

PAMPA_THRESHOLD = -6.0  # Jiang et al. 2023, J. Chem. Inf. Model. — CycPeptMPDB standard cutoff

# 2D descriptors to compute from SMILES
RDKIT_2D = {
    "MolWt":            Descriptors.MolWt,
    "MolLogP":          Descriptors.MolLogP,
    "TPSA":             Descriptors.TPSA,
    "NumHAcceptors":    rdMolDescriptors.CalcNumHBA,
    "NumHDonors":       rdMolDescriptors.CalcNumHBD,
    "NumRotatableBonds":rdMolDescriptors.CalcNumRotatableBonds,
    "FractionCSP3":     rdMolDescriptors.CalcFractionCSP3,
    "RingCount":        rdMolDescriptors.CalcNumRings,
    "NumHeavyAtoms":    rdMolDescriptors.CalcNumHeavyAtoms,
    "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings,
}


def compute_2d_descriptors(smiles_series: pd.Series) -> pd.DataFrame:
    rows = []
    for smi in smiles_series:
        try:
            mol = Chem.MolFromSmiles(str(smi))
            if mol is None:
                rows.append({k: np.nan for k in RDKIT_2D})
                continue
            rows.append({k: fn(mol) for k, fn in RDKIT_2D.items()})
        except Exception:
            rows.append({k: np.nan for k in RDKIT_2D})
    return pd.DataFrame(rows, index=smiles_series.index)


def build(pampa_csv: str, conformers_csv: str | None, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load PAMPA curated ───────────────────────────────────────────────────
    df = pd.read_csv(pampa_csv, low_memory=False)
    smiles_col = "SMILES_canonical" if "SMILES_canonical" in df.columns else "SMILES"
    print(f"PAMPA subset: {len(df)} compounds")

    # ── 2D descriptors ───────────────────────────────────────────────────────
    # CycPeptMPDB CSV already contains most RDKit 2D descriptors.
    # Only compute the ones that are missing to avoid duplicate columns.
    missing_2d = [k for k in RDKIT_2D if k not in df.columns]
    if missing_2d:
        print(f"Computing {len(missing_2d)} missing 2D descriptors: {missing_2d}")
        desc_2d = compute_2d_descriptors(df[smiles_col])
        desc_2d = desc_2d[missing_2d]  # only add missing columns
        df = df.reset_index(drop=True)
        desc_2d = desc_2d.reset_index(drop=True)
        df = pd.concat([df, desc_2d], axis=1)
    else:
        print(f"All 2D descriptors already present in CSV — skipping recomputation")

    # ── DB Δ features ────────────────────────────────────────────────────────
    if "H2O_3DPSA" in df.columns and "CHCl3_3DPSA" in df.columns:
        df["delta_3DPSA_db"] = df["H2O_3DPSA"] - df["CHCl3_3DPSA"]
        n_db = df["delta_3DPSA_db"].notna().sum()
        print(f"DB delta_3DPSA available for {n_db} / {len(df)} compounds")
    else:
        df["delta_3DPSA_db"] = np.nan
        print("Warning: H2O_3DPSA / CHCl3_3DPSA not found in PAMPA CSV")

    # ── Tier-1 conformer Δ features ──────────────────────────────────────────
    if conformers_csv and Path(conformers_csv).exists():
        print(f"Loading conformer descriptors from {conformers_csv} ...")
        conf_df = pd.read_csv(conformers_csv)
        conf_df = conf_df[conf_df["error"].isna()].copy()
        print(f"  Successful conformer runs: {len(conf_df)}")

        delta_cols = [c for c in conf_df.columns if c.startswith("delta_")
                      or c.startswith("psa3d") or c.startswith("hb_spread")
                      or c.startswith("aq_") or c.startswith("mem_")]
        conf_df = conf_df[["ID"] + delta_cols].copy()

        df = df.merge(conf_df, on="ID", how="left")
        n_conf = df["delta_psa3d"].notna().sum() if "delta_psa3d" in df.columns else 0
        print(f"  Merged Tier-1 Δ features for {n_conf} / {len(df)} compounds")
    else:
        print("No conformer descriptors file found — skipping Tier-1 Δ features")
        print("  Run conformer_engine.py first, then re-run this script")

    # ── Binary permeability label ─────────────────────────────────────────────
    df["permeable"] = (df["PAMPA"] >= PAMPA_THRESHOLD).astype(int)
    print(f"\nPermeable (>= {PAMPA_THRESHOLD}): {df['permeable'].sum()} / {df['permeable'].notna().sum()}")

    # ── Feature sets for analysis ─────────────────────────────────────────────
    feature_groups = {
        "2D_baseline": [
            "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
            "NumRotatableBonds", "FractionCSP3", "RingCount",
        ],
        "DB_delta": ["delta_3DPSA_db", "H2O_3DPSA", "CHCl3_3DPSA"],
        "Tier1_delta": [
            "delta_psa3d", "delta_hb", "delta_Rg", "delta_NPR1", "delta_NPR2",
            "delta_Asphericity", "psa3d_spread", "psa3d_std", "hb_spread",
        ],
        "combined": [
            "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
            "delta_3DPSA_db", "delta_psa3d", "delta_hb", "delta_Rg",
            "psa3d_spread", "delta_NPR1", "delta_NPR2",
        ],
    }

    import json
    (outdir / "feature_groups.json").write_text(json.dumps(feature_groups, indent=2))

    # ── Save full feature matrix ──────────────────────────────────────────────
    out_path = outdir / "feature_matrix.csv"
    df.to_csv(out_path, index=False)
    print(f"\nSaved feature matrix: {out_path}")
    print(f"  Shape: {df.shape}")

    # ── Coverage report ───────────────────────────────────────────────────────
    print("\n── Feature coverage ──")
    all_features = [f for grp in feature_groups.values() for f in grp]
    all_features = list(dict.fromkeys(all_features))  # deduplicate
    for feat in all_features:
        if feat in df.columns:
            n_valid = df[feat].notna().values.sum()
            print(f"  {feat:30s}: {n_valid:5d} / {len(df)}")
        else:
            print(f"  {feat:30s}: NOT COMPUTED")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build feature matrix for PAMPA analysis")
    parser.add_argument("--pampa",      "-p", default="data/pampa_curated.csv")
    parser.add_argument("--conformers", "-c", default="results/conformer_descriptors_raw.csv")
    parser.add_argument("--outdir",     "-o", default="results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build(
        pampa_csv=args.pampa,
        conformers_csv=args.conformers,
        outdir=Path(args.outdir),
    )
