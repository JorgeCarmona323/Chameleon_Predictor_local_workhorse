# env: chameleon-calc
"""
00_curate_data.py
-----------------
Load CycPeptMPDB CSV, validate SMILES, build PAMPA subset,
and curate the chameleonic reference set (CycloA + analogs).

Outputs (written to ../data/):
  pampa_curated.csv       — PAMPA subset with valid SMILES
  reference_set.csv       — CycloA + 4 analogs with known permeability

Usage:
  python curate_data.py [--input <path_to_csv>] [--outdir <data_dir>]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

# Reference compounds: CycloA (ID=1) + structurally verified analogs from literature
# IDs from CycPeptMPDB; original name given for traceability
REFERENCE_IDS = {
    1:    "Cyclosporine_A",         # CycloA, canonical chameleonic reference
    22:   "CycloA_same_1",
    932:  "CycloA_same_2",
    981:  "CycloA_same_3",
    1822: "CycloA_same_4",
}


def standardize_smiles(smi: str) -> str | None:
    """Parse, standardize, and return canonical isomeric SMILES. Returns None on failure."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        mol = rdMolStandardize.Normalizer().normalize(mol)
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def curate(input_csv: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {input_csv} ...")
    df = pd.read_csv(input_csv, low_memory=False)
    print(f"  Total rows: {len(df)}")

    # ── Canonicalize SMILES ──────────────────────────────────────────────────
    print("Canonicalizing SMILES ...")
    df["SMILES_canonical"] = df["SMILES"].apply(
        lambda s: standardize_smiles(str(s)) if pd.notna(s) else None
    )
    n_failed = df["SMILES_canonical"].isna().sum()
    print(f"  SMILES parse failures: {n_failed}")
    df = df[df["SMILES_canonical"].notna()].copy()
    print(f"  Rows after SMILES filter: {len(df)}")

    # ── PAMPA subset ─────────────────────────────────────────────────────────
    pampa = df[df["PAMPA"].notna()].copy()
    print(f"\nPAMPA subset: {len(pampa)} compounds")

    # Drop rows missing both 3D PSA columns (needed for ΔΨ baseline)
    pampa_with_3dpsa = pampa[pampa["CHCl3_3DPSA"].notna() & pampa["H2O_3DPSA"].notna()].copy()
    print(f"PAMPA with CHCl3/H2O 3DPSA: {len(pampa_with_3dpsa)} compounds")

    # Compute ΔΨ baseline from DB values (H2O - CHCl3 = aqueous exposure delta)
    pampa_with_3dpsa["delta_3DPSA_db"] = (
        pampa_with_3dpsa["H2O_3DPSA"] - pampa_with_3dpsa["CHCl3_3DPSA"]
    )

    # Binary permeability label — added BEFORE writing CSV so downstream scripts get it
    threshold = -6.0
    pampa["permeable"] = (pampa["PAMPA"] >= threshold).astype(int)
    n_perm = pampa["permeable"].sum()
    print(f"\nPermeable (PAMPA >= {threshold}): {n_perm} / {len(pampa)} ({100*n_perm/len(pampa):.1f}%)")

    # Save threshold to sidecar JSON so downstream scripts always use the same value
    import json
    (outdir / "config.json").write_text(json.dumps({"pampa_threshold": threshold}, indent=2))

    # Save PAMPA curated (all with SMILES, keep rows missing 3DPSA for Tier-1 conformer run)
    pampa.to_csv(outdir / "pampa_curated.csv", index=False)
    pampa_with_3dpsa.to_csv(outdir / "pampa_with_db_3dpsa.csv", index=False)
    print(f"Saved: {outdir / 'pampa_curated.csv'}  ({len(pampa)} rows)")
    print(f"Saved: {outdir / 'pampa_with_db_3dpsa.csv'}  ({len(pampa_with_3dpsa)} rows)")

    # ── Reference set ────────────────────────────────────────────────────────
    ref_ids = list(REFERENCE_IDS.keys())
    ref = df[df["ID"].isin(ref_ids)].copy()
    ref["Ref_Name"] = ref["ID"].map(REFERENCE_IDS)

    # Add any Structurally_Unique_ID==1 (CycloA group) as additional analogs
    if "Structurally_Unique_ID" not in df.columns:
        raise KeyError(
            "Column 'Structurally_Unique_ID' not found in CSV. "
            "Check CycPeptMPDB version — expected v1.2."
        )
    cycloA_group = df[df["Structurally_Unique_ID"] == 1].copy()
    ref = pd.concat([ref, cycloA_group], ignore_index=True).drop_duplicates(subset="ID")
    print(f"\nReference set: {len(ref)} compounds (CycloA + structurally same)")

    ref.to_csv(outdir / "reference_set.csv", index=False)
    print(f"Saved: {outdir / 'reference_set.csv'}")

    # ── Summary stats ────────────────────────────────────────────────────────
    print("\n── PAMPA permeability distribution ──")
    print(pampa["PAMPA"].describe().round(3).to_string())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Curate CycPeptMPDB for PAMPA analysis")
    parser.add_argument(
        "--input", "-i",
        default=str(Path(__file__).parents[1] / "CycPeptMPDB_Peptide_All (2).csv"),
        help="Path to CycPeptMPDB CSV",
    )
    parser.add_argument(
        "--outdir", "-o",
        default=str(Path(__file__).parents[1] / "data"),
        help="Output directory for curated CSVs",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    curate(args.input, Path(args.outdir))
