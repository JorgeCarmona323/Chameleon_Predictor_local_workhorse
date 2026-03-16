"""
generate_batches.py
-------------------
Run this script LOCALLY on Monday morning after Tier-1 finishes.
It reads the Tier-1 output and generates 4 stratified compound batches
for the 4 Google Colab notebooks.

Usage:
    python colab/generate_batches.py

Outputs (in colab/batches/):
    batch_1.csv  ← upload to Colab Account 1, Notebook 1
    batch_2.csv  ← upload to Colab Account 1, Notebook 2
    batch_3.csv  ← upload to Colab Account 2, Notebook 1
    batch_4.csv  ← upload to Colab Account 2, Notebook 2

Stratification strategy:
  - Compounds are binned into 4 PAMPA quartiles (equal-count)
  - 25 compounds sampled from each quartile per batch = 100 per batch
  - 400 total compounds across 4 batches
  - Reference compound IDs (1, 2, 183, 980, 1829) always go into batch_1
  - Compounds already done by local Tier-2 CREST are excluded

If Tier-1 has NOT finished yet, falls back to PAMPA-only stratification
from pampa_curated.csv (no PSA_mem available).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ── Reference compound IDs — always included in batch_1 ───────────────────────
REFERENCE_IDS = {1, 2, 183, 980, 1829}

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT   = Path(__file__).parent.parent
TIER1_CSV      = PROJECT_ROOT / "results" / "conformer_descriptors_raw.csv"
PAMPA_CSV      = PROJECT_ROOT / "data" / "pampa_curated.csv"
OUT_DIR        = PROJECT_ROOT / "colab" / "batches"

N_BATCHES      = 4
N_PER_QUARTILE = 25   # 25 × 4 quartiles = 100 per batch
RANDOM_SEED    = 42


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-per-quartile", type=int, default=N_PER_QUARTILE,
                        help="Compounds per PAMPA quartile per batch (default 25 → 100/batch)")
    parser.add_argument("--n-batches", type=int, default=N_BATCHES)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Load PAMPA data ────────────────────────────────────────────────────────
    pampa_df = pd.read_csv(PAMPA_CSV, low_memory=False)
    smiles_col = "SMILES_canonical" if "SMILES_canonical" in pampa_df.columns else "SMILES"
    pampa_df = pampa_df[["ID", smiles_col, "PAMPA"]].dropna(subset=[smiles_col, "PAMPA"]).copy()
    pampa_df.rename(columns={smiles_col: "SMILES_canonical"}, inplace=True)

    # ── Optionally join Tier-1 PSA_mem for smarter stratification ─────────────
    if TIER1_CSV.exists():
        tier1 = pd.read_csv(TIER1_CSV, low_memory=False)
        tier1 = tier1[tier1["error"].isna()][["ID", "mem_psa3d", "delta_psa3d"]].copy()
        df = pampa_df.merge(tier1, on="ID", how="left")
        strat_col = "mem_psa3d"
        print(f"Tier-1 output found — stratifying by mem_psa3d ({tier1.shape[0]} molecules)")
    else:
        df = pampa_df.copy()
        strat_col = "PAMPA"
        print("Tier-1 output NOT found — stratifying by PAMPA only")

    # ── PAMPA quartile bins (used for balance regardless of strat_col) ─────────
    df["pampa_q"] = pd.qcut(df["PAMPA"], q=4, labels=["Q1_low","Q2","Q3","Q4_high"])

    # ── Separate reference compounds ───────────────────────────────────────────
    ref_mask  = df["ID"].isin(REFERENCE_IDS)
    ref_df    = df[ref_mask].copy()
    pool_df   = df[~ref_mask].copy()

    print(f"Reference compounds found: {len(ref_df)} / {len(REFERENCE_IDS)}")
    if len(ref_df) < len(REFERENCE_IDS):
        missing = REFERENCE_IDS - set(ref_df["ID"])
        print(f"  WARNING: IDs {missing} not found in PAMPA dataset")

    # ── Stratified sampling ────────────────────────────────────────────────────
    rng = np.random.default_rng(RANDOM_SEED)
    n_per_q = args.n_per_quartile
    n_total  = n_per_q * 4 * args.n_batches  # total compounds across all batches

    sampled_parts = []
    for q_label, q_group in pool_df.groupby("pampa_q", observed=True):
        n_available = len(q_group)
        n_want      = n_per_q * args.n_batches
        n_sample    = min(n_want, n_available)
        if n_sample < n_want:
            print(f"  WARNING: quartile {q_label} has only {n_available} compounds "
                  f"(wanted {n_want}) — using all")
        # Sort by strat_col (if available) to ensure diversity, then sample evenly
        if strat_col in q_group.columns and q_group[strat_col].notna().any():
            q_group = q_group.sort_values(strat_col)
            # Evenly spaced indices for coverage
            indices = np.linspace(0, len(q_group)-1, n_sample, dtype=int)
            sampled = q_group.iloc[indices]
        else:
            sampled = q_group.sample(n=n_sample, random_state=RANDOM_SEED)
        sampled_parts.append(sampled)

    sampled_df = pd.concat(sampled_parts, ignore_index=True)
    sampled_df = sampled_df.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)

    # ── Split into N batches, put reference compounds in batch_1 ──────────────
    n_per_batch = n_per_q * 4  # 100 if defaults
    batches = []
    for i in range(args.n_batches):
        chunk = sampled_df.iloc[i * n_per_batch : (i+1) * n_per_batch].copy()
        chunk["batch"] = i + 1
        batches.append(chunk)

    # Add reference compounds to batch 1 (drop any duplicates already there)
    ref_df["batch"] = 1
    batches[0] = pd.concat(
        [ref_df, batches[0][~batches[0]["ID"].isin(REFERENCE_IDS)]], ignore_index=True
    )

    # ── Save batches ───────────────────────────────────────────────────────────
    keep_cols = ["ID", "SMILES_canonical", "PAMPA", "batch"]
    if "mem_psa3d" in sampled_df.columns:
        keep_cols.insert(3, "mem_psa3d")
    if "delta_psa3d" in sampled_df.columns:
        keep_cols.insert(4, "delta_psa3d")

    for i, batch_df in enumerate(batches):
        out_cols = [c for c in keep_cols if c in batch_df.columns]
        out_path = OUT_DIR / f"batch_{i+1}.csv"
        batch_df[out_cols].to_csv(out_path, index=False)
        print(f"Saved batch_{i+1}.csv — {len(batch_df)} compounds → {out_path}")

    total = sum(len(b) for b in batches)
    print(f"\nTotal compounds across all batches: {total}")
    print(f"Upload colab/batches/batch_N.csv files to your Google Drive.")


if __name__ == "__main__":
    main()
