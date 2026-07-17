#!/usr/bin/env python3
"""stretchiness_by_size.py — is vacuum ΔPSA a within-size-class signal or a size proxy?

Single-feature AUC-ROC of the vacuum "stretchiness" descriptors (accessible ΔPSA range /
flexibility) vs. PAMPA permeability, STRATIFIED by residue count (Monomer_Length: all / <9 /
≥9). Run on the full dataset and, with --sources, on a clean-protocol subset. This is the
cheap first test behind the multi-fidelity plan: if the stretchiness AUC survives *within* a
size bin it is a real capacity signal; if it only shows up pooled, it is a size proxy.

See docs/experiments/2026-07-17_stretchiness_by_size_auc.md for the first run + findings.

Usage:
  python scripts/ml/stretchiness_by_size.py --matrix results/feature_matrix.csv \
      --out results/stretchiness_by_size_full.csv
  # clean-protocol subset:
  python scripts/ml/stretchiness_by_size.py --matrix results/feature_matrix.csv \
      --sources 2016_Furukawa 2013_CHUGAI --out results/stretchiness_by_size_clean.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# stretchiness (capacity) features + a few references for context
DEFAULT_FEATURES = [
    "delta_psa3d", "psa3d_std", "psa3d_spread", "delta_psa3d_per_mw",
    "delta_hb", "delta_Rg", "MolLogP", "TPSA", "MolWt",
]
THRESHOLD = -6.0        # PAMPA LogPexp >= threshold => permeable


def single_feature_auc(df: pd.DataFrame, feat: str, label: str) -> dict:
    sub = df[[feat, label]].apply(pd.to_numeric, errors="coerce").dropna()
    y = sub[label].astype(int).values
    n, npos = len(sub), int(sub[label].sum())
    if npos == 0 or npos == n:          # need both classes
        return {"n": n, "permrate": round(npos / n, 3) if n else np.nan, "AUC": np.nan}
    auc = roc_auc_score(y, sub[feat].values)
    return {"n": n, "permrate": round(npos / n, 3), "AUC": round(float(auc), 3)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", default="results/feature_matrix.csv")
    ap.add_argument("--out", default="results/stretchiness_by_size.csv", type=Path)
    ap.add_argument("--features", nargs="+", default=DEFAULT_FEATURES)
    ap.add_argument("--sources", nargs="*", default=None,
                    help="restrict to these Source values (e.g. 2016_Furukawa 2013_CHUGAI)")
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    args = ap.parse_args()

    df = pd.read_csv(args.matrix, low_memory=False)
    if args.sources:
        df = df[df["Source"].isin(args.sources)].copy()
        print(f"Source filter {args.sources}: {len(df)} compounds")

    # permeability label: prefer existing column, else derive from PAMPA
    if "permeable" not in df.columns:
        df["permeable"] = (pd.to_numeric(df["PAMPA"], errors="coerce") >= args.threshold).astype(int)

    ml = pd.to_numeric(df["Monomer_Length"], errors="coerce")
    bins = {
        "all": df,
        "<9":  df[ml < 9],
        ">=9": df[ml >= 9],
    }
    print("bin sizes: " + "  ".join(f"{k}={len(v)}" for k, v in bins.items()))

    rows = []
    for feat in args.features:
        if feat not in df.columns:
            print(f"  (skip {feat}: not in matrix)")
            continue
        for bn, bdf in bins.items():
            r = single_feature_auc(bdf, feat, "permeable")
            rows.append({"feature": feat, "bin": bn, **r})

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    with pd.option_context("display.max_rows", None, "display.width", 120):
        print(out.to_string(index=False))
    print(f"\nsaved -> {args.out}")
    print("Read: a feature whose pooled AUC (bin=all) drops to ~0.5 within <9 AND >=9 is a "
          "size proxy, not a within-class signal.")


if __name__ == "__main__":
    main()
