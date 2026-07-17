#!/usr/bin/env python3
"""auc_by_subset.py — the "Figure 0" grouped AUC bar chart, done reproducibly.

Single-feature AUC-ROC of key descriptors on the CLEAN subset (Furukawa + Chugai) vs. the
FULL dataset, as a grouped bar chart. Replaces the notebook-buried version whose y-axis was
truncated at ~0.43 — which CLIPPED any bar below that (e.g. MolLogP on the clean subset =
0.317, NumHDonors on full = 0.406) and left their value labels floating with no visible bar.
Here the y-axis starts at 0, so every bar renders with its label attached.

Usage:
  python scripts/ml/auc_by_subset.py --matrix results/feature_matrix.csv \
      --out results/figures/auc_by_subset.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

CLEAN_SOURCES = ["2016_Furukawa", "2013_CHUGAI"]
THRESHOLD = -6.0

# (column, display label, group) — group drives color/shading
DESCRIPTORS = [
    ("delta_psa3d_per_mw", "ΔPSA / MolWt\n(per-MW normalized)", "3d"),
    ("delta_psa3d",        "Absolute ΔPSA\n(Tier-1 ensemble)",  "3d"),
    ("MolLogP",            "MolLogP\n(2D baseline)",            "2d"),
    ("NumHDonors",         "NumHDonors\n(2D baseline)",         "2d"),
    ("TPSA",               "TPSA\n(2D baseline)",               "2d"),
    ("delta_3DPSA_db",     "DB ΔPSA\n(single-structure)",       "db"),
]


def auc(df: pd.DataFrame, feat: str, label: str) -> float:
    sub = df[[feat, label]].apply(pd.to_numeric, errors="coerce").dropna()
    y = sub[label].astype(int).values
    if y.min() == y.max():
        return np.nan
    return float(roc_auc_score(y, sub[feat].values))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", default="results/feature_matrix.csv")
    ap.add_argument("--out", default="results/figures/auc_by_subset.png", type=Path)
    ap.add_argument("--sources", nargs="*", default=CLEAN_SOURCES)
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    args = ap.parse_args()

    df = pd.read_csv(args.matrix, low_memory=False)
    if "permeable" not in df.columns:
        df["permeable"] = (pd.to_numeric(df["PAMPA"], errors="coerce") >= args.threshold).astype(int)
    clean = df[df["Source"].isin(args.sources)]
    n_clean, n_full = len(clean), len(df)
    print(f"clean (={'+'.join(args.sources)}) n={n_clean}   full n={n_full}")

    labels, auc_clean, auc_full, groups = [], [], [], []
    for col, disp, grp in DESCRIPTORS:
        if col not in df.columns:
            print(f"  (skip {col}: not in matrix)")
            continue
        labels.append(disp); groups.append(grp)
        auc_clean.append(auc(clean, col, "permeable"))
        auc_full.append(auc(df, col, "permeable"))

    x = np.arange(len(labels)); w = 0.38
    fig, ax = plt.subplots(figsize=(13, 6))
    b1 = ax.bar(x - w/2, auc_clean, w, label=f"Clean subset (Furukawa + Chugai, n={n_clean})",
                color="#1f9bd6", edgecolor="white")
    b2 = ax.bar(x + w/2, auc_full, w, label=f"Full dataset (n={n_full})",
                color="#f56a3c", edgecolor="white")
    for bars in (b1, b2):
        for r in bars:
            h = r.get_height()
            if np.isfinite(h):
                ax.annotate(f"{h:.3f}", (r.get_x() + r.get_width()/2, h),
                            ha="center", va="bottom", fontsize=9, fontweight="bold")

    # shade the 3D ensemble descriptors
    for i, grp in enumerate(groups):
        if grp == "3d":
            ax.axvspan(i - 0.5, i + 0.5, color="#eef4fa", zorder=0)

    ax.axhline(0.5, ls="--", color="grey", lw=1, label="Chance (AUC = 0.5)")
    ax.set_ylim(0, 0.80)                      # <-- floor at 0: no clipped bars / floating labels
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("AUC-ROC")
    ax.set_title("AUC-ROC by descriptor: clean subset vs. full dataset\n"
                 "Per-MW normalization recovers within-size chameleonic signal; "
                 "single-structure baseline stays at chance", fontweight="bold")
    ax.legend(loc="upper right", framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"saved -> {args.out}")


if __name__ == "__main__":
    main()
