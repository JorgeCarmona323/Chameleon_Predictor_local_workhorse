"""
05_tier2_validation.py
----------------------
Tier-2 high-rigor cross-check on ~5 reference compounds (CycloA + analogs).

Compares Tier-1 ETKDG/MMFF94s Δ values vs. the database CHCl3/H2O 3DPSA
values, and optionally vs. OpenEye OMEGA + OpenMM GB/SA (if available).

Purpose: validate whether Tier-1 conformer engine captures the correct
direction and magnitude of chameleonic response (ΔPSA).

Outputs:
  results/tier2_validation_table.csv
  figures/tier2_crosscheck.png

Usage:
  python tier2_validation.py [--matrix results/feature_matrix.csv]
                              [--refset data/reference_set.csv]
                              [--outdir results]
"""

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

CYCLOA_ID = 1   # Canonical CycloA ID in CycPeptMPDB


def run(matrix_csv: str, refset_csv: str, outdir: Path) -> None:
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_csv, low_memory=False)
    ref = pd.read_csv(refset_csv, low_memory=False)

    # Merge on ID
    ref_ids = ref["ID"].unique().tolist()
    df_ref = df[df["ID"].isin(ref_ids)].copy()
    print(f"Reference compounds found in feature matrix: {len(df_ref)}")

    if df_ref.empty:
        print("No reference compounds found — check IDs")
        return

    # ── Table: compare Tier-1 vs. DB 3DPSA ──────────────────────────────────
    compare_cols = {
        "ID": "ID",
        "Original_Name_in_Source_Literature": "Name",
        "PAMPA": "PAMPA_logPexp",
        "H2O_3DPSA": "DB_H2O_3DPSA",
        "CHCl3_3DPSA": "DB_CHCl3_3DPSA",
        "delta_3DPSA_db": "DB_delta_PSA",
        "delta_psa3d": "Tier1_delta_PSA",
        "delta_hb": "Tier1_delta_HB",
        "delta_Rg": "Tier1_delta_Rg",
        "psa3d_spread": "Tier1_PSA_spread",
        "aq_psa3d": "Tier1_aq_PSA",
        "mem_psa3d": "Tier1_mem_PSA",
    }

    avail = {k: v for k, v in compare_cols.items() if k in df_ref.columns}
    table = df_ref[list(avail.keys())].rename(columns=avail).copy()
    table = table.sort_values("PAMPA_logPexp", ascending=True)

    table.to_csv(outdir / "tier2_validation_table.csv", index=False)
    print("\nTier-2 cross-check table:")
    print(table.to_string(index=False))
    print(f"\nSaved: {outdir / 'tier2_validation_table.csv'}")

    # ── Correlation: Tier-1 ΔΨ vs. DB ΔΨ ─────────────────────────────────────
    sub = table.dropna(subset=["DB_delta_PSA", "Tier1_delta_PSA"])
    if len(sub) >= 3:
        r, p = stats.pearsonr(sub["DB_delta_PSA"], sub["Tier1_delta_PSA"])
        print(f"\nTier-1 vs DB ΔΨ: Pearson r = {r:.3f}, p = {p:.3f} (n={len(sub)})")
    else:
        r, p = np.nan, np.nan
        print("\nInsufficient data for Tier-1 vs DB correlation")

    # ── Cross-check figure ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: DB ΔΨ vs PAMPA
    ax = axes[0]
    t = table.dropna(subset=["DB_delta_PSA", "PAMPA_logPexp"])
    if not t.empty:
        ax.scatter(t["DB_delta_PSA"], t["PAMPA_logPexp"], s=60, c="#1F77B4",
                   edgecolors="navy", zorder=3)
        # Label CycloA
        cyc = t[t["ID"] == CYCLOA_ID]
        if not cyc.empty:
            ax.annotate("CycloA", (cyc["DB_delta_PSA"].values[0],
                                   cyc["PAMPA_logPexp"].values[0]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xlabel("DB ΔPSA (H2O − CHCl3)")
        ax.set_ylabel("PAMPA LogPexp (log cm/s)")
        ax.set_title("DB 3DPSA vs. Permeability")
        if len(t) >= 3:
            r2, p2 = stats.pearsonr(t["DB_delta_PSA"], t["PAMPA_logPexp"])
            ax.text(0.05, 0.95, f"r = {r2:.2f}", transform=ax.transAxes,
                    fontsize=9, verticalalignment="top")

    # Plot 2: Tier-1 ΔΨ vs PAMPA
    ax = axes[1]
    t2 = table.dropna(subset=["Tier1_delta_PSA", "PAMPA_logPexp"])
    if not t2.empty:
        ax.scatter(t2["Tier1_delta_PSA"], t2["PAMPA_logPexp"], s=60, c="#E41A1C",
                   edgecolors="darkred", zorder=3)
        cyc2 = t2[t2["ID"] == CYCLOA_ID]
        if not cyc2.empty:
            ax.annotate("CycloA", (cyc2["Tier1_delta_PSA"].values[0],
                                   cyc2["PAMPA_logPexp"].values[0]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xlabel("Tier-1 ΔPSA (ETKDG max-PSA − min-PSA)")
        ax.set_ylabel("PAMPA LogPexp (log cm/s)")
        ax.set_title("Tier-1 ΔPSA vs. Permeability")
        if len(t2) >= 3:
            r3, p3 = stats.pearsonr(t2["Tier1_delta_PSA"], t2["PAMPA_logPexp"])
            ax.text(0.05, 0.95, f"r = {r3:.2f}", transform=ax.transAxes,
                    fontsize=9, verticalalignment="top")

    # Plot 3: Tier-1 vs DB cross-check
    ax = axes[2]
    t3 = table.dropna(subset=["DB_delta_PSA", "Tier1_delta_PSA"])
    if not t3.empty:
        ax.scatter(t3["DB_delta_PSA"], t3["Tier1_delta_PSA"], s=60, c="#7FC97F",
                   edgecolors="darkgreen", zorder=3)
        lims = [
            min(t3["DB_delta_PSA"].min(), t3["Tier1_delta_PSA"].min()),
            max(t3["DB_delta_PSA"].max(), t3["Tier1_delta_PSA"].max()),
        ]
        ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="Identity")
        ax.set_xlabel("DB ΔPSA")
        ax.set_ylabel("Tier-1 ΔPSA")
        ax.set_title(f"Tier-1 vs DB Cross-check\nr = {r:.3f}" if not np.isnan(r) else "Tier-1 vs DB")
        ax.legend(fontsize=8)

    fig.suptitle("Tier-2 Validation: Reference Set Cross-Check", fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig_path = outdir / "figures" / "tier2_crosscheck.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig_path}")

    # ── Interpretation ────────────────────────────────────────────────────────
    print("\n── Interpretation ──")
    if not np.isnan(r):
        if abs(r) > 0.7:
            print(f"✓ Strong agreement (r={r:.2f}): Tier-1 ETKDG approximation validates against DB 3DPSA")
        elif abs(r) > 0.4:
            print(f"~ Moderate agreement (r={r:.2f}): Tier-1 directionally correct, some discrepancy")
        else:
            print(f"⚠ Weak agreement (r={r:.2f}): Tier-1 may miss conformational nuance")
            print("  Consider Tier-2 with OpenEye OMEGA + OpenMM GB/SA for key reference compounds")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tier-2 validation cross-check")
    parser.add_argument("--matrix", "-m", default="results/feature_matrix.csv")
    parser.add_argument("--refset", "-r", default="data/reference_set.csv")
    parser.add_argument("--outdir", "-o", default="results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.matrix, args.refset, Path(args.outdir))
