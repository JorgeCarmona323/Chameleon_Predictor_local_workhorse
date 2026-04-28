"""
plot_tier2_results.py
---------------------
Local analysis and plotting for tier2_crest.py output.
Run after pulling results/tier2_crest_table.csv from the cluster.

Usage:
  python scripts/plot_tier2_results.py --csv results/tier2_crest_table.csv
  python scripts/plot_tier2_results.py --csv results/tier2_crest_table.csv --outdir results
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats


def load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    required = {"compound", "pampa", "permeable", "crest_delta_psa"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")
    return df


def plot(df: pd.DataFrame, outdir: Path) -> None:
    valid = df[df["crest_delta_psa"].notna()].copy()
    if valid.empty:
        print("No compounds with crest_delta_psa — nothing to plot.")
        return

    names   = valid["compound"].tolist()
    c_dpsa  = valid["crest_delta_psa"].tolist()
    pampa   = valid["pampa"].tolist()
    c_dhb   = valid["crest_delta_hb"].tolist() if "crest_delta_hb" in valid.columns else [None] * len(valid)
    db_dpsa = valid["db_delta_psa"].tolist() if "db_delta_psa" in valid.columns else [None] * len(valid)
    colors  = ["#D6604D" if p else "#4393C3" for p in valid["permeable"]]
    n       = len(valid)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(
        "Tier-2 CREST+ALPB Validation — Dual-Dielectric Conformer Sampling\n"
        "Water (ε=80) vs CHCl₃ (ε=4.8)  |  GFN2-xTB + Boltzmann weighting",
        fontsize=11, fontweight="bold",
    )

    # ── Panel A: ΔPSA comparison — DB static vs CREST ensemble ───────────────
    ax = axes[0, 0]
    x = np.arange(n)
    w = 0.35
    db_vals = [v if v is not None and not (isinstance(v, float) and np.isnan(v)) else 0.0
               for v in db_dpsa]
    db_has  = [v is not None and not (isinstance(v, float) and np.isnan(v))
               for v in db_dpsa]
    bars_db   = ax.bar(x - w/2, db_vals, width=w, label="DB static (CycPeptMPDB)",
                       color="#BEAED4", edgecolor="grey", linewidth=0.5)
    bars_crest = ax.bar(x + w/2, c_dpsa, width=w, label="CREST+ALPB ensemble",
                        color="#FDC086", edgecolor="grey", linewidth=0.5)
    # Hatch bars where DB value is unavailable
    for bar, has in zip(bars_db, db_has):
        if not has:
            bar.set_hatch("//")
            bar.set_alpha(0.3)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([n.split("(")[0].strip() for n in names], fontsize=8)
    ax.set_ylabel("ΔPSA (Å²) = PSA_aq − PSA_mem")
    ax.set_title("A. ΔPSA: DB static vs CREST+ALPB ensemble", fontweight="bold")
    ax.legend(fontsize=8)

    # ── Panel B: PAMPA vs CREST ΔPSA ─────────────────────────────────────────
    ax = axes[0, 1]
    for dpsa, pmp, col, name in zip(c_dpsa, pampa, colors, names):
        ax.scatter(dpsa, pmp, s=120, c=col, edgecolors="black", linewidths=0.8, zorder=4)
        ax.annotate(name.split("(")[0].strip(), (dpsa, pmp),
                    xytext=(5, 4), textcoords="offset points", fontsize=7.5)
    ax.axhline(-6.0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("CREST ΔPSA (Å²)")
    ax.set_ylabel("PAMPA (log cm/s)")
    ax.set_title("B. PAMPA vs CREST ΔPSA", fontweight="bold")
    ax.legend(handles=[
        mpatches.Patch(facecolor="#D6604D", label="Permeable"),
        mpatches.Patch(facecolor="#4393C3", label="Impermeable"),
        plt.Line2D([0], [0], color="grey", linestyle="--", label="−6.0 threshold"),
    ], fontsize=8)

    # ── Panel C: PAMPA vs CREST ΔHB ──────────────────────────────────────────
    ax = axes[1, 0]
    if any(v is not None for v in c_dhb):
        for dhb, pmp, col, name in zip(c_dhb, pampa, colors, names):
            if dhb is None:
                continue
            ax.scatter(dhb, pmp, s=120, c=col, edgecolors="black", linewidths=0.8, zorder=4)
            ax.annotate(name.split("(")[0].strip(), (dhb, pmp),
                        xytext=(5, 4), textcoords="offset points", fontsize=7.5)
        ax.axhline(-6.0, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel("CREST ΔHB (H-bonds mem − H-bonds aq)")
        ax.set_ylabel("PAMPA (log cm/s)")
        ax.set_title("C. PAMPA vs CREST ΔHB\n(intramolecular H-bond formation in membrane)",
                     fontweight="bold")
    else:
        ax.text(0.5, 0.5, "ΔHB data not available", transform=ax.transAxes,
                ha="center", va="center", color="grey")
        ax.set_title("C. PAMPA vs CREST ΔHB", fontweight="bold")

    # ── Panel D: CREST ΔPSA vs DB ΔPSA cross-check ───────────────────────────
    ax = axes[1, 1]
    paired = [(db, cr, col, name)
              for db, cr, col, name in zip(db_dpsa, c_dpsa, colors, names)
              if db is not None and not (isinstance(db, float) and np.isnan(db))]
    ax.set_xlabel("DB ΔPSA (static CycPeptMPDB)")
    ax.set_ylabel("CREST ΔPSA (ensemble, dual-dielectric)")
    ax.set_title("D. CREST vs DB Cross-Check\n(DB static misses CsA chameleonism)",
                 fontweight="bold")
    if paired:
        db_v, cr_v, col_v, name_v = zip(*paired)
        for db, cr, col, name in zip(db_v, cr_v, col_v, name_v):
            ax.scatter(db, cr, s=120, c=col, edgecolors="black", linewidths=0.8, zorder=4)
            ax.annotate(name.split("(")[0].strip(), (db, cr),
                        xytext=(5, 4), textcoords="offset points", fontsize=7.5)
        all_v = list(db_v) + list(cr_v)
        lim = [min(all_v) - 5, max(all_v) + 5]
        ax.plot(lim, lim, "k--", linewidth=0.8, alpha=0.5, label="y=x")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.legend(fontsize=8)
        if len(paired) >= 3:
            try:
                r_val, p_val = stats.pearsonr(db_v, cr_v)
                ax.text(0.05, 0.95, f"r = {r_val:.2f}  (p={p_val:.3f})",
                        transform=ax.transAxes, fontsize=9, va="top", fontweight="bold")
            except Exception:
                pass
    else:
        ax.text(0.5, 0.5, "No DB ΔPSA available for cross-check",
                transform=ax.transAxes, ha="center", va="center", color="grey")

    plt.tight_layout()
    out_path = outdir / "figures" / "tier2_crest_crosscheck.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def print_summary(df: pd.DataFrame) -> None:
    print(f"\n{'='*65}")
    print("Tier-2 CREST Summary")
    print(f"{'='*65}")
    cols = ["compound", "pampa", "permeable", "db_delta_psa",
            "crest_delta_psa", "crest_delta_psa_lowen",
            "crest_delta_hb", "aq_n_confs", "mem_n_confs"]
    avail = [c for c in cols if c in df.columns]
    print(df[avail].to_string(index=False))

    # CsA check
    csa = df[df["compound"].str.contains("Cyclosporin|CsA", case=False, na=False)]
    if not csa.empty and "crest_delta_psa" in csa.columns:
        val = csa["crest_delta_psa"].iloc[0]
        status = "PASS ✓" if val is not None and val > 50 else "CHECK — expected ~75 Å²"
        print(f"\nCsA ΔPSA validation: {val:.1f} Å²  →  {status}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot tier2_crest.py results locally")
    p.add_argument("--csv",    required=True, type=Path,
                   help="Path to tier2_crest_table.csv")
    p.add_argument("--outdir", default="results", type=Path,
                   help="Output directory for figures (default: results/)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    df = load(args.csv)
    print_summary(df)
    plot(df, args.outdir)
