"""
03_correlation_analysis.py
--------------------------
Correlate 3D Δ features vs. experimental PAMPA LogPexp.

Outputs (to results/):
  correlation_table.csv     — Pearson r, Spearman ρ, p-values per feature
  auc_roc_table.csv         — AUC-ROC per feature (binarized at PAMPA >= -6.0)
  feature_importance.csv    — Logistic regression coefficients (combined features)
  figures/correlation_heatmap.png
  figures/auc_roc_bar.png
  figures/scatter_top_features.png

Usage:
  python correlation_analysis.py [--matrix results/feature_matrix.csv]
                                  [--outdir results]
"""

import argparse
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

PAMPA_THRESHOLD = -6.0  # Jiang et al. 2023, J. Chem. Inf. Model. — CycPeptMPDB standard cutoff

# Features to evaluate (must exist in feature_matrix.csv)
FEATURE_GROUPS = {
    "2D Baseline": [
        "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
        "NumRotatableBonds", "FractionCSP3",
    ],
    "DB 3D (from CycPeptMPDB)": [
        "delta_3DPSA_db", "H2O_3DPSA", "CHCl3_3DPSA",
    ],
    "Tier-1 Δ (conformer engine)": [
        "delta_psa3d", "norm_delta_psa", "delta_hb", "delta_Rg",
        "delta_NPR1", "delta_NPR2", "delta_Asphericity",
        "psa3d_std",
    ],
    "Tier-1 Absolute (conformer engine)": [
        "aq_psa3d", "mem_psa3d",
        "aq_Rg", "mem_Rg",
        "aq_hb_count", "mem_hb_count",
        "aq_NPR1", "mem_NPR1", "aq_NPR2", "mem_NPR2",
    ],
}

# Map JSON short keys → display palette keys
GROUP_NAME_MAP = {
    "2D_baseline":  "2D Baseline",
    "DB_delta":     "DB 3D (from CycPeptMPDB)",
    "Tier1_delta":  "Tier-1 Δ (conformer engine)",
}

GROUP_PALETTE = {
    "2D Baseline":              "#7FC97F",   # green
    "DB 3D (from CycPeptMPDB)": "#BEAED4",   # purple
    "Tier-1 Δ (conformer engine)": "#FDC086", # orange
}


def load_and_validate(matrix_csv: str) -> pd.DataFrame:
    df = pd.read_csv(matrix_csv, low_memory=False)
    df = df[df["PAMPA"].notna()].copy()
    print(f"Loaded {len(df)} compounds with PAMPA values")
    df["permeable"] = (df["PAMPA"] >= PAMPA_THRESHOLD).astype(int)
    print(f"  Permeable: {df['permeable'].sum()} ({100*df['permeable'].mean():.1f}%)")
    return df


def compute_correlations(df: pd.DataFrame, feature_groups: dict) -> pd.DataFrame:
    """Pearson r and Spearman ρ of each feature vs. PAMPA (LogPexp)."""
    rows = []
    y = df["PAMPA"].values

    for group, features in feature_groups.items():
        for feat in features:
            if feat not in df.columns:
                continue
            mask = df[feat].notna() & df["PAMPA"].notna()
            x = df.loc[mask, feat].values
            y_sub = df.loc[mask, "PAMPA"].values
            n = mask.sum()

            if n < 10:
                rows.append({"Feature": feat, "Group": group, "N": n,
                             "Pearson_r": np.nan, "Pearson_p": np.nan,
                             "Spearman_rho": np.nan, "Spearman_p": np.nan})
                continue

            pr, pp = stats.pearsonr(x, y_sub)
            sr, sp = stats.spearmanr(x, y_sub)
            rows.append({
                "Feature": feat, "Group": group, "N": n,
                "Pearson_r": round(pr, 4), "Pearson_p": round(pp, 4),
                "Spearman_rho": round(sr, 4), "Spearman_p": round(sp, 4),
            })

    corr_df = pd.DataFrame(rows)
    corr_df["abs_Spearman"] = corr_df["Spearman_rho"].abs()
    corr_df = corr_df.sort_values("abs_Spearman", ascending=False)
    return corr_df


def compute_auc_roc(df: pd.DataFrame, feature_groups: dict) -> pd.DataFrame:
    """AUC-ROC for binarized PAMPA >= PAMPA_THRESHOLD."""
    rows = []
    y = df["permeable"].values

    for group, features in feature_groups.items():
        for feat in features:
            if feat not in df.columns:
                continue
            mask = df[feat].notna() & df["permeable"].notna()
            x = df.loc[mask, feat].values
            y_sub = df.loc[mask, "permeable"].values
            n = mask.sum()

            if n < 10 or y_sub.sum() == 0 or y_sub.sum() == n:
                # roc_auc_score raises ValueError on single-class input — guard explicitly
                rows.append({"Feature": feat, "Group": group, "N": n, "AUC_ROC": np.nan})
                continue

            # AUC for x and -x (take max, so directionality doesn't matter)
            try:
                auc = max(
                    roc_auc_score(y_sub, x),
                    roc_auc_score(y_sub, -x),
                )
            except Exception:
                auc = np.nan

            rows.append({"Feature": feat, "Group": group, "N": n, "AUC_ROC": round(auc, 4)})

    return pd.DataFrame(rows).sort_values("AUC_ROC", ascending=False)


def logistic_regression_importance(df: pd.DataFrame, features: list) -> pd.DataFrame:
    """Logistic regression on combined feature set → standardized coefficients.

    Scaling rationale:
      StandardScaler (mean=0, std=1) is used here because logistic regression
      coefficients are only comparable across features when they are on the same
      scale. RobustScaler would work too but StandardScaler is more interpretable
      for coefficient comparison. The raw correlation analysis (Pearson/Spearman)
      does not require scaling since it is rank/linear-based.
    """
    available = [f for f in features if f in df.columns]
    if not available:
        return pd.DataFrame()

    sub = df[available + ["permeable"]].dropna()
    if len(sub) < 50:
        print("  Insufficient data for logistic regression")
        return pd.DataFrame()

    X = sub[available].values
    y = sub["permeable"].values

    # StandardScaler: features need same scale for LR coefficient comparison
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    lr = LogisticRegression(C=1.0, max_iter=1000, random_state=42, solver="lbfgs")
    lr.fit(X_scaled, y)

    coef_df = pd.DataFrame({
        "Feature": available,
        "LR_coefficient": lr.coef_[0].round(4),
        "abs_coef": np.abs(lr.coef_[0]).round(4),
    }).sort_values("abs_coef", ascending=False)

    return coef_df


def plot_correlation_heatmap(corr_df: pd.DataFrame, df: pd.DataFrame,
                              outdir: Path) -> None:
    """Heatmap of Pearson r and Spearman ρ side-by-side."""
    valid = corr_df.dropna(subset=["Pearson_r", "Spearman_rho"])
    if valid.empty:
        return

    top_n = min(20, len(valid))
    top_feats = valid.head(top_n)["Feature"].tolist()

    r_vals    = valid.set_index("Feature").loc[top_feats, "Pearson_r"].values
    rho_vals  = valid.set_index("Feature").loc[top_feats, "Spearman_rho"].values

    data = np.column_stack([r_vals, rho_vals])
    fig, ax = plt.subplots(figsize=(5, max(4, top_n * 0.4 + 1)))

    sns.heatmap(
        data,
        annot=True, fmt=".3f",
        xticklabels=["Pearson r", "Spearman ρ"],
        yticklabels=top_feats,
        cmap="coolwarm", center=0, vmin=-1, vmax=1,
        linewidths=0.3, ax=ax,
        annot_kws={"size": 8},
    )
    ax.set_title(f"Feature vs. PAMPA LogPexp Correlation (top {top_n})", fontsize=11, fontweight="bold")
    ax.set_xlabel("Correlation metric")

    # Color y-tick labels by feature group
    feat_to_group = valid.set_index("Feature")["Group"].to_dict()
    for tick, feat in zip(ax.get_yticklabels(), top_feats):
        grp = feat_to_group.get(feat, "")
        tick.set_color(GROUP_PALETTE.get(grp, "black"))
        tick.set_fontweight("bold")

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=g) for g, c in GROUP_PALETTE.items()]
    fig.legend(handles=legend_elements, loc="lower center", fontsize=7,
               ncol=3, bbox_to_anchor=(0.5, -0.02), frameon=True)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)

    path = outdir / "figures" / "correlation_heatmap.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def plot_auc_bar(auc_df: pd.DataFrame, outdir: Path) -> None:
    valid = auc_df.dropna(subset=["AUC_ROC"])
    if valid.empty:
        return

    top_n = min(20, len(valid))
    top = valid.head(top_n)

    colors = []
    for g in top["Group"]:
        colors.append(GROUP_PALETTE.get(g, "#CCCCCC"))

    fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.35 + 1)))
    bars = ax.barh(range(top_n), top["AUC_ROC"].values, color=colors, edgecolor="grey", linewidth=0.5)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=0.8, label="Random (0.5)")
    ax.set_yticks(range(top_n))
    ax.set_yticklabels(top["Feature"].tolist(), fontsize=8)
    ax.set_xlabel("AUC-ROC")
    ax.set_title("Feature AUC-ROC vs. PAMPA Permeability\n(threshold = −6.0 log cm/s)", fontsize=11)
    ax.invert_yaxis()

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=g) for g, c in GROUP_PALETTE.items()]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    path = outdir / "figures" / "auc_roc_bar.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def plot_scatter_top(df: pd.DataFrame, corr_df: pd.DataFrame, outdir: Path,
                     n_top: int = 3) -> None:
    """Scatter plot of top-N and bottom-N features vs. PAMPA (by Spearman ρ)."""
    valid = corr_df.dropna(subset=["Spearman_rho"])
    top_feats = valid.head(n_top)["Feature"].tolist()
    bot_feats = valid.tail(n_top)["Feature"].tolist()
    top_feats = top_feats + bot_feats
    n_top = len(top_feats)

    cols = min(3, n_top)
    rows = (n_top + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3.5))
    axes = axes.flatten() if n_top > 1 else [axes]

    permeable_mask = df["permeable"] == 1

    for i, feat in enumerate(top_feats):
        ax = axes[i]
        if feat not in df.columns:
            ax.set_visible(False)
            continue

        sub = df[[feat, "PAMPA", "permeable"]].dropna()
        rho = corr_df.set_index("Feature").loc[feat, "Spearman_rho"]

        ax.scatter(sub.loc[sub["permeable"] == 0, feat],
                   sub.loc[sub["permeable"] == 0, "PAMPA"],
                   s=8, alpha=0.3, c="#1F77B4", label="Non-permeable", rasterized=True)
        ax.scatter(sub.loc[sub["permeable"] == 1, feat],
                   sub.loc[sub["permeable"] == 1, "PAMPA"],
                   s=15, alpha=0.6, c="#E41A1C", label="Permeable", rasterized=True)

        ax.axhline(PAMPA_THRESHOLD, color="grey", linestyle="--", linewidth=0.8)
        ax.set_xlabel(feat, fontsize=8)
        ax.set_ylabel("PAMPA (log cm/s)", fontsize=8)
        ax.set_title(f"ρ = {rho:.3f}", fontsize=9, fontweight="bold")

        if i == 0:
            ax.legend(fontsize=7, markerscale=1.5)

    for j in range(len(top_feats), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Top 3 vs. Bottom 3 Features — PAMPA LogPexp (Spearman ρ)", fontsize=12, fontweight="bold")
    plt.tight_layout()

    path = outdir / "figures" / "scatter_top_features.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def run(matrix_csv: str, outdir: Path, sources: list = None,
        min_residues: int = 0) -> None:
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    df = load_and_validate(matrix_csv)

    # ── Source filter ─────────────────────────────────────────────────────────
    if sources:
        before = len(df)
        df = df[df["Source"].isin(sources)].copy()
        print(f"Source filter {sources}: {len(df)} / {before} compounds")

    # ── Residue size filter ───────────────────────────────────────────────────
    if min_residues > 0 and "Monomer_Length" in df.columns:
        before = len(df)
        df = df[df["Monomer_Length"] >= min_residues].copy()
        print(f"Monomer_Length >= {min_residues}: {len(df)} / {before} compounds")

    # ── Compute norm_delta_psa inline (ΔPSA / aq_psa3d) ─────────────────────
    # Approximation of Yu et al. 2026 (ΔPSA / SASA_total): uses aq_psa3d as
    # denominator since total SASA is not in feature_matrix.csv.
    # Captures fractional burial of polar surface — removes MW confounding.
    if "delta_psa3d" in df.columns and "aq_psa3d" in df.columns:
        df["norm_delta_psa"] = (df["delta_psa3d"] / df["aq_psa3d"]).replace(
            [np.inf, -np.inf], np.nan
        )

    # ── Build suffix for output filenames ────────────────────────────────────
    suffix_parts = []
    if sources:
        suffix_parts.append("_".join(s.split("_")[0] for s in sources))
    if min_residues > 0:
        suffix_parts.append(f"res{min_residues}plus")
    suffix = f"_{len(df)}" if not suffix_parts else f"_{'_'.join(suffix_parts)}_{len(df)}"

    # Merge feature groups from JSON if it exists
    fg_path = outdir / "feature_groups.json"
    if fg_path.exists():
        fg_raw = json.loads(fg_path.read_text())
        DROP = {"psa3d_spread", "hb_spread"}
        feature_groups = {}
        for k, v in fg_raw.items():
            if k == "combined":
                continue
            display = GROUP_NAME_MAP.get(k, k)
            feature_groups[display] = [f for f in v if f not in DROP]
        # Inject new groups not in JSON
        feature_groups["Tier-1 Δ (conformer engine)"] = list(dict.fromkeys(
            feature_groups.get("Tier-1 Δ (conformer engine)", []) + ["norm_delta_psa"]
        ))
    else:
        feature_groups = FEATURE_GROUPS

    # Add absolute conformer group if columns present
    abs_cols = ["aq_psa3d", "mem_psa3d", "aq_Rg", "mem_Rg",
                "aq_hb_count", "mem_hb_count", "aq_NPR1", "mem_NPR1",
                "aq_NPR2", "mem_NPR2"]
    if any(c in df.columns for c in abs_cols):
        feature_groups["Tier-1 Absolute (conformer engine)"] = abs_cols

    n_label = f"n={len(df)}"
    if sources:
        n_label += f"  sources={sources}"
    if min_residues > 0:
        n_label += f"  ≥{min_residues} residues"

    print(f"\n── Correlation analysis ({n_label}) ──")
    corr_df = compute_correlations(df, feature_groups)
    corr_csv = outdir / f"correlation_table{suffix}.csv"
    corr_df.to_csv(corr_csv, index=False)
    print(corr_df.head(15).to_string(index=False))

    print(f"\n── AUC-ROC analysis ({n_label}) ──")
    auc_df = compute_auc_roc(df, feature_groups)
    auc_csv = outdir / f"auc_roc_table{suffix}.csv"
    auc_df.to_csv(auc_csv, index=False)
    print(auc_df.head(15).to_string(index=False))

    print("\n── Logistic regression (combined features) ──")
    all_features = [f for grp in feature_groups.values() for f in grp]
    lr_df = logistic_regression_importance(df, all_features)
    if not lr_df.empty:
        lr_df.to_csv(outdir / f"feature_importance{suffix}.csv", index=False)
        print(lr_df.head(10).to_string(index=False))

    print("\n── Generating figures ──")
    # Patch outdir for figures so subset runs don't overwrite full-dataset PNGs
    fig_outdir = outdir / "figures"

    def _suffixed_save(fig, basename):
        path = fig_outdir / f"{basename}{suffix}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")
        return path

    # Heatmap
    valid = corr_df.dropna(subset=["Pearson_r", "Spearman_rho"])
    if not valid.empty:
        top_n = min(20, len(valid))
        top_feats = valid.head(top_n)["Feature"].tolist()
        r_vals   = valid.set_index("Feature").loc[top_feats, "Pearson_r"].values
        rho_vals = valid.set_index("Feature").loc[top_feats, "Spearman_rho"].values
        data = np.column_stack([r_vals, rho_vals])
        fig, ax = plt.subplots(figsize=(5, max(4, top_n * 0.4 + 1)))
        import seaborn as sns
        sns.heatmap(data, annot=True, fmt=".3f",
                    xticklabels=["Pearson r", "Spearman ρ"],
                    yticklabels=top_feats, cmap="coolwarm", center=0,
                    vmin=-1, vmax=1, linewidths=0.3, ax=ax, annot_kws={"size": 8})
        feat_to_group = valid.set_index("Feature")["Group"].to_dict()
        for tick, feat in zip(ax.get_yticklabels(), top_feats):
            tick.set_color(GROUP_PALETTE.get(feat_to_group.get(feat, ""), "black"))
            tick.set_fontweight("bold")
        ax.set_title(f"Feature vs. PAMPA Correlation ({n_label})", fontsize=10,
                     fontweight="bold")
        plt.tight_layout()
        _suffixed_save(fig, "correlation_heatmap")

    # AUC bar — this is the key slide figure
    valid_auc = auc_df.dropna(subset=["AUC_ROC"])
    if not valid_auc.empty:
        top_n = min(20, len(valid_auc))
        top = valid_auc.head(top_n)
        colors = [GROUP_PALETTE.get(g, "#CCCCCC") for g in top["Group"]]
        fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.35 + 1)))
        ax.barh(range(top_n), top["AUC_ROC"].values, color=colors,
                edgecolor="grey", linewidth=0.5)
        ax.axvline(0.5, color="black", linestyle="--", linewidth=0.8,
                   label="Random (0.5)")
        ax.set_yticks(range(top_n))
        ax.set_yticklabels(top["Feature"].tolist(), fontsize=8)
        ax.set_xlabel("AUC-ROC")
        ax.set_title(f"Feature AUC-ROC vs. PAMPA  |  {n_label}\n"
                     f"(threshold = −6.0 log cm/s)", fontsize=11)
        ax.invert_yaxis()
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=c, label=g)
                           for g, c in GROUP_PALETTE.items()]
        ax.legend(handles=legend_elements, loc="lower right", fontsize=8)
        plt.tight_layout()
        auc_bar_path = fig_outdir / f"auc_roc_bar{suffix}.png"
        fig.savefig(auc_bar_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {auc_bar_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Correlate Δ features vs. PAMPA")
    parser.add_argument("--matrix",  "-m", default="results/feature_matrix.csv")
    parser.add_argument("--outdir",  "-o", default="results")
    parser.add_argument("--sources", "-s", nargs="+", default=None,
                        help="Filter to these Source values (e.g. 2016_Furukawa 2013_CHUGAI)")
    parser.add_argument("--min-residues", "-r", type=int, default=0,
                        help="Keep only compounds with Monomer_Length >= this value")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.matrix, Path(args.outdir),
        sources=args.sources,
        min_residues=args.min_residues)
