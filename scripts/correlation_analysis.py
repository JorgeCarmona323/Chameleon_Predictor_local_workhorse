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
        "delta_psa3d", "delta_hb", "delta_Rg",
        "delta_NPR1", "delta_NPR2", "delta_Asphericity",
        "psa3d_spread", "psa3d_std", "hb_spread",
    ],
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
    plt.tight_layout()

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
    palette = {"2D Baseline": "#7FC97F", "DB 3D (from CycPeptMPDB)": "#BEAED4",
               "Tier-1 Δ (conformer engine)": "#FDC086"}
    for g in top["Group"]:
        colors.append(palette.get(g, "#CCCCCC"))

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
    legend_elements = [Patch(facecolor=c, label=g) for g, c in palette.items()]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    path = outdir / "figures" / "auc_roc_bar.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def plot_scatter_top(df: pd.DataFrame, corr_df: pd.DataFrame, outdir: Path,
                     n_top: int = 6) -> None:
    """Scatter plot of top-N features vs. PAMPA."""
    valid = corr_df.dropna(subset=["Spearman_rho"])
    top_feats = valid.head(n_top)["Feature"].tolist()

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

    fig.suptitle("Top Features vs. PAMPA LogPexp", fontsize=12, fontweight="bold")
    plt.tight_layout()

    path = outdir / "figures" / "scatter_top_features.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def run(matrix_csv: str, outdir: Path) -> None:
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    df = load_and_validate(matrix_csv)

    # Merge feature groups from JSON if it exists
    fg_path = outdir / "feature_groups.json"
    if fg_path.exists():
        fg_raw = json.loads(fg_path.read_text())
        feature_groups = {k: v for k, v in fg_raw.items() if k != "combined"}
    else:
        feature_groups = FEATURE_GROUPS

    print("\n── Correlation analysis ──")
    corr_df = compute_correlations(df, feature_groups)
    corr_df.to_csv(outdir / "correlation_table.csv", index=False)
    print(corr_df.head(15).to_string(index=False))
    print(f"\nSaved: {outdir / 'correlation_table.csv'}")

    print("\n── AUC-ROC analysis ──")
    auc_df = compute_auc_roc(df, feature_groups)
    auc_df.to_csv(outdir / "auc_roc_table.csv", index=False)
    print(auc_df.head(15).to_string(index=False))
    print(f"\nSaved: {outdir / 'auc_roc_table.csv'}")

    print("\n── Logistic regression (combined features) ──")
    all_features = [f for grp in feature_groups.values() for f in grp]
    lr_df = logistic_regression_importance(df, all_features)
    if not lr_df.empty:
        lr_df.to_csv(outdir / "feature_importance.csv", index=False)
        print(lr_df.head(10).to_string(index=False))

    print("\n── Generating figures ──")
    plot_correlation_heatmap(corr_df, df, outdir)
    plot_auc_bar(auc_df, outdir)
    plot_scatter_top(df, corr_df, outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Correlate Δ features vs. PAMPA")
    parser.add_argument("--matrix", "-m", default="results/feature_matrix.csv")
    parser.add_argument("--outdir", "-o", default="results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.matrix, Path(args.outdir))
