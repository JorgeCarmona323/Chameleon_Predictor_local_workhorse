"""
04_umap_visualization.py
------------------------
UMAP of CycPeptMPDB PAMPA subset colored by LogPexp (continuous) and
permeability class (binary), with CycloA reference overlay.

Two UMAP panels:
  Panel A — 2D descriptors (MolWt, MolLogP, TPSA, HBA, HBD, RotBonds, CSP3)
  Panel B — 3D Δ features  (delta_3DPSA_db, delta_psa3d, delta_hb, delta_Rg,
                              delta_NPR1, delta_NPR2, psa3d_spread)

For each panel: cosine metric, Leiden clustering on UMAP kNN graph.
Hit enrichment per cluster reported.

Usage:
  python umap_visualization.py [--matrix results/feature_matrix.csv]
                                [--outdir results]
"""

import argparse
import warnings
from pathlib import Path

try:
    import igraph as ig
    import leidenalg
    _LEIDEN_AVAILABLE = True
except ImportError:
    _LEIDEN_AVAILABLE = False
    print("WARNING: leidenalg/igraph not available — Leiden clustering will be skipped.")
    print("  Fix: pip install leidenalg igraph --no-binary leidenalg")

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import umap
from matplotlib.colors import Normalize
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import RobustScaler

warnings.filterwarnings("ignore")

PAMPA_THRESHOLD = -6.0  # Jiang et al. 2023, J. Chem. Inf. Model. — CycPeptMPDB standard cutoff
RANDOM_STATE = 42

# UMAP parameters — optimized for macrocycle chemical space
UMAP_PARAMS = dict(
    n_neighbors=30,
    min_dist=0.15,
    n_components=2,
    metric="cosine",
    random_state=RANDOM_STATE,
    low_memory=False,
)

FEATURE_PANELS = {
    "Panel_A_2D": [
        "MolWt", "MolLogP", "TPSA", "NumHAcceptors",
        "NumHDonors", "NumRotatableBonds", "FractionCSP3", "RingCount",
    ],
    "Panel_B_3D_delta": [
        "delta_3DPSA_db", "delta_psa3d", "delta_hb", "delta_Rg",
        "delta_NPR1", "delta_NPR2", "psa3d_spread", "psa3d_std",
    ],
    "Panel_C_combined": [
        "MolWt", "MolLogP", "TPSA",
        "delta_3DPSA_db", "delta_psa3d", "delta_hb", "delta_Rg",
        "psa3d_spread", "delta_NPR1",
    ],
}

# CycloA IDs in CycPeptMPDB (from reference_set)
CYCLOA_IDS = {1, 22, 932, 981, 1822, 1862, 2356, 7188, 7353}


def run_leiden(reducer: umap.UMAP) -> np.ndarray:
    """Run Leiden clustering on UMAP's internal kNN graph.

    Use upper-triangle only to avoid duplicate edges from the symmetric matrix
    (each undirected edge appears as both (i,j) and (j,i) in the full COO).
    """
    from scipy.sparse import triu
    cx = triu(reducer.graph_).tocoo()
    sources = cx.row.tolist()
    targets = cx.col.tolist()
    weights = cx.data.tolist()
    g = ig.Graph(
        n=reducer.graph_.shape[0],
        edges=list(zip(sources, targets)),
        edge_attrs={"weight": weights},
    )
    partition = leidenalg.find_partition(
        g, leidenalg.ModularityVertexPartition,
        weights="weight", seed=RANDOM_STATE,
    )
    return np.array(partition.membership)


def hit_enrichment(labels: np.ndarray, permeable: np.ndarray) -> pd.DataFrame:
    """Compute hit enrichment ratio per cluster.

    Both labels and permeable must be numpy arrays with positional alignment
    (not pandas Series) to avoid index-alignment bugs.
    """
    permeable = np.asarray(permeable)  # ensure positional indexing
    total_perm = permeable.sum()
    total = len(permeable)
    bg_rate = total_perm / total if total > 0 else 0

    rows = []
    for lab in np.unique(labels):
        mask = labels == lab
        n_cluster = mask.sum()
        n_hits = permeable[mask].sum()
        cluster_rate = n_hits / n_cluster if n_cluster > 0 else 0
        enrichment = cluster_rate / bg_rate if bg_rate > 0 else 0
        rows.append({
            "cluster": int(lab),
            "n_cluster": int(n_cluster),
            "n_permeable": int(n_hits),
            "hit_rate": round(float(cluster_rate), 4),
            "enrichment_ratio": round(float(enrichment), 3),
        })
    return pd.DataFrame(rows).sort_values("enrichment_ratio", ascending=False)


def make_umap_panel(
    df_panel: pd.DataFrame,
    panel_name: str,
    features: list,
    outdir: Path,
) -> dict:
    """Fit UMAP, run Leiden, plot, return metrics dict."""

    available = [f for f in features if f in df_panel.columns and df_panel[f].notna().sum() > 50]
    if not available:
        print(f"  {panel_name}: insufficient features, skipping")
        return {}

    # Drop rows missing any feature
    sub = df_panel[available + ["PAMPA", "permeable", "ID"]].dropna().copy()
    print(f"\n── {panel_name} ──")
    print(f"  Features: {available}")
    print(f"  Compounds: {len(sub)}")

    X = sub[available].values
    y_cont = sub["PAMPA"].values
    y_bin  = sub["permeable"].values

    # ── Scale first — critical before any distance-based method ─────────────
    # RobustScaler: centers to median, scales by IQR.
    # Correct for dense continuous descriptors with different magnitudes
    # (e.g., MolWt ~600-1400 Da vs. MolLogP ~0-8 vs. delta_psa3d ~0-100 Å²).
    # Prevents large-magnitude features from dominating cosine distances.
    scaler = RobustScaler()
    X_scaled = scaler.fit_transform(X)

    # ── PCA dimensionality reduction (only if n_features > 10) ──────────────
    # PCA on already-scaled dense continuous data is the correct choice here.
    # TruncatedSVD is for sparse matrices (fingerprints); PCA is for scaled
    # dense continuous descriptors like our Δ features.
    # Silhouette score is then computed on PCA coords — NOT on 2D UMAP coords,
    # because UMAP distorts inter-cluster distances in the 2D projection.
    if X_scaled.shape[1] > 10:
        n_pca = min(50, X_scaled.shape[1] - 1, X_scaled.shape[0] - 1)
        pca = PCA(n_components=n_pca, random_state=RANDOM_STATE)
        X_reduced = pca.fit_transform(X_scaled)
        var_exp = pca.explained_variance_ratio_.cumsum()[-1]
        print(f"  PCA({n_pca}): {100*var_exp:.1f}% cumulative variance explained")
        # Find elbow: n_components for 90% variance
        n_90 = np.searchsorted(pca.explained_variance_ratio_.cumsum(), 0.90) + 1
        print(f"  PCA components for 90% variance: {n_90}")
    else:
        X_reduced = X_scaled
        print(f"  Skipping PCA ({X_scaled.shape[1]} features < 10 — using scaled features directly)")

    # UMAP fit
    print(f"  Fitting UMAP ...")
    reducer = umap.UMAP(**UMAP_PARAMS)
    embedding = reducer.fit_transform(X_reduced)

    # Leiden clustering on UMAP kNN graph (falls back to single cluster if unavailable)
    if _LEIDEN_AVAILABLE:
        print("  Running Leiden clustering ...")
        leiden_labels = run_leiden(reducer)
        n_clusters = len(np.unique(leiden_labels))
        print(f"  Leiden: {n_clusters} clusters")
    else:
        leiden_labels = np.zeros(len(sub), dtype=int)
        n_clusters = 1
        print("  Leiden unavailable — assigning all compounds to cluster 0")

    # Silhouette on X_reduced (NOT on 2D embedding)
    if n_clusters > 1:
        sil = silhouette_score(X_reduced, leiden_labels, metric="cosine")
        print(f"  Silhouette (cosine, SVD coords): {sil:.4f}")
    else:
        sil = np.nan

    # Hit enrichment
    enrich_df = hit_enrichment(leiden_labels, y_bin)
    enrich_df["panel"] = panel_name
    enrich_df.to_csv(outdir / f"{panel_name}_hit_enrichment.csv", index=False)
    print(f"  Top enriched clusters:\n{enrich_df.head(5).to_string(index=False)}")

    # ── Figure 1: LogPexp continuous colormap ─────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Subplot 1: colored by PAMPA LogPexp
    ax = axes[0]
    norm = Normalize(vmin=np.percentile(y_cont, 5), vmax=np.percentile(y_cont, 95))
    cmap = plt.cm.RdYlGn

    sc = ax.scatter(
        embedding[:, 0], embedding[:, 1],
        c=y_cont, cmap=cmap, norm=norm,
        s=8, alpha=0.5, rasterized=True,
    )
    plt.colorbar(sc, ax=ax, label="PAMPA LogPexp (log cm/s)")

    # CycloA overlay
    cycloA_mask = sub["ID"].isin(CYCLOA_IDS)
    if cycloA_mask.sum() > 0:
        ax.scatter(
            embedding[cycloA_mask, 0], embedding[cycloA_mask, 1],
            s=120, marker="*", c="black", zorder=10,
            edgecolors="white", linewidths=0.8, label=f"CycloA (n={cycloA_mask.sum()})",
        )
        ax.legend(fontsize=9)

    ax.axhline(0, color="grey", linewidth=0.3, alpha=0.5)
    ax.axvline(0, color="grey", linewidth=0.3, alpha=0.5)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"{panel_name}\nColored by PAMPA LogPexp")

    # Subplot 2: colored by Leiden cluster
    ax2 = axes[1]
    cluster_cmap = plt.cm.tab20
    colors_cluster = [cluster_cmap(l % 20) for l in leiden_labels]
    ax2.scatter(
        embedding[:, 0], embedding[:, 1],
        c=colors_cluster, s=8, alpha=0.5, rasterized=True,
    )
    # Overlay permeable compounds
    perm_mask = y_bin == 1
    ax2.scatter(
        embedding[perm_mask, 0], embedding[perm_mask, 1],
        s=35, alpha=0.85, marker="^", c="#E41A1C",
        edgecolors="darkred", linewidths=0.5, zorder=5,
        label=f"Permeable (n={perm_mask.sum()})",
    )
    if cycloA_mask.sum() > 0:
        ax2.scatter(
            embedding[cycloA_mask, 0], embedding[cycloA_mask, 1],
            s=120, marker="*", c="black", zorder=10,
            edgecolors="white", linewidths=0.8, label="CycloA",
        )
    ax2.legend(fontsize=9)
    ax2.set_xlabel("UMAP 1")
    ax2.set_ylabel("UMAP 2")
    ax2.set_title(f"{panel_name}\nLeiden clusters (n={n_clusters})")

    fig.suptitle(
        f"{panel_name} | metric=cosine | Leiden | sil={sil:.3f}" if not np.isnan(sil)
        else f"{panel_name} | metric=cosine | Leiden",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig_path = outdir / "figures" / f"{panel_name}_umap.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fig_path}")

    # ── Figure 2: LogPexp distribution by cluster (violin) ────────────────────
    sub_plot = sub.copy()
    sub_plot["cluster"] = leiden_labels
    sub_plot["embedding_x"] = embedding[:, 0]
    sub_plot["embedding_y"] = embedding[:, 1]

    # Save embedding with metadata
    emb_path = outdir / f"{panel_name}_embedding.csv"
    sub_plot[["ID", "PAMPA", "permeable", "cluster", "embedding_x", "embedding_y"]].to_csv(
        emb_path, index=False
    )

    return {
        "panel": panel_name,
        "n_compounds": len(sub),
        "n_features": len(available),
        "n_clusters": n_clusters,
        "silhouette": round(float(sil), 4) if not np.isnan(sil) else None,
        "max_enrichment": float(enrich_df["enrichment_ratio"].max()),
        "permeable_total": int(perm_mask.sum()),
        "cycloA_found": int(cycloA_mask.sum()),
    }


def run(matrix_csv: str, outdir: Path) -> None:
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_csv, low_memory=False)
    df = df[df["PAMPA"].notna()].copy()
    df["permeable"] = (df["PAMPA"] >= PAMPA_THRESHOLD).astype(int)
    print(f"Loaded {len(df)} compounds with PAMPA values")

    summary_rows = []
    for panel_name, features in FEATURE_PANELS.items():
        result = make_umap_panel(df, panel_name, features, outdir)
        if result:
            summary_rows.append(result)

    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(outdir / "umap_panel_summary.csv", index=False)
        print(f"\n── UMAP Panel Summary ──\n{summary.to_string(index=False)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="UMAP visualization for PAMPA analysis")
    parser.add_argument("--matrix", "-m", default="results/feature_matrix.csv")
    parser.add_argument("--outdir", "-o", default="results")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.matrix, Path(args.outdir))
