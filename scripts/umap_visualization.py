"""
04_umap_visualization.py
------------------------
Dual-track clustering workflow for CycPeptMPDB PAMPA subset.

Common ground:
  RobustScaler → PCA → cosine distance matrix (input to both tracks)

Track A — K-Medoids (structural archetypes):
  Forces data into N_KMEDOIDS clusters using cosine distance.
  Each medoid is a real molecule representing its archetype.
  Goal: identify which 3D shape archetypes are permeable winners.

Track B — UMAP + HDBSCAN (natural signal):
  UMAP projects PCA-space into 2D using cosine metric.
  HDBSCAN on 2D UMAP coordinates finds density-based clusters.
  Noise points (label=-1) are explicitly modeled — not forced into clusters.
  Goal: find permeability islands where data naturally clumps.

Visualization per panel (3 subplots, same UMAP layout):
  Plot 1 — K-Medoid cluster IDs (medoids marked ★)
  Plot 2 — HDBSCAN cluster IDs (noise in grey)
  Plot 3 — PAMPA LogPexp (the clincher)

Convergence analysis:
  Where K-Medoids and HDBSCAN agree on high-permeability regions
  = double-validated permeability islands.
  If HDBSCAN labels a molecule "noise" but K-Medoids puts it in a
  permeable archetype → the molecule has the right average shape
  but lacks the specific 3D density to cross the membrane.

Usage:
  python umap_visualization.py [--matrix results/feature_matrix.csv]
                                [--outdir results]
                                [--k 8]
"""

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from matplotlib.colors import Normalize
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import RobustScaler

warnings.filterwarnings("ignore")

# ── Optional dependencies ──────────────────────────────────────────────────────
try:
    from sklearn_extra.cluster import KMedoids
    _KMEDOIDS_AVAILABLE = True
except ImportError:
    _KMEDOIDS_AVAILABLE = False
    print("WARNING: sklearn-extra not installed — K-Medoids unavailable.")
    print("  Fix: pip install scikit-learn-extra")

try:
    import hdbscan as hdbscan_lib
    _HDBSCAN_AVAILABLE = True
except ImportError:
    _HDBSCAN_AVAILABLE = False
    print("WARNING: hdbscan not installed — HDBSCAN unavailable.")
    print("  Fix: pip install hdbscan")

# ── Constants ──────────────────────────────────────────────────────────────────
PAMPA_THRESHOLD = -6.0
RANDOM_STATE    = 42
CYCLOA_IDS      = {1, 22, 932, 981, 1822, 1862, 2356, 7188, 7353}
N_KMEDOIDS      = 8    # override with --k argument

UMAP_PARAMS = dict(
    n_neighbors  = 30,
    min_dist     = 0.15,
    n_components = 2,
    metric       = "cosine",
    random_state = RANDOM_STATE,
    low_memory   = False,
)

HDBSCAN_PARAMS = dict(
    min_cluster_size       = 50,
    min_samples            = 10,
    cluster_selection_method = "eom",
    metric                 = "euclidean",  # on 2D UMAP coordinates
)

FEATURE_PANELS = {
    "Panel_A_2D": [
        "MolWt", "MolLogP", "TPSA", "NumHAcceptors",
        "NumHDonors", "NumRotatableBonds", "RingCount",
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
    "Panel_C_normalized": [
        "MolWt", "MolLogP", "TPSA",
        "delta_3DPSA_db", "delta_psa3d_per_mw", "delta_hb", "delta_Rg",
        "psa3d_spread", "delta_NPR1",
    ],
}


# ── Enrichment helpers ─────────────────────────────────────────────────────────

def enrichment_table(labels: np.ndarray, permeable: np.ndarray,
                     noise_label: int = None) -> pd.DataFrame:
    """Permeability enrichment ratio per cluster label."""
    bg_rate = permeable.mean()
    rows = []
    for lab in sorted(np.unique(labels)):
        mask = labels == lab
        rate = permeable[mask].mean()
        rows.append({
            "cluster":     int(lab),
            "label":       "noise" if lab == noise_label else f"C{lab}",
            "n":           int(mask.sum()),
            "n_permeable": int(permeable[mask].sum()),
            "perm_rate":   round(float(rate), 3),
            "enrichment":  round(float(rate / bg_rate), 3) if bg_rate > 0 else np.nan,
        })
    return pd.DataFrame(rows).sort_values("enrichment", ascending=False)


def convergence_analysis(km_labels: np.ndarray, hdb_labels: np.ndarray,
                         permeable: np.ndarray) -> pd.DataFrame:
    """
    Cross-tabulate K-Medoid vs HDBSCAN clusters.
    Flags double-validated permeability islands where both tracks
    agree on above-background permeable enrichment (>1.2× background).
    """
    bg_rate = permeable.mean()
    rows = []
    for km_lab in np.unique(km_labels):
        km_mask = km_labels == km_lab
        km_rate = permeable[km_mask].mean()

        # Dominant HDBSCAN cluster within this K-Medoid region (ignoring noise)
        hdb_in_km = hdb_labels[km_mask]
        non_noise  = hdb_in_km[hdb_in_km != -1]
        if len(non_noise):
            dominant_hdb = int(pd.Series(non_noise).mode().iloc[0])
            hdb_mask = hdb_labels == dominant_hdb
            hdb_rate = permeable[hdb_mask].mean()
        else:
            dominant_hdb = -1
            hdb_rate     = np.nan

        # Noise analysis: molecules K-Medoids calls permeable but HDBSCAN calls noise
        noise_in_km = (hdb_in_km == -1).sum()

        both_enriched = (
            km_rate > bg_rate * 1.2
            and dominant_hdb != -1
            and not np.isnan(hdb_rate)
            and hdb_rate > bg_rate * 1.2
        )
        rows.append({
            "km_cluster":            int(km_lab),
            "km_n":                  int(km_mask.sum()),
            "km_perm_rate":          round(float(km_rate), 3),
            "km_enrichment":         round(float(km_rate / bg_rate), 3),
            "dominant_hdb_cluster":  dominant_hdb,
            "hdb_perm_rate":         round(float(hdb_rate), 3) if not np.isnan(hdb_rate) else None,
            "noise_in_km_region":    int(noise_in_km),
            "double_validated":      both_enriched,
        })
    return pd.DataFrame(rows).sort_values("km_perm_rate", ascending=False)


# ── UMAP stability check ───────────────────────────────────────────────────────

STABILITY_SEEDS   = [42, 1, 7, 99, 314]   # 5 seeds → 10 pairwise ARI comparisons
STABILITY_ARI_MIN = 0.85                   # minimum acceptable pairwise ARI

def check_umap_stability(X_red: np.ndarray, outdir: Path, panel_name: str) -> float:
    """
    Run UMAP + HDBSCAN with STABILITY_SEEDS seeds and compute pairwise ARI
    on non-noise points only.

    Rationale: HDBSCAN clusters are computed on the 2D UMAP coordinates, so
    they inherit any instability in the UMAP layout.  If the permeability
    islands shift between seeds, HDBSCAN cluster boundaries shift with them
    and any scientific claim based on those clusters becomes unreliable.
    ARI is evaluated only on non-noise points because noise assignment
    (label=-1) naturally varies across seeds — penalising that would
    artificially tank ARI for otherwise stable cluster cores.

    Returns the minimum pairwise ARI across all seed pairs.
    Prints a WARNING if any pair is below STABILITY_ARI_MIN (0.85).
    """
    if not _HDBSCAN_AVAILABLE:
        print("  [Stability] Skipped — hdbscan not installed.")
        return np.nan

    print(f"\n  [Stability] Running UMAP×{len(STABILITY_SEEDS)} seeds "
          f"(ARI threshold = {STABILITY_ARI_MIN}) ...")

    all_labels = []
    all_valid  = []   # boolean mask: non-noise for each seed run

    for seed in STABILITY_SEEDS:
        params = {**UMAP_PARAMS, "random_state": seed}
        emb    = umap.UMAP(**params).fit_transform(X_red)
        lbl    = hdbscan_lib.HDBSCAN(**HDBSCAN_PARAMS).fit_predict(emb)
        all_labels.append(lbl)
        all_valid.append(lbl != -1)

    # Pairwise ARI restricted to points that are non-noise in BOTH runs
    ari_values = []
    n = len(STABILITY_SEEDS)
    for i in range(n):
        for j in range(i + 1, n):
            shared = all_valid[i] & all_valid[j]
            if shared.sum() < 20:
                continue   # too few shared non-noise points to be meaningful
            ari = adjusted_rand_score(
                all_labels[i][shared],
                all_labels[j][shared],
            )
            ari_values.append(ari)

    if not ari_values:
        print("  [Stability] Could not compute ARI — all HDBSCAN runs returned noise only.")
        return np.nan

    min_ari  = float(np.min(ari_values))
    mean_ari = float(np.mean(ari_values))
    print(f"  [Stability] Pairwise ARI — min={min_ari:.3f}  mean={mean_ari:.3f}  "
          f"n_pairs={len(ari_values)}")

    if min_ari < STABILITY_ARI_MIN:
        print(
            f"  *** WARNING: UMAP layout unstable for {panel_name} ***\n"
            f"  Min pairwise ARI={min_ari:.3f} < {STABILITY_ARI_MIN}.\n"
            f"  HDBSCAN cluster boundaries shift between seeds — do not\n"
            f"  interpret specific cluster assignments as scientific evidence.\n"
            f"  Consider increasing n_neighbors (currently "
            f"{UMAP_PARAMS['n_neighbors']}) or min_cluster_size."
        )
    else:
        print(f"  [Stability] PASS — layout stable (min ARI {min_ari:.3f} >= {STABILITY_ARI_MIN})")

    # Save ARI matrix to CSV for reporting
    ari_path = outdir / "figures" / f"{panel_name}_umap_stability_{len(X_red)}.csv"
    ari_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    idx  = 0
    for i in range(n):
        for j in range(i + 1, n):
            shared = all_valid[i] & all_valid[j]
            if shared.sum() < 20:
                rows.append({"seed_i": STABILITY_SEEDS[i], "seed_j": STABILITY_SEEDS[j],
                             "ARI": None, "n_shared_non_noise": int(shared.sum())})
            else:
                rows.append({"seed_i": STABILITY_SEEDS[i], "seed_j": STABILITY_SEEDS[j],
                             "ARI": round(ari_values[idx], 4),
                             "n_shared_non_noise": int(shared.sum())})
                idx += 1
    pd.DataFrame(rows).to_csv(ari_path, index=False)
    print(f"  Stability ARI table saved: {ari_path}")

    return min_ari


# ── Main panel function ────────────────────────────────────────────────────────

def make_dual_track_panel(
    df: pd.DataFrame,
    panel_name: str,
    features: list,
    outdir: Path,
    n_kmedoids: int = N_KMEDOIDS,
    extra_cols: list = None,
) -> dict:
    """
    Full dual-track pipeline for one feature panel.
    Returns a metrics dict for the summary table.
    """
    available = [
        f for f in features
        if f in df.columns and df[f].notna().sum() > 50
    ]
    if len(available) < 2:
        print(f"  {panel_name}: only {len(available)} features available — skipping")
        return {}

    # extra_cols: pass-through columns for Track D (e.g. MolWt); not used as features
    base_cols = available + ["PAMPA", "permeable", "ID"]
    extra = [c for c in (extra_cols or []) if c in df.columns and c not in base_cols]
    sub = df[base_cols + extra].dropna(subset=base_cols).copy()
    print(f"\n── {panel_name} ──")
    print(f"  Features  : {available}")
    print(f"  Compounds : {len(sub):,}")

    X      = sub[available].values
    y_cont = sub["PAMPA"].values
    y_bin  = sub["permeable"].values.astype(int)

    # ── RobustScaler ─────────────────────────────────────────────────────────
    # RobustScaler (median=0, IQR=1) is more appropriate than StandardScaler
    # for cyclic peptides: the PAMPA subset contains outlier compounds with
    # extreme MW, logP, or ring counts that would inflate StandardScaler's mean
    # and standard deviation, distorting cosine distances for the majority.
    # RobustScaler is resistant to these outliers while still neutralizing
    # unit differences across heterogeneous 3D descriptors.
    X_scaled = RobustScaler().fit_transform(X)

    # ── PCA intentionally omitted ─────────────────────────────────────────────
    # Given the curated nature of the 12-feature descriptor panel, PCA was
    # omitted to preserve the physical interpretability of individual 3D
    # gatekeepers during clustering. With ≤12 features, PCA would compress
    # physically meaningful axes (e.g., ΔPSA, ΔHB, ΔRg) into abstract
    # components, making it impossible to ask "which descriptor drove this
    # cluster's permeability?" — the core scientific question.
    X_red = X_scaled

    # ── Track A: K-Medoids on PCA-cosine distance ────────────────────────────
    print(f"\n  [Track A] K-Medoids (k={n_kmedoids}, metric=cosine) ...")
    km_labels    = np.zeros(len(sub), dtype=int)
    medoid_idx   = []
    sil_km       = np.nan

    if _KMEDOIDS_AVAILABLE:
        km = KMedoids(
            n_clusters   = n_kmedoids,
            metric       = "cosine",
            method       = "alternate",   # fast for large n
            random_state = RANDOM_STATE,
        )
        km_labels  = km.fit_predict(X_red)
        medoid_idx = list(km.medoid_indices_)
        if len(np.unique(km_labels)) > 1:
            sil_km = silhouette_score(X_red, km_labels, metric="cosine")
        print(f"  K-Medoids done. Silhouette (cosine, PCA coords): "
              f"{sil_km:.4f}" if not np.isnan(sil_km) else "  K-Medoids done.")
        # Report medoid permeability
        for idx in medoid_idx:
            cid = int(sub.iloc[idx]["ID"])
            perm = "permeable" if y_bin[idx] else "impermeable"
            print(f"    Medoid cluster {km_labels[idx]}: ID={cid}  PAMPA={y_cont[idx]:.2f}  {perm}")
    else:
        print("  K-Medoids skipped (scikit-learn-extra not installed)")

    # ── UMAP stability check (before committing to main run) ─────────────────
    # Runs UMAP+HDBSCAN across 5 seeds and computes pairwise ARI on non-noise
    # points.  If min ARI < 0.85 the layout is unstable and HDBSCAN cluster
    # boundaries should not be used as scientific evidence.
    stability_ari = check_umap_stability(X_red, outdir, panel_name)

    # ── Track B: UMAP (same cosine/n_neighbors as conceptual kNN graph) ──────
    print(f"\n  [Track B] UMAP (cosine, n_neighbors={UMAP_PARAMS['n_neighbors']}) ...")
    reducer   = umap.UMAP(**UMAP_PARAMS)
    embedding = reducer.fit_transform(X_red)
    print("  UMAP done.")

    # HDBSCAN on 2D UMAP coordinates
    print(f"  [Track B] HDBSCAN (min_cluster_size={HDBSCAN_PARAMS['min_cluster_size']}) ...")
    hdb_labels = np.full(len(sub), -1, dtype=int)
    sil_hdb    = np.nan

    if _HDBSCAN_AVAILABLE:
        clusterer  = hdbscan_lib.HDBSCAN(**HDBSCAN_PARAMS)
        hdb_labels = clusterer.fit_predict(embedding)
        n_hdb      = len(set(hdb_labels) - {-1})
        n_noise    = (hdb_labels == -1).sum()
        print(f"  HDBSCAN: {n_hdb} clusters + {n_noise} noise "
              f"({100*n_noise/len(sub):.1f}%)")
        non_noise_mask = hdb_labels != -1
        if len(set(hdb_labels[non_noise_mask])) > 1:
            sil_hdb = silhouette_score(
                X_red[non_noise_mask], hdb_labels[non_noise_mask], metric="cosine"
            )
            print(f"  Silhouette (HDBSCAN non-noise, cosine, PCA coords): {sil_hdb:.4f}")
    else:
        print("  HDBSCAN skipped (hdbscan not installed)")

    # ── Enrichment tables ────────────────────────────────────────────────────
    km_enrich  = enrichment_table(km_labels,  y_bin)
    hdb_enrich = enrichment_table(hdb_labels, y_bin, noise_label=-1)
    conv_df    = convergence_analysis(km_labels, hdb_labels, y_bin)

    print("\n  K-Medoids enrichment (top 5):")
    print(km_enrich.head(5).to_string(index=False))
    print("\n  HDBSCAN enrichment (excl. noise, top 5):")
    print(hdb_enrich[hdb_enrich["cluster"] != -1].head(5).to_string(index=False))

    dv = conv_df[conv_df["double_validated"]]
    if len(dv):
        print(f"\n  Double-validated permeability islands ({len(dv)}):")
        print(dv.to_string(index=False))
    else:
        print("\n  No double-validated islands found.")

    # ── Figure: 3 or 4 subplots (4 when MolWt Track D is requested) ──────────
    has_mw = "MolWt" in sub.columns
    n_subplots = 4 if has_mw else 3
    fig_width  = 28 if has_mw else 22
    fig, axes = plt.subplots(1, n_subplots, figsize=(fig_width, 6))
    cycloA_mask = sub["ID"].isin(CYCLOA_IDS).values

    # --- Subplot 1: K-Medoids ---
    ax1 = axes[0]
    cmap_km = plt.cm.tab10
    for lab in np.unique(km_labels):
        mask     = km_labels == lab
        pct_perm = 100 * y_bin[mask].mean()
        ax1.scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=[cmap_km(int(lab) % 10)], s=8, alpha=0.5, rasterized=True,
            label=f"K{lab} n={mask.sum()} ({pct_perm:.0f}%)",
        )
    if medoid_idx:
        ax1.scatter(
            embedding[medoid_idx, 0], embedding[medoid_idx, 1],
            s=200, marker="*", c="black", zorder=10,
            edgecolors="white", linewidths=0.8, label="Medoids",
        )
    sil_str = f"sil={sil_km:.3f}" if not np.isnan(sil_km) else ""
    ax1.set_title(f"Track A — K-Medoids (k={n_kmedoids})\n{sil_str}")
    ax1.set_xlabel("UMAP 1"); ax1.set_ylabel("UMAP 2")
    ax1.legend(fontsize=6, ncol=2, loc="upper right")

    # --- Subplot 2: HDBSCAN ---
    ax2 = axes[1]
    cmap_hdb    = plt.cm.tab20
    unique_hdb  = sorted(set(hdb_labels))
    for lab in unique_hdb:
        mask     = hdb_labels == lab
        pct_perm = 100 * y_bin[mask].mean()
        if lab == -1:
            ax2.scatter(
                embedding[mask, 0], embedding[mask, 1],
                c="lightgrey", s=5, alpha=0.3, rasterized=True,
                label=f"Noise n={mask.sum()}",
            )
        else:
            ax2.scatter(
                embedding[mask, 0], embedding[mask, 1],
                c=[cmap_hdb(int(lab) % 20)], s=8, alpha=0.6, rasterized=True,
                label=f"C{lab} n={mask.sum()} ({pct_perm:.0f}%)",
            )
    sil_str2 = f"sil={sil_hdb:.3f}" if not np.isnan(sil_hdb) else ""
    n_hdb_clusters = len(unique_hdb) - (1 if -1 in unique_hdb else 0)
    ax2.set_title(f"Track B — HDBSCAN ({n_hdb_clusters} clusters)\n{sil_str2}")
    ax2.set_xlabel("UMAP 1"); ax2.set_ylabel("UMAP 2")
    ax2.legend(fontsize=6, ncol=2, loc="upper right")

    # --- Subplot 3: PAMPA LogPexp ---
    ax3 = axes[2]
    norm = Normalize(
        vmin=np.percentile(y_cont, 5),
        vmax=np.percentile(y_cont, 95),
    )
    sc = ax3.scatter(
        embedding[:, 0], embedding[:, 1],
        c=y_cont, cmap=plt.cm.RdYlGn, norm=norm,
        s=8, alpha=0.5, rasterized=True,
    )
    plt.colorbar(sc, ax=ax3, label="PAMPA LogPexp (log cm/s)")
    if cycloA_mask.sum() > 0:
        ax3.scatter(
            embedding[cycloA_mask, 0], embedding[cycloA_mask, 1],
            s=150, marker="*", c="black", zorder=10,
            edgecolors="white", linewidths=0.8,
            label=f"CycloA (n={cycloA_mask.sum()})",
        )
        ax3.legend(fontsize=9)
    ax3.set_title("The Clincher — PAMPA LogPexp\n(validate cluster chemistry)")
    ax3.set_xlabel("UMAP 1"); ax3.set_ylabel("UMAP 2")

    # --- Subplot 4: Track D — Molecular Weight (optional) ---
    if has_mw:
        ax4 = axes[3]
        mw_vals = sub["MolWt"].values
        norm_mw = Normalize(
            vmin=np.percentile(mw_vals, 5),
            vmax=np.percentile(mw_vals, 95),
        )
        sc_mw = ax4.scatter(
            embedding[:, 0], embedding[:, 1],
            c=mw_vals, cmap=plt.cm.plasma, norm=norm_mw,
            s=8, alpha=0.5, rasterized=True,
        )
        plt.colorbar(sc_mw, ax=ax4, label="Molecular Weight (Da)")
        if cycloA_mask.sum() > 0:
            ax4.scatter(
                embedding[cycloA_mask, 0], embedding[cycloA_mask, 1],
                s=150, marker="*", c="cyan", zorder=10,
                edgecolors="white", linewidths=0.8,
                label=f"CycloA (n={cycloA_mask.sum()})",
            )
            ax4.legend(fontsize=9)
        # Overlay permeability boundary: mark permeable points with a thin ring
        perm_mask = y_bin.astype(bool)
        ax4.scatter(
            embedding[perm_mask, 0], embedding[perm_mask, 1],
            s=18, facecolors="none", edgecolors="limegreen",
            linewidths=0.4, alpha=0.35, zorder=5,
            label=f"Permeable (n={perm_mask.sum()})",
        )
        ax4.legend(fontsize=7, loc="upper right")
        # Annotate median MW for permeable vs impermeable
        med_perm   = np.median(mw_vals[perm_mask])
        med_imperm = np.median(mw_vals[~perm_mask])
        ax4.set_title(
            f"Track D — Molecular Weight\n"
            f"Med permeable={med_perm:.0f} Da  |  impermeable={med_imperm:.0f} Da"
        )
        ax4.set_xlabel("UMAP 1"); ax4.set_ylabel("UMAP 2")

    n_label = f"{len(sub):,}".replace(",", "")
    fig.suptitle(
        f"{panel_name}  (n={len(sub):,})  |  "
        f"Dual-Track: K-Medoids (archetypes) + HDBSCAN (natural signal)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig_path = outdir / "figures" / f"{panel_name}_umap_{n_label}.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Figure saved: {fig_path}")

    # ── Save per-compound embedding + cluster labels ─────────────────────────
    save_cols = ["ID", "PAMPA", "permeable"] + (["MolWt"] if has_mw else [])
    sub_out = sub[save_cols].copy()
    sub_out["embedding_x"]      = embedding[:, 0]
    sub_out["embedding_y"]      = embedding[:, 1]
    sub_out["kmedoids_cluster"] = km_labels
    sub_out["hdbscan_cluster"]  = hdb_labels
    sub_out.to_csv(outdir / f"{panel_name}_embedding_{n_label}.csv", index=False)

    fig_dir = outdir / "figures"
    km_enrich.to_csv(fig_dir / f"{panel_name}_kmedoids_enrichment_{n_label}.csv", index=False)
    hdb_enrich.to_csv(fig_dir / f"{panel_name}_hdbscan_enrichment_{n_label}.csv", index=False)
    conv_df.to_csv(fig_dir / f"{panel_name}_convergence_{n_label}.csv", index=False)

    return {
        "panel":                    panel_name,
        "n_compounds":              len(sub),
        "n_features":               len(available),
        "n_kmedoids":               n_kmedoids if _KMEDOIDS_AVAILABLE else 0,
        "sil_kmedoids":             round(float(sil_km), 4) if not np.isnan(sil_km) else None,
        "n_hdbscan_clusters":       n_hdb_clusters,
        "n_hdbscan_noise":          int((hdb_labels == -1).sum()),
        "sil_hdbscan":              round(float(sil_hdb), 4) if not np.isnan(sil_hdb) else None,
        "double_validated_islands": int(conv_df["double_validated"].sum()),
        "cycloA_found":             int(cycloA_mask.sum()),
        "umap_stability_min_ari":   round(stability_ari, 4) if not np.isnan(stability_ari) else None,
        "umap_stability_pass":      bool(stability_ari >= STABILITY_ARI_MIN) if not np.isnan(stability_ari) else None,
    }


# ── Entry point ────────────────────────────────────────────────────────────────

def run(matrix_csv: str, outdir: Path, n_kmedoids: int,
        sources: list = None, panels: list = None) -> None:
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(matrix_csv, low_memory=False)
    df = df[df["PAMPA"].notna()].copy()
    df["permeable"] = (df["PAMPA"] >= PAMPA_THRESHOLD).astype(int)

    if sources:
        df = df[df["Source"].isin(sources)].copy()
        print(f"Source filter: {sources}")

    print(f"Loaded {len(df):,} compounds with PAMPA values")
    print(f"Permeable: {df['permeable'].sum():,} ({100*df['permeable'].mean():.1f}%)")

    active_panels = {k: v for k, v in FEATURE_PANELS.items()
                     if panels is None or k in panels}

    summary_rows = []
    for panel_name, features in active_panels.items():
        result = make_dual_track_panel(
            df, panel_name, features, outdir, n_kmedoids,
            extra_cols=["MolWt"],
        )
        if result:
            summary_rows.append(result)

    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        n_total = summary["n_compounds"].iloc[0] if len(summary) else "all"
        suffix  = f"_{n_total}" if sources else ""
        summary.to_csv(outdir / f"umap_panel_summary{suffix}.csv", index=False)
        print(f"\n── Summary ──\n{summary.to_string(index=False)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dual-track UMAP visualization")
    parser.add_argument("--matrix",  "-m", default="results/feature_matrix.csv")
    parser.add_argument("--outdir",  "-o", default="results")
    parser.add_argument("--k",       "-k", type=int, default=N_KMEDOIDS,
                        help=f"Number of K-Medoid archetypes (default: {N_KMEDOIDS})")
    parser.add_argument("--sources", "-s", nargs="+", default=None,
                        help="Filter to these Source values (e.g. 2016_Furukawa 2013_CHUGAI)")
    parser.add_argument("--panels",  "-p", nargs="+", default=None,
                        help="Run only these panels (e.g. Panel_C_combined)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.matrix, Path(args.outdir), args.k,
        sources=args.sources, panels=args.panels)
