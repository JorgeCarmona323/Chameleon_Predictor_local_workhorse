# 2026-03-18 — Panel C (1,566 compounds) + Track D Molecular Weight

## Question

The 1.5k Panel B looked cleaner than the 7k Panel B. Does Panel C on the same 1.5k subset sharpen further? And does the permeable cluster in Panel C correspond to higher-MW compounds — which would provide evidence that (1) per-residue ΔPSA normalization is needed, and (2) the current AUC signal is driven by large chameleonic macrolides, not small hexapeptides?

## What we did

### Pipeline change: `--sources`, `--panels`, Track D

Added three features to `scripts/umap_visualization.py`:

- `--sources` CLI flag — filter `feature_matrix.csv` to specific `Source` values before running UMAP. Allows source-stratified runs without a separate data file.
- `--panels` CLI flag — run only named panels (avoids re-running all three when iterating).
- **Track D subplot** — when `MolWt` is in the feature set (Panel C includes it), a 4th subplot is automatically added to the figure. Same UMAP embedding as Tracks A/B/C, colored by MW using the `plasma` colormap. Permeable compounds (PAMPA ≥ −6.0) are overlaid as thin limegreen rings. Median MW for permeable vs. impermeable is annotated in the title.
- Filenames now include `_n` suffix (e.g. `_1566`) so subset runs don't overwrite full-dataset outputs.

Run command:

```bash
python scripts/umap_visualization.py \
  --matrix results/feature_matrix.csv \
  --outdir results \
  --sources 2016_Furukawa 2013_CHUGAI \
  --panels Panel_C_combined \
  --k 8
```

### Subset

| Filter | n | Permeable (PAMPA ≥ −6) | Perm rate |
|--------|---|------------------------|-----------|
| 2016_Furukawa + 2013_CHUGAI | 1,566 | 1,187 | 75.8% |

MW range: 607–1,778 Da. Monomer length: 6–15 residues.

## Results

### Stability

Min pairwise ARI = **0.995** across 5 seeds (threshold 0.85). Extremely stable — one dominant structure in combined 2D+3D feature space.

### HDBSCAN clusters

| Cluster | n | Permeable | Perm rate | Enrichment |
|---------|---|-----------|-----------|------------|
| C0 | 883 | 775 | **87.8%** | 1.16× |
| C1 | 651 | 386 | 59.3% | 0.78× |
| Noise | 32 | — | — | — |

### K-Medoids enrichment (top 3)

| Cluster | n | Perm rate |
|---------|---|-----------|
| K3 | 168 | **97.6%** |
| K2 | 210 | 87.6% |
| K5 | 325 | 84.0% |

### Track D — Molecular Weight

| Population | Median MW |
|------------|-----------|
| Permeable (PAMPA ≥ −6) | **1,180 Da** |
| Impermeable | 820 Da |
| Ratio | 1.44× |

### Figure

![Panel C 1566 + Track D](../figures/Panel_C_combined_umap_1566.png)

The limegreen permeable rings in Track D cluster almost entirely in the high-MW (plasma hot) region. The impermeable-enriched C1 cluster occupies the low-MW (purple/blue) region.

## Interpretation

**The permeable cluster is large macrolides.** The 1.44× MW gap between permeable and impermeable compounds on the same UMAP is direct evidence that the AUC = 0.744 signal in the 1.5k subset is driven by CsA-class macrolides (≥900 Da, ≥9 residues) undergoing genuine chameleonic switching — not by hexapeptides achieving permeability through N-methylation or lipophilicity alone.

**Hexapeptides may be permeable but not chameleonic.** The C1 cluster (59.3% permeable, lower MW) is not impermeable — it just isn't enriched. Small cyclic peptides from the Lokey group tradition likely achieve PAMPA permeability through reduced HBD count (N-methylation) rather than ΔPSA switching. Absolute ΔPSA is the wrong descriptor for them because they lack the conformational freedom to generate a large switch.

**Absolute ΔPSA is a size proxy in the current model.** A 15-residue peptide that barely switches produces a larger raw ΔPSA (Å²) than a 9-residue peptide that fully buries its polar surface. The Yu et al. 2026 normalization (ΔPSA/SASA_total) removes this confound.

## Next experiment

1. Add `delta_psa3d_per_mw`, `delta_psa3d_per_sasa`, `delta_psa3d_per_residue` to `conformer_engine.py`
2. Re-run AUC on the 1,566-compound subset with normalized descriptors
3. Apply ≥9-residue filter (Monomer_Length ≥ 9) and compare AUC with/without it — if the normalized descriptor recovers or improves on 0.744 for large compounds while failing on small ones, that confirms chameleonic behavior as a size-gated phenomenon
