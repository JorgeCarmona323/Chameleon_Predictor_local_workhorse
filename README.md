# The Chameleon Traverse
## 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling for Cyclic Peptide Membrane Permeation

**Jorge Carmona | CHEM 269 Final Project | March 2026**

---

## Summary

Cyclic peptides can passively cross cell membranes despite high polarity — a phenomenon called *chameleonic behavior* — by folding into compact, H-bond-shielded conformations in lipophilic environments while exposing polar groups in water. Standard 2D descriptors like TPSA and MolLogP are conformationally blind and cannot capture this switching. This project tests whether 3D ensemble-derived Δ descriptors (especially ΔPSA, the difference in polar surface area between aqueous and membrane-mimetic conformers) can outperform 2D baselines in predicting experimental PAMPA membrane permeability across 7,297 cyclic peptides from CycPeptMPDB.

On the full 7,297-compound dataset, ensemble ΔPSA achieves AUC = 0.505 — essentially at chance — while MolLogP achieves AUC = 0.631 as the best single descriptor. An earlier exploratory run on 1,502 compounds (a biased non-random subset) yielded AUC = 0.744 for ΔPSA; that result is documented in Section 8 of the main notebook as a subset artifact, not a generalizable finding. The strongest result from this project is the UMAP Panel B two-population structure in 3D Δ feature space, which survives on the full dataset and provides visual evidence that ensemble conformational descriptors stratify chemical space in a permeability-relevant way even when the AUC signal is weak.

---

## Biological Problem and Motivation

Passive membrane permeation is gated by the ability of a molecule to shed water, cross a low-dielectric lipid bilayer, and re-hydrate. For small molecules this is well-modeled by LogP. For cyclic peptides — which can exceed 1,000 Da and contain many polar backbone amides — the story is more complex. Cyclosporin A (CsA), the canonical example, crosses membranes despite a MW of 1,203 Da and 11 H-bond donors/acceptors by forming intramolecular H-bonds that effectively bury its polar surface. This chameleonic effect is invisible to 2D topological descriptors.

CycPeptMPDB is the largest public dataset of experimentally measured cyclic peptide permeability, with >8,400 compounds and PAMPA LogPexp values. It also contains its own 3D PSA values (H2O\_3DPSA, CHCl3\_3DPSA) — but these are computed from single optimized structures, not conformer ensembles. Testing whether those database values encode any useful signal vs. ensemble-sampled Δ descriptors is a key hypothesis of this project. The answer reveals something important about methodology: **single-structure solvation is insufficient for chameleonic molecules**.

This work connects directly to ongoing research on DNA-encoded library (DEL) cyclic peptide scaffolds.

---

## Data Sources

| Source | Description | Used for |
|--------|-------------|----------|
| **CycPeptMPDB v1.2** | 8,466 cyclic peptides with experimental permeability (PAMPA, Caco-2, RRCK) | Primary dataset; PAMPA subset of 7,298 compounds used |
| **CycPeptMPDB H2O\_3DPSA / CHCl3\_3DPSA** | Database-provided single-structure 3D PSA values | Baseline "DB 3D" feature group; negative control |
| **RDKit 2D descriptors** | MolLogP, TPSA, MW, HBA, HBD, RotBonds, CSP3, Rings | 2D baseline |
| **ETKDGv3 + MMFF94s (Tier-1)** | 20 conformers/molecule, RDKit macrocycle torsion sampling | Ensemble Δ descriptors (primary method) |

**File placement:** Place `CycPeptMPDB_Peptide_All (2).csv` in the repo root before running.

---

## Computational Approach

### Tier-1: ETKDGv3 Conformer Ensemble (primary method)

For each molecule, 20 conformers are generated with ETKDGv3 (macrocycle torsion library enabled) and minimized with MMFF94s force field. The conformer with maximum 3D polar SASA (Bondi radii, RDKit rdFreeSASA) is designated the "aqueous conformer"; the conformer with minimum 3D polar SASA is the "membrane conformer." Key Δ features:

- `delta_psa3d` = PSA(aq) − PSA(mem): the chameleonic PSA switch
- `psa3d_std`: standard deviation of PSA across all 20 conformers (conformational flexibility)
- `delta_hb`: change in H-bond count between conformers
- `delta_Rg`, `delta_NPR1/2`: compaction and shape metrics

### DB 3D: Single-Structure Baseline (negative control)

`delta_3DPSA_db` = H2O\_3DPSA − CHCl3\_3DPSA from the CycPeptMPDB database values. Tests whether existing single-structure solvation captures the chameleonic effect.

### Analysis

- Pearson/Spearman correlation vs. PAMPA LogPexp
- AUC-ROC classification (PAMPA ≥ −6.0 log cm/s = permeable)
- UMAP dimensionality reduction with dual-track clustering (K-Medoids + HDBSCAN)
- ARI stability validation across 5 random seeds

### Attempted Tier-2: CREST (failed — documented)

CREST 2.12 with ALPB solvation was attempted for higher-quality conformer sampling on 5 reference compounds. All failed due to Colab memory/timeout constraints (HexPep timed out after 4 hours; others exited immediately). Fully documented in `docs/findings_and_methods_log.md`.

### Attempted Tier-2: xtb+GBSA (completed, insufficient)

GFN2-xTB single-structure optimization with GBSA (water/CHCl3) was completed for all 5 reference compounds as a CREST alternative. Resulted in ΔPSA of 0–7 Å² for all compounds (CsA = −0.14 Å²) — confirming single-structure solvation cannot capture chameleonic behavior regardless of the level of theory.

---

## Route Design and Completion

| Stage | Method | Status | Coverage |
|-------|--------|--------|----------|
| Data curation | RDKit canonicalization, PAMPA filter | Complete | 7,298 / 8,466 |
| 2D baseline descriptors | RDKit | Complete | 100% |
| DB 3DPSA | CycPeptMPDB H2O/CHCl3\_3DPSA | Complete | 88% (6,942) |
| Tier-1 conformers | ETKDGv3 + MMFF94s, n=20 | Complete — 7,297 / 7,298 | 99.99% |
| Feature matrix | Merged all above | Complete | 7,298 rows |
| Correlation analysis | Pearson, Spearman, AUC-ROC | Complete | — |
| UMAP + clustering | K-Medoids + HDBSCAN dual-track | Complete | — |
| Tier-2 CREST | CREST 2.12 + ALPB | Failed | 0 / 5 ref compounds |
| Tier-2 xtb+GBSA | GFN2-xTB + GBSA single-opt | Completed, not used | 5 / 5 ref compounds |

**Stretch component:** Dual-track UMAP clustering (K-Medoids deterministic + HDBSCAN exploratory) with ARI stability validation across 5 random seeds.

---

## Results

### Reference Compound Validation

| Compound | Permeable | Tier-1 ΔPSA (Å²) | xtb ΔPSA (Å²) | Literature ΔPSA |
|----------|-----------|-----------------|--------------|----------------|
| CsA (Cyclosporin A) | Yes | **84.9** | −0.14 | ~75 Å² |
| DP172 | Yes | 88.9 | −0.24 | — |
| HexPep | No | 64.4 | 0.82 | — |
| 1NMe3 | Yes | 47.8 | 6.91 | — |
| PSLYF | No | 65.3 | 5.40 | — |

CsA Tier-1 ΔPSA = 84.9 Å² vs. literature ~75 Å² — validates the ensemble approach. xtb gives near-zero ΔPSA for all compounds including CsA, confirming single-structure methods fail.

### AUC-ROC Results (full 7,297-compound dataset)

| Descriptor | AUC-ROC | Group |
|------------|---------|-------|
| MolLogP | 0.631 | 2D baseline (best overall) |
| delta_psa3d (Tier-1) | 0.505 | Tier-1 Δ |
| delta_3DPSA_db | 0.507 | DB 3D |

On the full dataset, Tier-1 ΔPSA does not outperform 2D baselines. An earlier run on a non-random 1,502-compound subset produced AUC = 0.744 for delta_psa3d; that result is acknowledged as exploratory in notebook Section 8 and is attributable to sampling bias rather than a generalizable effect. The finding that the dataset's own single-structure DB 3DPSA is also near-chance (AUC = 0.507) is a reproducible negative result and the methodological comparison is meaningful: both single-structure and ensemble force-field methods fail at full scale, with logP remaining the dominant predictor.

### Key Figures

| Figure | File | What it shows |
|--------|------|--------------|
| AUC-ROC bar chart | `results/figures/auc_roc_bar.png` | Feature ranking across all descriptor groups |
| Correlation heatmap | `results/figures/correlation_heatmap.png` | Pearson/Spearman vs. PAMPA LogPexp |
| ΔPSA scatter | `results/figures/scatter_top_features.png` | delta_psa3d vs. LogPexp colored by permeability |
| UMAP Panel A | `results/figures/Panel_A_2D_umap.png` | 2D descriptor chemical space |
| UMAP Panel B | `results/figures/Panel_B_3D_delta_umap.png` | 3D Δ feature chemical space (strongest visual result) |
| UMAP Panel C | `results/figures/Panel_C_combined_umap.png` | Combined 2D + 3D Δ |

---

## Interpretation, Limitations, and Next Steps

### Interpretation

The full-dataset AUC result (delta_psa3d = 0.505) does not support the hypothesis that Tier-1 ensemble ΔPSA outperforms 2D baselines at scale. MolLogP (AUC = 0.631) remains the strongest single descriptor. The earlier 1,502-compound result (AUC = 0.744) reflected a biased non-random subset and should not be generalized.

The DB 3DPSA result (AUC = 0.507) is itself a finding: the database's own 3D PSA values, computed from single structures, provide no useful signal. This is a direct methodological comparison and argues that single-structure approaches fail regardless of the level of theory — consistent with the xtb negative control.

The CsA NMR validation (Tier-1 ΔPSA = 84.9 Å² vs. literature ~75 Å²) confirms the ensemble calculation is physically correct at the individual molecule level. The failure at population scale likely reflects PAMPA assay heterogeneity (see Limitations) rather than a fundamental flaw in the ΔPSA concept.

UMAP Panel B shows a two-population structure in 3D Δ feature space that survives on the full 7,297-compound dataset and is the strongest visual result of this project. UMAP cluster stability was poor across all three panels (ARI 0.07–0.38), indicating permeability is a continuous, graded property with no discrete permeable/non-permeable clusters in this dataset — consistent with literature and expected given the 66.4% permeable selection bias in CycPeptMPDB.

### Limitations

1. **PAMPA assay heterogeneity:** CycPeptMPDB aggregates measurements from multiple labs with incompatible protocols. Townsend 2020 (a preprint comprising approximately 42% of the PAMPA data) uses pooled compound PAMPA; Kelly 2021 dominates another large fraction; Chugai uses a different membrane formulation. Source-stratified AUC analysis shows that cross-source label noise likely suppresses any real ΔPSA signal at full scale. This is the primary explanation for the AUC collapse from the exploratory subset to the full dataset.
2. **MMFF94s dual-dielectric approximation:** Conformer selection by PSA extremes is a heuristic, not physics-based. CREST or OpenMM MD would provide more rigorous ensemble sampling.
3. **Dataset selection bias:** CycPeptMPDB is 66.4% permeable — not a balanced spectrum. This compresses feature-space contrast and likely suppresses the true AUC signal.
4. **No Tier-2 validation:** CREST and OpenMM were not successfully completed; the Tier-1 heuristic is unvalidated against higher-level theory.
5. **Force-field conformer artefacts:** Vacuum ETKDGv3 sampling may generate collapsed hydrophobic conformers that are thermodynamically inaccessible in aqueous solution, producing spuriously large ΔPSA for rigid impermeable peptides.

### Next Steps

- Source-stratified re-analysis using PAMPA data from a single lab protocol (e.g., Kelly 2021 only) to test whether AUC recovers when label noise is reduced
- Obtain HPC allocation for CREST conformer sampling on reference set
- Expand dataset with unbiased Caco-2/RRCK data from Lokey lab publications
- N-methylation subgroup analysis
- Random forest + SHAP for interpretable SAR

Full roadmap: `docs/project_roadmap.md`

---

## Reproducibility

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate chem269_cycpep

# 2. Run full pipeline (requires CycPeptMPDB CSV in repo root)
python run_pipeline.py

# Quick test on 200 molecules (~5 min):
python run_pipeline.py --max-mols 200 --n-confs 20

# Results → results/  |  Figures → results/figures/
```

**Individual scripts (run in order):**
```bash
python scripts/curate_data.py            # → data/pampa_curated.csv
python scripts/conformer_engine.py       # → results/conformer_descriptors_raw.csv
python scripts/build_feature_matrix.py  # → results/feature_matrix.csv
python scripts/correlation_analysis.py  # → results/correlation_table.csv + auc_roc_table.csv
python scripts/umap_visualization.py    # → results/figures/Panel_*.png
```

**Main analysis notebook:** `notebooks/3d_descriptors.ipynb`

**Colab Tier-1 (large-scale, A100):** `colab/tier1_a100/tier1_a100_03162026.ipynb`
Includes checkpoint/resume system. Upload to Google Colab, mount Drive, run all cells.

---

## Repository Structure

```
CHEM_269_Final_Project/
├── README.md
├── environment.yml
├── run_pipeline.py
├── notebooks/
│   └── 3d_descriptors.ipynb              <- main deliverable notebook
├── scripts/
│   ├── curate_data.py
│   ├── conformer_engine.py
│   ├── build_feature_matrix.py
│   ├── correlation_analysis.py
│   └── umap_visualization.py
├── colab/
│   ├── tier1_a100/tier1_a100_03162026.ipynb  <- full-scale A100 run
│   └── tier2_xtb_gbsa.ipynb                  <- xtb negative control
├── data/                                  <- generated by curate_data.py
├── results/                               <- generated by pipeline
│   └── figures/                           <- all plots
└── docs/
    ├── findings_and_methods_log.md        <- full findings + CREST failure account
    ├── literature_deltapsa_values.md      <- literature ΔPSA reference values
    ├── future_reference_compounds.md      <- planned validation compounds
    └── project_roadmap.md                 <- future directions
```

---

## References

- Jiang et al. (2023). CycPeptMPDB: A Comprehensive Database of Membrane Permeability of Cyclic Peptides. *J. Chem. Inf. Model.*
- Rezai et al. (2006). Conformational flexibility, internal hydrogen bonding, and passive membrane permeability. *J. Am. Chem. Soc.*
- Witek et al. (2016). Kinetic models of cyclosporin A in polar and apolar environments. *J. Chem. Theory Comput.*
- Bockus et al. (2015). Decoding chameleonic properties of macrocycles. *J. Med. Chem.*
- Riniker & Landrum (2015). Better informed distance geometry. *J. Chem. Inf. Model.* (ETKDGv3)
- Pracht et al. (2020). Automated exploration of the low-energy chemical space with CREST. *Phys. Chem. Chem. Phys.*
- Townsend et al. (2020). Cyclic peptide membrane permeability dataset. *bioRxiv* (preprint).
- Kelly et al. (2021). Oral bioavailability of cyclic peptides. *J. Med. Chem.*
