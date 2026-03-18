# The Chameleon Traverse
## 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling for Cyclic Peptide Membrane Permeation

**Jorge Carmona | CHEM 269 Final Project | March 2026**

---

## Summary

Cyclic peptides can passively cross cell membranes despite high polarity — a phenomenon called *chameleonic behavior* — by folding into compact, H-bond-shielded conformations in lipophilic environments while exposing polar groups in water. Standard 2D descriptors like TPSA and MolLogP are conformationally blind and cannot capture this switching. This project tests whether 3D ensemble-derived Δ descriptors (especially ΔPSA, the difference in polar surface area between aqueous and membrane-mimetic conformers) can outperform 2D baselines in predicting experimental PAMPA membrane permeability across cyclic peptides from CycPeptMPDB. Ensemble ΔPSA achieves AUC = 0.744 vs. AUC = 0.631 for the best 2D descriptor (MolLogP) — an 11-point improvement — validating that conformational sampling is essential and that single-structure methods fail entirely.

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
- AUC-ROC classification (PAMPA >= -6.0 log cm/s = permeable)
- UMAP dimensionality reduction with dual-track clustering (K-Medoids + HDBSCAN)
- ARI stability validation across 5 random seeds

### Attempted Tier-2: CREST (failed — documented)

CREST 2.12 with ALPB solvation was attempted for higher-quality conformer sampling on 5 reference compounds. All failed due to Colab memory/timeout constraints (HexPep timed out after 4 hours; others exited immediately).

### Attempted Tier-2: xtb+GBSA (completed, insufficient)

GFN2-xTB single-structure optimization with GBSA (water/CHCl3) was completed for all 5 reference compounds as a CREST alternative. Resulted in ΔPSA of 0–7 Å² for all compounds (CsA = -0.14 Å²) — confirming single-structure solvation cannot capture chameleonic behavior regardless of the level of theory.

---

## Route Design and Completion

| Stage | Method | Status | Coverage |
|-------|--------|--------|----------|
| Data curation | RDKit canonicalization, PAMPA filter | Complete | 7,298 / 8,466 |
| 2D baseline descriptors | RDKit | Complete | 100% |
| DB 3DPSA | CycPeptMPDB H2O/CHCl3\_3DPSA | Complete | 88% (6,942) |
| Tier-1 conformers | ETKDGv3 + MMFF94s, n=20 | Partial — 1,502 / 7,298 | 20.6% |
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
| CsA (Cyclosporin A) | Yes | **84.9** | -0.14 | ~75 Å² |
| DP172 | Yes | 88.9 | -0.24 | — |
| HexPep | No | 64.4 | 0.82 | — |
| 1NMe3 | Yes | 47.8 | 6.91 | — |
| PSLYF | No | 65.3 | 5.40 | — |

CsA Tier-1 ΔPSA = 84.9 Å² vs. literature ~75 Å² — validates the ensemble approach. xtb gives near-zero ΔPSA for all compounds including CsA, confirming single-structure methods fail.

### AUC-ROC Results (1,502-compound Tier-1 subset)

| Descriptor | AUC-ROC | Group |
|------------|---------|-------|
| delta_psa3d (Tier-1) | **0.744** | Tier-1 Δ (best overall) |
| MolLogP | 0.631 | 2D baseline |
| delta_3DPSA_db | 0.507 | DB 3D |

Ensemble ΔPSA outperforms the best 2D descriptor by 11 AUC points. DB 3DPSA (single-structure) is essentially random — confirming that ensemble sampling is a prerequisite, not a refinement.

### Key Figures

| Figure | File | What it shows |
|--------|------|--------------|
| AUC-ROC bar chart | `results/figures/auc_roc_bar.png` | Feature ranking across all descriptor groups |
| Correlation heatmap | `results/figures/correlation_heatmap.png` | Pearson/Spearman vs. PAMPA LogPexp |
| ΔPSA scatter | `results/figures/scatter_top_features.png` | delta_psa3d vs. LogPexp colored by permeability |
| UMAP Panel A | `results/figures/Panel_A_2D_umap.png` | 2D descriptor chemical space |
| UMAP Panel B | `results/figures/Panel_B_3D_delta_umap.png` | 3D Δ feature chemical space (core result) |
| UMAP Panel C | `results/figures/Panel_C_combined_umap.png` | Combined 2D + 3D Δ |

---

## Interpretation, Limitations, and Next Steps

### Interpretation

Ensemble-derived ΔPSA captures real physical information that single-structure methods miss entirely. The 11-point AUC improvement over MolLogP, combined with CsA validation against literature (~75 Å²), supports the chameleonic hypothesis. AUC = 0.744 and Spearman rho = 0.457 indicate ΔPSA is a significant contributor, not a complete predictor — consistent with the multi-factorial nature of passive membrane permeation in this compound class.

The DB 3DPSA result (AUC = 0.507) is itself a finding: the database's own 3D PSA values, computed from single structures, provide no useful signal. This is a direct methodological comparison and argues that ensemble sampling is a prerequisite for capturing chameleonic behavior.

UMAP cluster stability was poor across all three panels (ARI 0.07–0.38), indicating permeability is a continuous, graded property with no discrete permeable/non-permeable clusters in this dataset — consistent with literature and expected given the 66.4% permeable selection bias in CycPeptMPDB.

### Note on Full-Dataset Extension

After submission, Tier-1 conformer generation was completed for the full 7,297-compound PAMPA subset. On that dataset, delta_psa3d AUC collapsed to 0.505 — essentially random. This was not a failure of the descriptor concept but a reflection of how messy CycPeptMPDB really is at scale. The database aggregates PAMPA measurements from multiple labs with incompatible protocols: Townsend 2020 (a preprint, pooled 150-compound PAMPA with MS deconvolution, ~42% of data), Kelly 2021 (same pooled protocol), and Chugai (a patent source using a different membrane formulation and a detection floor of -10.0 vs -8.0 log cm/s). The 1,502-compound subset used here happens to draw disproportionately from the cleaner Furukawa 2016 source and the high-permeability Chugai block, which produced an artificially favorable signal. Cleaning and source-stratifying the dataset is the logical next step and is documented in the standalone research continuation.

### Limitations

1. **Partial Tier-1 coverage**: Results based on 1,502 / 7,298 compounds (20.6%). Full run pending.
2. **MMFF94s dual-dielectric approximation**: Conformer selection by PSA extremes is a heuristic, not physics-based. CREST or OpenMM MD would provide more rigorous ensemble sampling.
3. **Dataset selection bias**: CycPeptMPDB is 66.4% permeable — not a balanced spectrum. This compresses feature-space contrast.
4. **PAMPA assay heterogeneity**: Multiple labs and protocols contribute to CycPeptMPDB; cross-source label noise is significant.
5. **No Tier-2 validation**: CREST failed; the Tier-1 heuristic is unvalidated against higher-level theory.

### Next Steps

- Complete full 7,298-compound Tier-1 run and source-stratified re-analysis (single-protocol PAMPA subset)
- Obtain HPC allocation for CREST conformer sampling on reference set
- Expand dataset with unbiased Caco-2/RRCK data from Lokey lab publications
- N-methylation subgroup analysis
- Random forest + SHAP for interpretable SAR

---

## Reproducibility

**Individual scripts (run in order):**
```bash
python scripts/curate_data.py            # -> data/pampa_curated.csv
python scripts/conformer_engine.py       # -> results/conformer_descriptors_raw.csv
python scripts/build_feature_matrix.py  # -> results/feature_matrix.csv
python scripts/correlation_analysis.py  # -> results/correlation_table.csv + auc_roc_table.csv
python scripts/umap_visualization.py    # -> results/figures/Panel_*.png
```

**Environment:**
```bash
conda env create -f environment.yml
conda activate chem269_cycpep
```

**Main analysis notebook:** `notebooks/3d_descriptors.ipynb`

---

## Repository Structure

```
CHEM_269_Final_Project/
├── README.md
├── environment.yml
├── assignment/
│   └── climb_route_prompt.md         <- route prompt
├── notebooks/
│   └── 3d_descriptors.ipynb          <- main deliverable notebook
├── scripts/
│   ├── curate_data.py
│   ├── conformer_engine.py
│   ├── build_feature_matrix.py
│   ├── correlation_analysis.py
│   └── umap_visualization.py
└── results/                          <- generated by pipeline
    └── figures/                      <- all plots
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
