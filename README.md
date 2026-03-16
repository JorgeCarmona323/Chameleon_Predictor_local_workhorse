# CHEM 269 Final Project
## 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling to Decode Cyclic Peptide Membrane Permeation

**Jorge Carmona | March 2026 | CHEM 269**

---

## Scientific Question

Cyclic peptides exhibit **chameleonic behavior** — they can adopt different conformations in aqueous (ε=78) vs. membrane-mimetic (ε=4) environments by forming intramolecular H-bonds that shield polar groups. This enables passive membrane permeation even in molecules that violate Lipinski's rules. However, 2D descriptors (TPSA, MolLogP) are insensitive to this conformational switching.

**Can 3D Δ features computed from a conformer ensemble quantify chameleonic potential and correlate with experimental PAMPA LogPexp across ~8,000 compounds in CycPeptMPDB?**

---

## Quickstart (single command)

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate chem269_cycpep

# 2. Run full pipeline
python run_pipeline.py

# Results → results/  |  Figures → results/figures/
```

For a quick test on 200 molecules (~5 min):
```bash
python run_pipeline.py --max-mols 200 --n-confs 20
```

To skip conformer generation and use database 3DPSA values only:
```bash
python run_pipeline.py --skip-conformers
```

---

## Dataset

**CycPeptMPDB v1.2** — Publicly available at https://cycpeptmpdb.com/

- 8,466 cyclic peptides with experimental permeability measurements
- Primary target: PAMPA subset (~7,298 compounds with LogPexp)
- Permeability threshold: PAMPA LogPexp ≥ −6.0 log cm/s → "permeable"
- Database already contains `H2O_3DPSA` and `CHCl3_3DPSA` for ~7,451 compounds — used as baseline and Tier-2 cross-check

Place the CSV file at:
```
CycPeptMPDB_Peptide_All (2).csv   ← root of this repo
```

---

## Pipeline Architecture

```
CycPeptMPDB CSV (8,466 compounds)
         │
         ▼
┌─────────────────────────────────────────────────────┐
│  scripts/curate_data.py                             │
│  → SMILES canonicalization (RDKit LargestFragment)  │
│  → PAMPA subset filter (~7,298 compounds)           │
│  → Reference set curation (CycloA + analogs)        │
│  → Output: data/pampa_curated.csv                   │
└─────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────┐
│  scripts/conformer_engine.py  (Tier-1)              │
│  → ETKDGv3 (macrocycle torsions enabled)            │
│  → 50 conformers/molecule, MMFF94s minimization     │
│  → Aqueous conformer = max-PSA conformer            │
│  → Membrane conformer = min-PSA conformer           │
│  → Δ descriptors: ΔPSA, ΔHB, ΔRg, ΔNPR, PSA-spread│
│  → Output: results/conformer_descriptors_raw.csv    │
└─────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────┐
│  scripts/build_feature_matrix.py                    │
│  → Merge: Tier-1 Δ + DB 3DPSA + RDKit 2D baseline  │
│  → Binary label: permeable (PAMPA ≥ −6.0)          │
│  → Output: results/feature_matrix.csv               │
└─────────────────────────────────────────────────────┘
         │
         ├──────────────────────────┐
         ▼                          ▼
┌──────────────────┐     ┌──────────────────────────┐
│ correlation_     │     │ umap_visualization.py    │
│ analysis.py      │     │ → Panel A: 2D descriptors│
│ → Pearson/Spear  │     │ → Panel B: 3D Δ features │
│ → AUC-ROC        │     │ → Panel C: Combined      │
│ → LR importance  │     │ → Leiden clustering      │
│ → Tables + figs  │     │ → CycloA overlay         │
└──────────────────┘     └──────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────┐
│  scripts/tier2_validation.py                        │
│  → Tier-1 ΔPSA vs DB CHCl3/H2O 3DPSA cross-check   │
│  → Reference set (CycloA + analogs) scatter plots   │
│  → Output: results/tier2_validation_table.csv       │
└─────────────────────────────────────────────────────┘
```

---

## Feature Groups

| Group | Features | Coverage |
|-------|---------|---------|
| **DB 3D** (CycPeptMPDB) | `H2O_3DPSA`, `CHCl3_3DPSA`, `delta_3DPSA_db` | ~88% of PAMPA subset |
| **Tier-1 Δ** (conformer engine) | `delta_psa3d`, `delta_hb`, `delta_Rg`, `delta_NPR1`, `delta_NPR2`, `psa3d_spread`, `psa3d_std` | Computed for all PAMPA entries |
| **2D Baseline** (RDKit) | `MolWt`, `MolLogP`, `TPSA`, `HBA`, `HBD`, `RotBonds`, `CSP3`, `Rings` | 100% |

**Key Δ feature definitions:**
- `delta_psa3d` = PSA(max-PSA conformer) − PSA(min-PSA conformer)  ← chameleonic PSA spread
- `delta_hb` = HB(min-PSA conf) − HB(max-PSA conf)  ← H-bonds formed upon "membrane entry"
- `delta_Rg` = Rg(max-PSA conf) − Rg(min-PSA conf)  ← compaction upon membrane entry
- `psa3d_spread` = max(PSA) − min(PSA) across all conformers  ← conformational flexibility

---

## Outputs

All outputs written to `results/`:

```
results/
├── conformer_descriptors_raw.csv    # Per-molecule Δ descriptors (Tier-1)
├── feature_matrix.csv               # Full merged feature matrix
├── feature_groups.json              # Feature group definitions
├── correlation_table.csv            # Pearson r + Spearman ρ vs. PAMPA
├── auc_roc_table.csv                # AUC-ROC per feature
├── feature_importance.csv           # Logistic regression coefficients
├── umap_panel_summary.csv           # UMAP silhouette + enrichment per panel
├── Panel_A_2D_embedding.csv         # UMAP coordinates (2D panel)
├── Panel_B_3D_delta_embedding.csv   # UMAP coordinates (3D Δ panel)
├── Panel_A_2D_hit_enrichment.csv    # Leiden cluster enrichment (2D panel)
├── Panel_B_3D_delta_hit_enrichment.csv
├── tier2_validation_table.csv       # Reference set cross-check
└── figures/
    ├── fig1_data_overview.png
    ├── correlation_heatmap.png
    ├── auc_roc_bar.png
    ├── scatter_top_features.png
    ├── Panel_A_2D_umap.png
    ├── Panel_B_3D_delta_umap.png
    ├── Panel_C_combined_umap.png
    └── tier2_crosscheck.png
```

---

## Tier-2 High-Rigor Validation (Optional)

For the 5 reference compounds (CycloA + analogs):
```bash
python scripts/tier2_validation.py --matrix results/feature_matrix.csv \
                                    --refset data/reference_set.csv
```

If you have OpenEye OMEGA + academic license, you can replace Tier-1 conformers with OMEGA macrocycle conformers. The `tier2_validation.py` script will compare both approaches.

---

## Thesis Connection

This pipeline is designed to transfer directly to our in-house DEL cyclic peptide library:
- Replace `CycPeptMPDB_Peptide_All (2).csv` with the internal DEL library CSV
- The `conformer_engine.py` and `correlation_analysis.py` scripts work on any CSV with a `SMILES` column

---

## Limitations

1. **Tier-1 dual-dielectric approximation:** MMFF94s does not model explicit dielectric environments. Aqueous/membrane conformer selection by PSA extremes is a heuristic, not physics-based. Tier-2 (OMEGA + GB/SA) validates this.
2. **PAMPA assay heterogeneity:** Multiple labs/protocols contribute to CycPeptMPDB. PAMPA is our primary target; Caco-2/RRCK analyzed separately.
3. **Statistical power:** With ~7,000 PAMPA compounds, even weak correlations (|ρ| > 0.05) will be statistically significant. Report effect sizes, not just p-values.
4. **ETKDGv3 failure rate:** Large macrocycles (>12 residues) may fail embedding. Expected ~5-15% failure rate; reported in pipeline output.

---

## File Structure

```
CHEM_269_Final_Project/
├── CycPeptMPDB_Peptide_All (2).csv       ← place here (not tracked by git)
├── environment.yml                        ← conda environment
├── run_pipeline.py                        ← single-command entry point
├── README.md
├── notebooks/
│   └── 3d_descriptors.ipynb              ← main deliverable notebook
├── scripts/
│   ├── curate_data.py
│   ├── conformer_engine.py
│   ├── build_feature_matrix.py
│   ├── correlation_analysis.py
│   ├── umap_visualization.py
│   └── tier2_validation.py
├── data/                                  ← generated by curate_data.py
├── results/                               ← generated by pipeline
└── figures/                               ← symlink to results/figures
```

---

## References

- Jiang et al. (2023). CycPeptMPDB: A Comprehensive Database of Membrane Permeability of Cyclic Peptides. *J. Chem. Inf. Model.*
- Rezai et al. (2006). Conformational flexibility, internal hydrogen bonding, and passive membrane permeability. *J. Am. Chem. Soc.*
- Witek et al. (2016). Kinetic models of cyclosporin A in polar and apolar environments. *J. Chem. Theory Comput.*
- Bockus et al. (2015). Decoding chameleonic properties of macrocycles. *J. Med. Chem.*
- Riniker et al. (2015). Better informed distance geometry: Using what we know to improve conformation generation. *J. Chem. Inf. Model.* (ETKDGv3 foundation)
