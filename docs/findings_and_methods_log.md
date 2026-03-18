# Findings and Methods Log
## CHEM 269 Final Project — The Chameleon Traverse
**Jorge Carmona | Updated: 2026-03-17**

---

## Overview

This log tracks key methodological decisions, findings, and negative results in chronological order. It is intended as a complete audit trail for reproducibility and honest reporting.

---

## 2026-03-16 — Full 7,297-Compound Tier-1 Run Completed

### What was done

The Tier-1 ETKDGv3 + MMFF94s conformer pipeline was run to completion on the full PAMPA subset of CycPeptMPDB using Google Colab A100 GPU (CPU-parallel via joblib). The run used 20 conformers per molecule with the checkpoint/resume system in `conformer_engine.py`.

- Total PAMPA subset: 7,298 compounds
- Successfully completed: 7,297 (one compound failed embedding — `embed_failed`, likely a complex macrocycle with RDKit ETKDGv3 limitations)
- Coverage: 99.99%
- Runtime: approximately 13–15 hours on A100 with `n_confs=20`

### Full-dataset AUC results

After merging conformer descriptors and running `correlation_analysis.py` on all 7,297 compounds:

| Descriptor | AUC-ROC | Notes |
|------------|---------|-------|
| MolLogP | 0.631 | Best single descriptor (2D baseline) |
| delta_psa3d (Tier-1) | 0.505 | Ensemble ΔPSA — at chance on full dataset |
| delta_3DPSA_db | 0.507 | DB single-structure PSA — also at chance |

These results supersede the earlier 1,502-compound exploratory run (AUC = 0.744 for delta_psa3d), which was a non-random biased subset. The earlier result is acknowledged in notebook Section 8 as exploratory.

### psa3d_spread removed from analysis features

`psa3d_spread` (max PSA − min PSA across all 20 conformers) is mathematically equivalent to `delta_psa3d` (which is defined as PSA(max-PSA conformer) − PSA(min-PSA conformer)). Both quantities are identical by construction. `psa3d_spread` was removed from the correlation analysis feature set in `correlation_analysis.py` (added to the `DROP` set) to avoid presenting redundant descriptors as independent features. The column remains in `conformer_descriptors_raw.csv` for reproducibility.

---

## 2026-03-17 — PAMPA Assay Heterogeneity Finding

### Discovery

During interpretation of the full-dataset AUC collapse (from 0.744 on 1,502 compounds to 0.505 on 7,297 compounds), the source composition of the CycPeptMPDB PAMPA dataset was examined.

### Key finding

The PAMPA measurements in CycPeptMPDB are aggregated from at least three incompatible experimental protocols:

1. **Townsend 2020 (bioRxiv preprint)** — approximately 42% of the PAMPA data. Uses a pooled compound format (multiple peptides tested in the same well) rather than individual compound testing. This is non-standard for PAMPA and introduces compound-compound interference effects that do not reflect true passive permeability.

2. **Kelly 2021** — a large fraction of the remaining data. Uses individual compound testing with a standard lipid bilayer membrane formulation (DOPC in hexadecane). More comparable to pharmaceutical-standard PAMPA.

3. **Chugai dataset** — uses a proprietary membrane formulation that differs from both of the above. The LogPexp values from Chugai are not directly comparable to those from Kelly 2021 or Townsend 2020.

### Implications

Aggregating PAMPA measurements from incompatible protocols introduces substantial label noise. A compound measured at −5.5 log cm/s in one protocol might score −7.2 in another, flipping its binary permeable/impermeable label. This cross-protocol label noise is the most likely explanation for:

- The AUC collapse from the 1,502-compound subset (which may have been dominated by one source) to the full dataset
- The failure of both Tier-1 ΔPSA and DB 3DPSA to exceed chance at full scale
- The strong residual AUC for MolLogP (0.631): logP is the most physically stable predictor across protocols because it captures hydrophobicity, which is correlated with permeability regardless of the specific PAMPA membrane used

### Decision for course submission

A full source-stratified re-analysis (e.g., restricting to Kelly 2021 data only) was not completed for the course submission due to time constraints. The heterogeneity is acknowledged as the primary limitation in README.md and in the notebook Limitations section. This is the highest-priority follow-up analysis.

---

## 2026-02-XX — 1,502-Compound Exploratory Run (Superseded)

### Context

The first Tier-1 run was completed on 1,502 compounds — approximately the first 20% of the PAMPA subset, processed in alphabetical order by CycPeptMPDB compound ID. This was not a random or stratified sample.

### Results (exploratory subset only)

- delta_psa3d AUC = 0.744
- Spearman ρ = 0.457

### Why this result does not generalize

The 1,502 compounds processed first by ID ordering were not representative of the full chemical space. Compounds with lower IDs (earlier database entries) tend to be earlier literature compounds that may be systematically different in scaffold type, publication source, or assay protocol compared to later-added compounds. The biased sampling artificially inflated the apparent predictive power of delta_psa3d.

This result is preserved in notebook Section 8 with the label "exploratory subset — do not generalize."

---

## 2026-02-XX — CREST Tier-2 Attempt (Failed)

### What was attempted

CREST 2.12 with ALPB implicit solvation was attempted for all 5 reference compounds (CsA, DP172, HexPep, 1NMe3, PSLYF) on Google Colab (T4 and A100 instances).

### Outcome

All 5 attempts failed:

- **HexPep**: timed out after 4 hours without completing a single conformer
- **CsA, DP172, 1NMe3, PSLYF**: process exited immediately with memory allocation errors

### Root cause

CREST's GFN2-xTB Hamiltonian scales as O(N²) with atom count. CsA has 88 heavy atoms (203 with H). At the tight ALPB solvation level with default CREST iMTD-GC sampling, the memory requirement exceeds the Colab instance limit (~25 GB) before the first molecular dynamics step.

### Fallback: xtb + GBSA (completed)

GFN2-xTB single-structure optimization with GBSA solvation (water and CHCl3) was completed for all 5 reference compounds. Results:

| Compound | xtb ΔPSA (Å²) | Tier-1 ΔPSA (Å²) |
|----------|--------------|-----------------|
| CsA | −0.14 | 84.9 |
| DP172 | −0.24 | 88.9 |
| HexPep | 0.82 | 64.4 |
| 1NMe3 | 6.91 | 47.8 |
| PSLYF | 5.40 | 65.3 |

The xtb result is itself informative: single-structure semiempirical optimization with implicit solvation gives near-zero ΔPSA for all compounds, including CsA. This confirms that the chameleonic effect requires conformational ensemble sampling — it cannot be captured by optimizing a single structure in each solvent, regardless of the level of theory. This is the negative control that validates the design premise of Tier-1.

---

## 2026-02-XX — CsA NMR Validation

Tier-1 ΔPSA for CsA = 84.9 Å². Literature value from Witek et al. (2016, J. Chem. Theory Comput.) using ROESY in apolar/polar solvents: ~75–80 Å² (Definition B: ensemble PSA(polar) − PSA(apolar)).

The Tier-1 value is within ~10% of the literature benchmark, validating the ETKDGv3 + MMFF94s approach for the canonical chameleonic reference compound.

**Important note on ΔPSA definitions:** Three incompatible ΔPSA definitions exist in the literature (see `docs/literature_deltapsa_values.md`). The Tier-1 calculation uses Definition B (ensemble max-PSA minus min-PSA conformer). The Doak 2016 value of 174 Å² uses Definition A (2D TPSA minus single-structure nonpolar PSA) and is not directly comparable.

---

## 2026-02-XX — DB 3DPSA Negative Control Confirmed

AUC-ROC for `delta_3DPSA_db` (CycPeptMPDB database H2O\_3DPSA − CHCl3\_3DPSA) = 0.507 on the full dataset. This is the expected result: the database values come from single optimized structures and cannot capture the chameleonic conformational switch. They represent the aqueous and CHCl3 single-energy-minimized structures, which differ very little in PSA because the MMFF force field without explicit solvation does not drive the conformational collapse that buries polar groups in nonpolar media.

This result is robust to dataset size (consistent across 1,502-compound and 7,297-compound runs) and serves as the primary methodological comparison: conformer ensemble sampling is necessary but not sufficient for population-scale PAMPA prediction given the current dataset's assay heterogeneity.

---

## 2026-02-XX — UMAP Panel B Two-Population Finding

UMAP Panel B (3D Δ feature space) consistently shows two populations across the full 7,297-compound dataset: a compact high-ΔPSA cluster (enriched in permeable compounds) and a diffuse low-ΔPSA region. This structure survives across all 5 random seeds tested in the stability analysis, though the HDBSCAN cluster boundary positions are not stable (ARI 0.07–0.38 pairwise).

**Interpretation:** The two-population structure is a real feature of the 3D Δ descriptor space — conformationally flexible compounds (high ΔPSA, high delta_Rg) are systematically separated from conformationally rigid ones. That this does not translate to a high AUC-ROC is consistent with the PAMPA label noise hypothesis: the chameleonic compounds are correctly identified in conformational space, but their permeability labels are noisy across sources.

The poor ARI stability indicates that specific cluster assignments should not be used as scientific evidence — but the visual separation in Panel B (PAMPA LogPexp subplot) is robust and the strongest claim supported by the full-dataset analysis.
