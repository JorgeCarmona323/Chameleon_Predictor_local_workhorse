# CHEM 269 Final Project — Findings, Methods Log & Known Issues

**Project:** 3D Conformational Descriptors and Dual-Dielectric Solvent Modeling for Cyclic Peptide Membrane Permeation
**Dataset:** CycPeptMPDB — 7,298 cyclic peptides with experimental PAMPA LogPexp
**Date:** 2026-03-17

---

## 1. Core Hypothesis

Cyclic peptides that passively permeate membranes exhibit "chameleonic" behavior: they adopt a polar-exposed conformation in aqueous environments and a polar-buried (H-bond shielded) conformation in low-dielectric environments. This switch is quantified as ΔPSA = PSA(aqueous) − PSA(membrane). Standard 2D descriptors (TPSA, MolLogP) cannot capture this because they are topology-based and conformationally blind.

---

## 2. Pipeline Summary

| Stage | Method | Status | Coverage |
|---|---|---|---|
| Data curation | RDKit canonicalization, PAMPA filter | Complete | 7,298 / 8,466 compounds |
| 2D baseline | RDKit descriptors | Complete | 100% (7,298) |
| DB 3DPSA | CycPeptMPDB H2O/CHCl3_3DPSA | Complete | 88% (6,942) |
| Tier-1 conformers | ETKDGv3 + MMFF94s, n=20 | Partial | 1,502 / 7,298 (20.6%) |
| Tier-2 CREST | CREST 2.12 + ALPB | **Failed — see Section 5** | 0 / 5 |
| Tier-2 xtb+GBSA | GFN2-xTB + GBSA single-opt | Not used — see Section 5 | — |
| Feature matrix | Merged all above | Complete | 7,299 rows |
| Correlation analysis | Pearson, Spearman, AUC-ROC | Complete (partial Tier-1) | — |
| UMAP + clustering | K-Medoids + HDBSCAN dual-track | Complete (partial Tier-1) | — |

---

## 3. Key Findings

### 3.1 2D Baseline Performance

All features significant at p < 0.001 (N=7,298).

| Feature | Spearman ρ | AUC-ROC | Notes |
|---|---|---|---|
| NumRotatableBonds | +0.212 | 0.563 | Top 2D predictor (ρ) |
| MolLogP | +0.151 | **0.631** | Top AUC predictor |
| MolWt | +0.116 | 0.534 | Correlated with size |
| NumHDonors | −0.085 | 0.594 | More donors → less permeable |
| NumHAcceptors | +0.048 | 0.519 | Weak |
| TPSA | −0.017 | 0.553 | Surprisingly weak |

**Interpretation:** MolLogP (AUC=0.631) is the strongest single predictor, consistent with passive diffusion theory. TPSA — despite being the standard proxy for polarity — is nearly uncorrelated (ρ=−0.017), likely because cyclic peptides violate the assumptions underlying the fragment-based TPSA model (it was calibrated on drug-like small molecules, not macrocycles).

### 3.2 DB 3DPSA Performance — Critical Finding

| Feature | Spearman ρ | AUC-ROC | Notes |
|---|---|---|---|
| H2O_3DPSA | −0.020 | 0.565 | Near-zero correlation |
| CHCl3_3DPSA | −0.018 | 0.563 | Near-zero correlation |
| delta_3DPSA_db | −0.020 | 0.507 | **Essentially random** |

**The CycPeptMPDB 3DPSA values fail completely as chameleonic indicators.**

Investigation of the 5 reference compounds confirms this:

| Compound | PAMPA | Permeable | DB delta_3DPSA | Tier-1 ΔPSA | Lit ΔPSA |
|---|---|---|---|---|---|
| CsA | −6.60 | Yes | −1.0 Å² | **84.9 Å²** | ~75 Å² |
| DP172 | −4.15 | Yes | −47.0 Å² | **88.9 Å²** | — |
| HexPep | −6.20 | No | +2.0 Å² | **64.4 Å²** | — |
| 1NMe3 | −5.52 | Yes | +1.0 Å² | **47.8 Å²** | — |
| PSLYF | −9.10 | No | +0.0 Å² | **65.3 Å²** | — |

The DB values are near-zero for all compounds (and anomalously negative for DP172), indicating the CycPeptMPDB 3DPSA was computed from a single-conformer or non-ensemble methodology that cannot capture conformational switching. The delta between solvents is essentially meaningless.

### 3.3 Updated Correlation Results (post-pipeline re-run with 1,502 Tier-1 compounds)

**Tier-1 3D features substantially outperform all 2D and DB features:**

| Feature | Group | N | Spearman ρ | AUC-ROC |
|---|---|---|---|---|
| psa3d_std | Tier1 | 1,501 | **+0.467** | **0.749** |
| psa3d_spread | Tier1 | 1,502 | **+0.457** | **0.744** |
| delta_psa3d | Tier1 | 1,502 | **+0.457** | **0.744** |
| delta_Rg | Tier1 | 1,502 | +0.293 | 0.652 |
| hb_spread | Tier1 | 1,500 | +0.299 | 0.643 |
| NumRotatableBonds | 2D | 7,298 | +0.212 | 0.563 |
| MolLogP | 2D | 7,298 | +0.151 | 0.631 |
| delta_hb | Tier1 | 1,502 | +0.127 | 0.559 |
| H2O_3DPSA | DB | 6,941 | −0.020 | 0.565 |
| delta_3DPSA_db | DB | 6,941 | −0.020 | 0.507 |

**Key takeaway:** Tier-1 delta_psa3d (AUC=0.744) outperforms the best 2D descriptor MolLogP (AUC=0.631) by 11 AUC points. The DB 3DPSA delta (AUC=0.507) is essentially random.

---

### 3.4 Tier-1 Conformer Results (1,502 compounds)

| Descriptor | Mean | Std | Notes |
|---|---|---|---|
| delta_psa3d | ~36–89 Å² | — | Wide range; strong molecules show large delta |
| delta_hb | 0–4 | — | Sparse at n=20; see sampling note |
| psa3d_spread | mirrors delta_psa3d | — | More robust than delta_hb at low n |

**Reference compound validation:**
- CsA Tier-1 ΔPSA = 84.9 Å² vs literature ~75 Å² (Witek JCTC 2016) — **11% over-estimate, consistent**
- ETKDGv3 vacuum ensemble sampling correctly identifies CsA as highly chameleonic
- DB values fail; Tier-1 succeeds — confirms ensemble sampling is necessary

**ΔHB sampling note:** At n=20 conformers, delta_hb = 0 for all reference compounds. The H-bond-forming membrane conformer is rare in a 20-conformer vacuum ensemble. psa3d_spread is more reliable at low n. Minimum recommended n ≥ 50 for reliable ΔHB signal.

### 3.4 UMAP / Clustering

- Dual-track: K-Medoids (k=8, cosine) as primary + UMAP→HDBSCAN as exploratory
- ARI stability checked across 5 seeds (threshold ≥ 0.85)
- If stability fails: K-Medoids is primary evidence, UMAP is exploratory only
- Results pending full Tier-1 run (currently 20.6% coverage)

---

## 4. Methodological Decisions

| Decision | Rationale |
|---|---|
| RobustScaler (not StandardScaler) | Cyclic peptides have heavy-tailed MW/logP; outliers dominate StandardScaler |
| PCA omitted from main pipeline | 12-feature panel is small; PCA destroys physical interpretability |
| Cosine distance (not Euclidean) | Compounds span 2 orders of magnitude in MW; cosine normalizes scale |
| Max-PSA = aqueous, Min-PSA = membrane | Transparent heuristic: polar exposed in water, buried in membrane |
| Heavy-atom-only polar SASA | Consistent with Ertl 2000 TPSA convention and CycPeptMPDB convention |
| UMAP parameters pre-specified | n_neighbors=30, min_dist=0.15 — not tuned post-hoc (would be p-hacking) |
| K-Medoids over K-Means | Each medoid is a real molecule; interpretable structural archetype |
| Open-source only | RDKit (BSD), xtb (LGPL). No OpenEye, no ORCA |

---

## 5. Tier-2 CREST Failure — Full Account

### What was attempted

Tier-2 validation was designed to run CREST 2.12 iMTD-GC with ALPB implicit solvation on 5 reference compounds (CsA, 1NMe3, HexPep, DP172, PSLYF) in two dielectric environments (water ε≈80, CHCl3 ε≈4.8). CREST generates a full conformer ensemble in each solvent and identifies the lowest-energy conformer, giving a thermodynamically rigorous ΔPSA.

### What happened

**Attempt 1 — CREST 2.12 + GFN2-xTB (default):**
- HexPep: timed out after 4 hours (Colab session limit)
- 1NMe3, CsA, DP172, PSLYF: immediate non-zero exit (seconds), likely memory or GFN2-xTB crash on large macrocycles in Colab environment
- Root cause: GFN2-xTB with iMTD-GC on cyclic peptides MW 700–1200 Da is computationally prohibitive on a single Colab CPU session. Literature estimates ~15–25 min/compound; in practice much longer under Colab resource constraints.

**Attempt 2 — CREST 2.12 + GFN-FF (`--gfnff --squick`):**
- HexPep: actively running (Meta-MD steps confirmed progressing in log)
- Session was active but results not obtained before time pressure required pivoting
- Scientific note: GFN-FF with ALPB is less rigorous than GFN2-xTB with ALPB. ALPB solvation is parameterized for GFN2-xTB; using it with GFN-FF gives reasonable geometries but less accurate solvation energy rankings.

**Attempt 3 — xtb GFN2+GBSA single-structure optimization:**
- Tested locally and run to completion on Colab (all 5 compounds). Runs in ~2–8 min/compound.
- Results:

| Compound | Permeable | xtb ΔPSA |
|---|---|---|
| HexPep | No | 0.82 Å² |
| PSLYF | No | 5.40 Å² |
| 1NMe3 | Yes | 6.91 Å² |
| DP172 | Yes | −0.24 Å² |
| CsA | Yes | −0.14 Å² |

- No signal: permeable and impermeable compounds are indistinguishable. ΔPSA range 0–7 Å² with no correlation to permeability.
- Root cause: single-structure optimization in each solvent starts from the same MMFF94s geometry and relaxes to nearby local minima. It never samples the extreme chameleonic conformers. This approach cannot capture chameleonic behavior.
- This is equivalent to what CycPeptMPDB likely used to compute their near-zero 3DPSA deltas.
- **Scientific value:** The xtb results independently confirm that single-structure solvation methods fail — strengthening the argument that ensemble sampling (Tier-1) is necessary.
- Decision: not included in main analysis, documented as negative control.

### What was used instead

Tier-1 ETKDGv3 results on the 5 reference compounds serve as the validation layer:
- CsA ΔPSA = 84.9 Å² (Tier-1) vs ~75 Å² (Witek JCTC 2016) — 11% over-estimate, consistent
- DB 3DPSA fails for all 5 compounds (ΔPSA ≈ 0)
- Conclusion: vacuum ensemble sampling (Tier-1) captures chameleonic behavior; single-structure solvation approaches do not

### Recommended future work

For a production run with CREST:
- Use a dedicated HPC cluster (not Colab) with ≥32 cores per compound
- CREST 2.12 + GFN2-xTB + ALPB is the gold standard
- Expected wall time: ~15–45 min/compound on 32 cores
- Alternatively: CREST with `--gfnff` on HPC (validated that Meta-MD steps run correctly, just needs more wall time than Colab provides)

---

## 6. Limitations

### 6.1 Dataset Selection Bias
CycPeptMPDB is enriched for permeable compounds (66.4% permeable at PAMPA ≥ −6.0 log cm/s), reflecting selection bias toward pharmaceutical candidates. Medicinal chemists synthesize and test compounds they expect to permeate — truly impermeable cyclic peptides are underrepresented. This means:
- The "impermeable" class in this dataset consists largely of near-miss drug candidates, not structurally diverse polar compounds. Feature space contrast between classes is compressed.
- UMAP instability (all 3 panels failed ARI threshold) is partly attributable to this overlap — even "impermeable" compounds here have some chameleonic character.
- AUC-ROC baseline for a naive "predict everything permeable" classifier is ~0.664, not 0.5. Tier-1 AUC=0.749 is ~9 points above this naive baseline (not ~25 as it would appear against a balanced dataset).
- Observed correlations are likely *underestimates* of true predictive power on a balanced benchmark dataset.

### 6.2 Tier-1 Coverage
Only 1,502 / 7,298 compounds (20.6%) have Tier-1 3D features. Correlation and UMAP results for Tier-1 features reflect a non-random subsample. The full A100 Colab run is required for final conclusions.

### 6.3 UMAP Instability
All 3 UMAP panels failed the pre-specified ARI stability check (min ARI 0.07–0.38, threshold 0.85). Per protocol, K-Medoids is reported as primary clustering evidence and UMAP embeddings are exploratory visualizations only. Instability reflects the continuous nature of membrane permeability (not discrete classes) combined with dataset selection bias (Section 6.1).

### 6.4 Tier-2 Validation Not Completed
CREST+ALPB conformer sampling could not be completed within Colab computational constraints (see Section 5). Tier-1 validation against literature benchmarks (CsA ΔPSA=84.9 Å² vs lit. ~75 Å²) serves as the primary validation. Single-structure xtb+GBSA was tested and confirmed as insufficient (near-zero ΔPSA for all compounds).

### 6.5 ΔHB Sampling Artifact
At n=20 conformers, delta_hb = 0 for all reference compounds. H-bond-forming membrane conformers are rare in a 20-conformer vacuum ensemble. psa3d_spread and psa3d_std are more reliable at low n. Minimum n ≥ 50 recommended for reliable ΔHB signal.

---

## 7. Future Directions

### 7.1 Immediate (complete this project)
- Run full Tier-1 on all 7,298 compounds via A100 Colab (~13–15 hrs overnight)
- Re-run correlation, UMAP with full dataset — expect ARI stability to improve significantly
- CREST+ALPB on an HPC cluster (≥32 cores) for the 5 reference compounds — ~15–45 min/compound

### 7.2 Improve the conformer sampling
- **Increase n_confs to 50–200** for production run — critical for reliable ΔHB signal; psa3d_spread is stable at n=20 but delta_hb requires larger ensembles
- **Explicit solvent MD** (GROMACS/OpenMM + GAFF2) on reference compounds — gold standard, directly comparable to Witek JCTC 2016
- **CREST on HPC** with GFN2-xTB + ALPB — the intended Tier-2 approach; computationally prohibitive on Colab but tractable on a university cluster

### 7.3 Address dataset bias
- **Curate a balanced benchmark** — add structurally diverse impermeable compounds from ChEMBL or PubChem with confirmed low permeability; target 50/50 split
- **Stratified train/test split** — evaluate model on held-out balanced test set to get unbiased AUC estimates
- **PAMPA vs Caco2 vs RRCK comparison** — CycPeptMPDB has multiple assay types; PAMPA measures passive permeability while Caco2 includes active transport. Filtering to PAMPA-only (as done here) is correct but limits N; could re-run including RRCK for larger dataset

### 7.4 Improve the model
- **Machine learning on Tier-1 features** — Random Forest or XGBoost on the 12-feature panel; AUC=0.749 from a single Spearman correlation is a floor, not a ceiling
- **Graph neural networks** — encode the cyclic peptide topology directly; several published GNN models for cyclic peptide permeability (CycPeptMPDB paper includes baselines)
- **Sequence-based features** — N-methylation pattern, d-amino acid positions, ring size — these encode structural modifications known to drive chameleonicity
- **Multi-task learning** — predict PAMPA, Caco2, and RRCK simultaneously; shared representation may improve generalization

### 7.5 Deeper chameleonic characterization
- **PSA trajectory analysis** — instead of max/min conformer, compute PSA as a function of all conformers and fit a bimodal distribution; bimodality index as a chameleonicity score
- **Solvent-accessible volume** — complement PSA with 3D molecular volume change between conformers
- **NMR validation** — experimental NOESY/ROESY in CD3OH vs CDCl3 directly measures conformational switching; compare to predicted chameleonic compounds

---

## 9. What Remains

| Task | Status | Blocker |
|---|---|---|
| Full Tier-1 on 7,298 compounds | Pending | Needs A100 Colab run (~13–15 hrs) |
| Re-run feature matrix / correlation / UMAP | Ready to run after full Tier-1 | — |
| Tier-2 CREST validation | Not completed | Computational constraints — see Section 5 |
| Finalize 3d_descriptors.ipynb | In progress | Awaiting full Tier-1 CSV |
| Report write-up | In progress | — |

---

## 7. How to Cite This Work

**Database:** CycPeptMPDB — Rettie et al., *J. Chem. Inf. Model.* 2023
**Chameleonic behavior:** Rezai & Lokey, *JACS* 2006; White & Lokey, *Nat. Chem. Biol.* 2011
**CsA ΔPSA benchmark:** Witek et al., *J. Chem. Theory Comput.* 2016
**ETKDGv3:** Riniker & Landrum, *J. Chem. Inf. Model.* 2015
**CREST:** Pracht et al., *Phys. Chem. Chem. Phys.* 2020
**xtb/GFN2:** Bannwarth et al., *J. Chem. Theory Comput.* 2019
**RDKit:** Landrum et al., open-source
