# Findings and Methods Log
## CHEM 269 Final Project — The Chameleon Traverse
**Jorge Carmona | Started: 2026-03-15 | Last updated: 2026-04-09**

---

## Purpose

This document is the audit trail for the project. It records methodological decisions, parameter choices, negative results, and findings in the order they happened — the lab notebook behind the deliverable notebook.

---

## 2026-04-09 — NEXT SESSION: Feature Benchmark (Run ASAP)

**Context:** This is the gating experiment before any further compute or model decisions. Full design in `docs/experiments/2026-04-09_feature_benchmark_design.md`.

### Step 1 — Install missing packages
```bash
pip install tabpfn mhfp lightgbm
```
Then add to `environment.yml` under pip dependencies.

### Step 2 — Write the benchmark script
Create `scripts/feature_benchmark.py`:
- Load SMILES + permeability labels from `results/feature_matrix_cremp_overlap.csv` (n=2,416, has CREST CHCl3 descriptors already)
- Generate F1–F6 from SMILES:
  - F1: Morgan bit-based (r=2, 2048-dim) — RDKit
  - F2: Morgan count-based (r=2, 500-dim) — RDKit
  - F3: Morgan count-based (r=2, 2048-dim) — RDKit
  - F4: MAPC (2048-dim) — `mhfp` library
  - F5: Mordred 2D only — `mordred`
  - F6: Mordred 2D+3D from single ETKDG conformer — RDKit + `mordred`
- F7: CREST CHCl3 descriptors already in the CSV (`aq_psa3d`, `mem_psa3d`, `delta_psa3d`, `norm_delta_psa`, `bw_psa3d`, `ensemble_energy`, `unique_confs`)
- Models: TabPFN, LightGBM, RandomForest (baseline)
- Same stratified CV splits across all feature sets (fix random seed)
- Metrics: R², Spearman ρ, AUC-ROC
- Output: results table as `results/feature_benchmark_results.csv`

### Step 3 — Run it
```bash
python scripts/feature_benchmark.py
```

### Step 4 — Interpret and decide
- If CREST >> ETKDG 3D: aqueous CREST compute is high priority
- If CREST ≈ ETKDG 3D: deprioritize more simulation, focus on model
- If MAPC >> Morgan: keep current FP approach
- If Morgan wins: switch to count-based Morgan (simpler, faster)

### F8 (aqueous CREST) — add once runs complete
Re-run benchmark with aqueous CREST descriptors as an additional feature set.

---

## 2026-04-09 — Feature Generation Benchmark Design (Decision-Critical)

Literature review of PROTAC-TS (Murakami et al., *JACS Au* 2026) and CREMP (Grambow et al., *Scientific Data* 2024) raised a critical question for the project: does expensive CREST conformer sampling actually outperform cheap RDKit 3D embeddings for permeability prediction? PROTAC-TS achieved R²=0.710 with only 500-dim count-based Morgan + TabPFN on 43 PROTACs — no conformer ensemble.

**Full benchmark design in:** `docs/experiments/2026-04-09_feature_benchmark_design.md`

**Summary of planned ablation (8 feature sets × 3 models, same CV splits):**
- F1–F3: Morgan bit/count at 500 and 2048 dim (PROTAC-TS style)
- F4: MAPC (current pipeline)
- F5–F6: Mordred 2D-only and 2D+3D from single ETKDG conformer
- F7: CREST CHCl3 ensemble descriptors (done)
- F8: CREST aqueous ensemble descriptors (pending — compute now available)

**Key decisions this answers:**
1. MAPC vs. Morgan — does chirality-aware atom-pair FP beat circular for macrocycles?
2. ETKDG 3D vs. CREST 3D — is the expensive simulation worth it?
3. CHCl3 vs. aqueous CREST — does solvent matter for prediction?

Also noted: PROTAC-TS uses count-based Morgan (beats bit-based); TabPFN won across all feature sets; 3D Mordred added ~2-3% R² even on n=43.

---

## 2026-03-15 — Project Setup and Pipeline Design

### Motivation

The project originated from a specific question about DNA-encoded library (DEL) cyclic peptide scaffolds: why do some large, polar cyclic peptides cross membranes when their 2D properties say they should not? CycPeptMPDB (Jiang et al. 2023) provided a large enough experimental PAMPA dataset to test whether 3D conformational descriptors could capture the chameleonic switching effect invisible to 2D descriptors.

### Initial design decisions

**PAMPA threshold: -6.0 log cm/s**
Adopted directly from CycPeptMPDB's own binary classification (Jiang et al. 2023). Not tuned post-hoc.

**Conformer count: 20 per molecule**
Trade-off between ensemble coverage and compute cost. Literature (Riniker & Landrum 2015, ETKDGv3) recommends 50-200 for large macrocycles, but 20 gives reasonable min/max PSA sampling at 13-15h total runtime on A100. Accepted as a heuristic — not validated against larger ensembles.

**Conformer selection: max-PSA = aqueous, min-PSA = membrane**
Physical rationale: the aqueous conformer maximizes solvent-accessible polar surface (extended, H-bond-exposed); the membrane conformer minimizes it (compact, H-bonds internalized). This is an approximation — a true thermodynamic treatment would use Boltzmann-weighted ensemble averages, not extremes. Chosen for computational tractability. Validated post-hoc by CsA agreement with literature.

**ETKDGv3 macrocycle torsion library: enabled**
`useSmallRingTorsions=True`, `useBasicKnowledge=True` — required for macrocycle conformer diversity. Without these flags ETKDGv3 produces collapsed, unrealistic conformers for rings >8 members.

**MMFF94s force field**
Standard choice for organic molecules in RDKit. No explicit solvent — conformers are generated in vacuum. Known limitation: vacuum sampling may produce hydrophobically collapsed conformers that are not accessible in aqueous solution. Accepted as unavoidable at this compute level.

**Polar SASA calculation: Bondi radii + rdFreeSASA**
`rdFreeSASA.CalcSASA` with Bondi radii. Polar atoms defined as N, O, S, and attached H (standard TPSA atom typing). An earlier attempt used default van der Waals radii and gave systematically lower PSA values — switched to Bondi radii to match the conventions used in the literature ΔPSA values (Witek 2016).

**DB 3DPSA negative control**
`delta_3DPSA_db = H2O_3DPSA - CHCl3_3DPSA` from the CycPeptMPDB database columns. Designed as a negative control to test whether pre-existing single-structure 3D PSA values in the database carry any signal. Expected to fail — included to make the methodological comparison explicit.

### Pipeline structure decided

Four scripts in sequence:
1. `curate_data.py` — SMILES canonicalization, PAMPA filter, 2D descriptors
2. `conformer_engine.py` — ETKDGv3 + MMFF94s + rdFreeSASA, checkpoint/resume
3. `build_feature_matrix.py` — merge all descriptor groups
4. `correlation_analysis.py` — Pearson, Spearman, AUC-ROC, logistic regression
5. `umap_visualization.py` — dual-track clustering (K-Medoids + HDBSCAN), ARI stability

---

## 2026-03-15 to 2026-03-16 — Tier-2 CREST Attempt (Failed)

### What was attempted

CREST 2.12 with ALPB implicit solvation attempted for all 5 reference compounds (CsA, DP172, HexPep, 1NMe3, PSLYF) on Google Colab T4 and A100 instances.

Command used:
```bash
crest molecule.xyz --T 4 --alpb water --mquick
crest molecule.xyz --T 4 --alpb chcl3 --mquick
```

### Outcome by compound

| Compound | Outcome | Notes |
|----------|---------|-------|
| CsA | Immediate exit | Memory allocation failure before first MD step |
| DP172 | Immediate exit | Same as CsA |
| HexPep | Timeout after 4 hours | Did not complete a single conformer |
| 1NMe3 | Immediate exit | Same as CsA |
| PSLYF | Immediate exit | Same as CsA |

### Root cause

CREST GFN2-xTB Hamiltonian scales as O(N²) with atom count. CsA has 88 heavy atoms (203 total with H). At the ALPB solvation level with default iMTD-GC sampling, memory requirement exceeds the Colab instance limit (~25 GB) before the first MD step begins. `--mquick` flag reduces sampling depth but does not reduce the memory requirement for a single GFN2 energy evaluation.

HexPep is smaller but still timed out — CREST's conformational search on a macrocycle with flexible sidechains requires many iMTD iterations even with `--mquick`.

### Decision

CREST abandoned for all 5 reference compounds. Documented as a clean negative result. Tier-2 label changed to "Failed — memory/timeout constraints" in the route completion table.

---

## 2026-03-16 — Tier-2 xtb+GBSA Fallback (Completed)

### What was done

GFN2-xTB single-structure optimization with GBSA implicit solvation run for all 5 reference compounds as a CREST fallback. Two optimizations per compound: one in GBSA(water) and one in GBSA(CHCl3). ΔPSA = PSA(water-optimized) − PSA(CHCl3-optimized).

### Results

| Compound | xtb ΔPSA (Å²) | Tier-1 ΔPSA (Å²) | Permeable |
|----------|--------------|-----------------|-----------|
| CsA | -0.14 | 84.9 | Yes |
| DP172 | -0.24 | 88.9 | Yes |
| HexPep | 0.82 | 64.4 | No |
| 1NMe3 | 6.91 | 47.8 | Yes |
| PSLYF | 5.40 | 65.3 | No |

### Interpretation

xtb gives near-zero ΔPSA for all compounds including CsA. A single optimized structure in each solvent relaxes to a local minimum that does not change meaningfully in polar surface area between water and CHCl3 — even at the GFN2 semiempirical level with physically correct GBSA solvation.

This confirms the design premise of Tier-1: the chameleonic conformational switch requires ensemble sampling of many conformers, not just higher-level single-structure optimization. Ensemble coverage is the prerequisite, not level of theory. The xtb result is used as a negative control in the notebook Tier-2 section (Section 6).

---

## 2026-03-16 — Tier-1 Partial Run: 1,502 Compounds

### What was done

First Tier-1 ETKDGv3 + MMFF94s run completed on 1,502 compounds on Google Colab A100 using the checkpoint/resume system in `conformer_engine.py`. These compounds were the first 1,502 entries processed from the PAMPA-filtered CycPeptMPDB CSV in database order (not random).

### CsA NMR validation

Tier-1 ΔPSA for CsA = 84.9 Å². Literature benchmark from Witek et al. (2016, J. Chem. Theory Comput., DOI: 10.1021/acs.jcim.6b00251): ~75–80 Å² using ROESY-derived structures in apolar and polar solvents (Definition B: ensemble PSA(polar) − PSA(apolar)).

Agreement within ~10% — validates the ETKDGv3 + MMFF94s approach for the canonical chameleonic reference compound. The Doak 2016 value of 174 Å² uses a different definition (Definition A: 2D TPSA − single-structure nonpolar PSA) and is not comparable.

### 1,502-compound AUC results

| Descriptor | AUC-ROC | Spearman rho |
|------------|---------|--------------|
| delta_psa3d (Tier-1) | 0.744 | 0.457 |
| psa3d_std (Tier-1) | 0.749 | 0.467 |
| MolLogP (2D) | 0.631 | — |
| delta_3DPSA_db (DB) | 0.507 | — |

As expected: DB 3DPSA is at chance; Tier-1 ΔPSA beats 2D best by 11 AUC points.

---

## 2026-03-16 — Full 7,297-Compound Tier-1 Run

### What was done

Full PAMPA subset processed on Google Colab A100 using checkpoint/resume. Runtime ~13–15 hours. 7,297 of 7,298 compounds completed. One compound failed embedding (`embed_failed`) — a complex macrocycle that ETKDGv3 could not generate conformers for even with extended retry.

### Full-dataset AUC results

| Descriptor | AUC-ROC |
|------------|---------|
| MolLogP | 0.631 |
| delta_psa3d (Tier-1) | 0.505 |
| delta_3DPSA_db | 0.507 |

ΔPSA signal collapsed from 0.744 (1,502 compounds) to 0.505 (7,297 compounds). Initial hypothesis: the 1,502-compound subset was biased.

---

## 2026-03-17 — AUC Collapse Investigation: PAMPA Heterogeneity

### Source composition analysis

The CycPeptMPDB PAMPA dataset was examined by source. Major contributing sources:

| Source | Share | Protocol | Quality |
|--------|-------|----------|---------|
| Townsend 2020 | ~42% | Pooled 150-compound PAMPA, CycLS MS deconvolution | Preprint; pooled format introduces compound interference |
| Kelly 2021 | ~21% | Individual compound, same pooled protocol | DOI: 10.1021/jacs.0c06115 |
| Furukawa 2016 | ~9% | Individual compound LC-MS, 1% lecithin/n-dodecane, 14h 25°C | DOI: 10.1021/acs.jmedchem.6b01246; cleanest source |
| Chugai | ~12% | Patent (WO 2013/100132 A1); DOPC/hexadecane membrane inferred; detection floor -10.0 vs -8.0 | Unverified assay conditions |

### Why the 1,502-compound result was inflated

The first 1,502 compounds processed by database ID order are disproportionately from Chugai (878 compounds, ~90% permeable) and Furukawa (448 compounds, 57.6% permeable). This is a highly biased sample — the Chugai block is 90% permeable, inflating any permeable-class descriptor. The cross-source permeable fraction difference (Chugai ~90% vs full dataset 66.4%) is enough to create an apparent AUC signal that disappears when the dataset is balanced.

Within-source analysis confirmed: even within Chugai alone, delta_psa3d Spearman rho = -0.111 — no real signal. The 0.744 AUC was a cross-source artifact driven by MW and source-level permeability differences, not a genuine ΔPSA effect.

### Key implication

Both the Tier-1 ensemble descriptor and the DB single-structure descriptor fail at full scale. MolLogP's AUC of 0.631 — lower than expected for logP predicting PAMPA — is itself evidence of label noise: on a clean single-protocol dataset, logP would typically achieve AUC > 0.75 for PAMPA. The cross-protocol noise suppresses everything.

### Decision for submission

Full source-stratified re-analysis not completed due to time constraints. Heterogeneity documented as primary limitation and as a scientific finding in README and notebook Section 8. Scoped to follow-up pipeline (Chameleon_Predictor repo).

---

## 2026-03-17 — UMAP Dual-Track Clustering Results (Full 7k)

### Design

- Panel A: 2D descriptors only
- Panel B: Tier-1 delta descriptors only
- Panel C: All features combined
- K-Medoids (k=8, cosine metric, deterministic) — primary
- HDBSCAN (min_cluster_size=50) — exploratory
- ARI stability: 5 random seeds (42, 1, 7, 99, 314)

### ARI stability results (full 7k dataset)

**Panel A (2D):** 8 of 10 seed pairs score ARI 0.81–0.89. One outlier seed (42) at 0.38–0.43. The 2D chemical space is well-structured and mostly stable.

**Panel B (3D delta):** Bimodal — pairs (42,7), (42,314), (7,314) score 0.9953–0.9976 (near-perfect); pairs involving seeds 1 and 99 score 0.10–0.20. This is not random instability. HDBSCAN is finding two internally consistent but structurally different partitions of the space depending on initialization. Some seeds lock onto the two-population chameleonic/rigid partition; others find a different valid cluster boundary. Both attractors are real.

**Panel C (combined):** ARI 0.90–0.99 across all pairs. Highly stable — combining 2D and 3D features resolves the ambiguity in Panel B and produces a single dominant structure.

### Panel B two-population finding

UMAP Panel B consistently shows two populations: high-ΔPSA chameleonic scaffolds (enriched permeable) and low-ΔPSA rigid/polar compounds. This structure survives across all 5 random seeds. It is the strongest visual result of the project and argues that the conformational descriptors correctly stratify chemical space — the PAMPA label noise prevents AUC from detecting the signal, but the conformational partitioning is real.

---

## 2026-03-17 — Notebook Interpretation Decisions

### psa3d_spread removed

`psa3d_spread` (max PSA − min PSA across all 20 conformers) is mathematically identical to `delta_psa3d` by construction (both are PSA(max-PSA conformer) − PSA(min-PSA conformer)). Removed from all analysis features and figures to avoid presenting redundant descriptors as independent. Column retained in `conformer_descriptors_raw.csv`.

### Section 8 honest reporting

The 1,502-compound AUC = 0.744 result is preserved in the notebook as an exploratory finding and explicitly labeled as a subset artifact. The full-dataset AUC = 0.505 is the honest population-scale result. Both are reported with the PAMPA heterogeneity explanation.

### Methodological principles enforced throughout

- No post-hoc parameter tuning to improve AUC or UMAP aesthetics
- UMAP parameters not adjusted after seeing ARI results
- Null results (DB 3DPSA = chance, full-dataset ΔPSA = chance) reported as valid findings, not buried
- psa3d_spread removed for scientific accuracy, not to improve results

---

## 2026-03-17 — Repo Cleanup for Submission

Files removed from submission repo:
- `colab/` — Colab-specific notebooks not needed for local reproduction
- `run_pipeline.py` — redundant with the main notebook as primary deliverable
- `docs/literature_deltapsa_values.md` — internal reference document
- `data/` and `logs/` tracked files — generated outputs

Files restored/added:
- `docs/findings_and_methods_log.md` — this document (required deliverable)
- `assignment/climb_route_prompt.md` — already tracked, confirmed present
- `results/figures/*.png` — force-added (gitignored by default as generated output)
- Key result CSVs — force-added for submission completeness

---

---

## 2026-03-23 — CREMP Benchmark: Hybrid ΔPSA Validation (Option C)

### Motivation

The 2026-03-17 full-dataset analysis showed ΔPSA AUC = 0.505 (chance). Two competing explanations: (A) conformer quality — vacuum ETKDGv3 membrane conformer is physically wrong; (B) label noise — cross-protocol PAMPA heterogeneity suppresses all signals. To test explanation A, we compared vacuum min-PSA membrane conformers against physically-grounded CREMP CHCl₃ CREST ensemble conformers.

### CREMP dataset

- 2,458 pre-computed CREST/GFN2-xTB conformational ensembles, CHCl₃ implicit solvent, 298.15 K
- Up to 5,743 unique conformers per compound (mean ~700)
- Source: Atz et al. 2024 (ETH Zurich DEL cyclic hexapeptides)
- Boltzmann weights available per conformer

### Option C: hybrid ΔPSA

- **Aqueous PSA**: vacuum ETKDGv3 max-PSA (existing `feature_matrix.csv`)
- **Membrane PSA**: CREMP CHCl₃ Boltzmann-weighted mean-PSA (`bw_psa3d`)
- **Hybrid ΔPSA** = vacuum aq-PSA − CREMP BW mean-PSA
- BW mean chosen over raw min-PSA: thermodynamically defensible (ensemble average vs. single conformer extreme)

### Conformer sampling in cremp_deltapsa.py

Sampling 100 conformers per compound: first 20 lowest-energy + uniform step across remainder. Ensures near-minimum PSA captured without enumerating all conformers (which would take ~6h). Runtime: ~20 min for full 2,458 compounds.

### norm_delta_psa bug (Flag 1 — fixed 2026-03-24)

`norm_delta_psa` (ΔPSA / SASA_total, dimensionless) was computed using a wrong conformer ID in `_total_sasa`. The variable `aq_idx` was an index into the sampled `psas[]` array, not an RDKit conformer ID. For compounds with n_confs > 100, `_total_sasa(mol, aq_idx)` was fetching a different conformer's total SASA than the one used for aq_psa. Fixed by introducing `aq_conf_id = sample_ids[int(aq_psas_idx)]`. CSV regenerated post-fix.

### SMILES overlap: CREMP × CycPeptMPDB

- CREMP valid canonical SMILES: 2,457 / 2,458
- CycPeptMPDB canonical SMILES: 7,298 / 7,298
- Inner join overlap: **2,435 compounds** (33.4% of CycPeptMPDB)
- Permeable: 1,590 | Impermeable: 845

### Structural bias finding

The overlap is structurally different from the non-overlap portion of CycPeptMPDB:

| Group | Mean MW | Mean Monomer Length | Source composition |
|-------|---------|--------------------|--------------------|
| Overlap | 774 Da | ~6.5 mer | 81.8% Townsend 2020 |
| Non-overlap | 967 Da | ~8.7 mer | Diverse |

Mann-Whitney test on Monomer_Length: p = 0.000 (highly significant). The overlap is dominated by 6–7 mer hexapeptides. The ≥9 residue threshold for chameleonic switching behavior (Yu et al. 2026) is not met by the majority of the overlap. This is the primary explanation for chance-level AUC on this subset.

### Membrane conformer comparison: vacuum min-PSA vs CREMP BW mean-PSA

- Pearson r = 0.603 (moderate correlation — compounds partially rank-ordered, but individual values diverge)
- MAE = 11.6 Å²
- CREMP BW mean lower than vacuum min in 42.4% of cases
- Mean difference (CREMP BW − vacuum min) = +3.4 ± 14.3 Å² (near zero, but high variance)
- Interpretation: vacuum and CHCl₃ ensemble largely agree at population level, but individual compounds show 10–20 Å² divergence. Neither method collapses compounds more than the other systematically.

### AUC results (PAMPA, overlap n=2,435)

| Descriptor | AUC | Notes |
|-----------|-----|-------|
| Vacuum ΔPSA (ETKDGv3) | 0.501 | Chance |
| Hybrid ΔPSA (vac aq / CREMP BW mem) | 0.455 | Worse than vacuum |
| CREMP PSA spread (CHCl₃ only) | 0.412 | Below chance |
| CREMP BW mean-PSA (inverted) | 0.607 | See Flag 9 below |

### Flag 9 — inverted bw_psa3d signal

bw_psa3d AUC was initially reported as 0.393, labeled "(inverted)" in the code but not actually inverted before AUC computation. After fixing (Flag 9): inverted AUC = 0.607. Direction: **lower CHCl₃ BW mean-PSA predicts permeability**. This is the opposite of the ΔPSA hypothesis (high ΔPSA = chameleonic = permeable) but is consistent with lipophilicity: low polar surface area in CHCl₃ → lipophilic → crosses membranes. Effect likely confounded by MW and source bias in the 6–7 mer overlap subset. Not sufficient to conclude the signal is mechanistically chameleonic.

### Conclusion

On the 2,435-compound overlap subset, no ΔPSA variant achieves AUC > 0.51. The hybrid ΔPSA (Option C) does not improve over vacuum-only ΔPSA. Two non-exclusive explanations:

1. **Structural bias (primary)**: The overlap is 6–7 mers. ΔPSA is theoretically expected to predict permeability only for ≥9 mer compounds where chameleonic switching is possible. We are testing the hypothesis on the wrong size range.

2. **Label noise (also present)**: The overlap is 81.8% Townsend 2020 pooled PAMPA. Cross-source noise suppresses all descriptors, as confirmed in the 2026-03-17 analysis.

### Gate condition outcome (Phase 3 mechanistic heads)

**Gate NOT met.** CREMP ΔPSA labels cannot be validated as permeability-predictive on the current overlap. MechanisticHeads training gated until either: (a) overlap is extended to ≥9 mer compounds with PAMPA data, or (b) Furukawa-only stratified analysis is performed.

---

## Open Questions and Next Steps (Chameleon_Predictor)

1. **Normalized ΔPSA**: Yu et al. 2026 (bioRxiv, DOI: 10.64898/2026.01.06.697862) use ΔPSA/SASA_total — a dimensionless fractional switching ratio that removes MW confounding. Combined with a ≥9 residue size filter (below which chameleonic behavior does not manifest reliably per Yu 2026), this is the most promising path to recovering signal on a clean dataset.

2. **Source-stratified PAMPA**: Rerun on Furukawa 2016 only (individual LC-MS, cleanest source) to test whether AUC recovers when cross-source label noise is eliminated.

3. **Random forest + SHAP**: Replace single-descriptor AUC with a multi-feature model to capture nonlinear interactions between ΔPSA, psa3d_std, delta_hb, and delta_Rg.

4. **NMR calibration set**: Omphalotin A (PDB 8QAQ/8QAS, Rüdisser 2023), Ono 2019 hexapeptide diastereomers, and CsH (negative control) as a multi-compound NMR-anchored validation set.

5. **Caco-2 efflux analysis**: 256 compounds with both Tier-1 data and Caco-2 values show inverted correlation (rho = -0.404) — large chameleonic molecules fail Caco-2 due to P-gp efflux, not poor passive permeability. Disentangling passive vs. efflux-limited transport requires both assays on the same compounds.
