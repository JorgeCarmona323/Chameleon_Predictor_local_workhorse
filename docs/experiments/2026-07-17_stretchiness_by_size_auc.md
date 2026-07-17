# 2026-07-17 — Is vacuum ΔPSA a within-size-class signal, or a size proxy?

**Experiment.** Test whether the vacuum "stretchiness" descriptors (accessible ΔPSA range /
flexibility) predict PAMPA permeability **within a size class**, or whether their apparent
signal is really a between-size-class (size proxy) effect. Directly tests the hypothesis
behind the multi-fidelity plan ([[project_multifidelity_vacuum_implicit]]) and the ΔPSA
caveat in `2026-07-17_umap_report_v2_pi_response.md`.

## Method

Single-feature AUC-ROC vs. the permeability label (`permeable`, PAMPA ≥ −6.0 log cm/s),
computed as the Mann-Whitney statistic (rank-based, average-rank ties). Stratified by
`Monomer_Length`: **all / <9 / ≥9 residues**. Run on the full dataset and on the clean-source
subset (2016_Furukawa + 2013_CHUGAI, the internally consistent PAMPA protocols).

- Data: `results/archive/feature_matrix.csv` (7,298 compounds; all 2D + 3D + DB features).
- Outputs: `results/stretchiness_by_size_full.csv`, `results/stretchiness_by_size_clean.csv`.
- Reproducible script: `scripts/ml/stretchiness_by_size.py` (pandas + sklearn on the HPC/Colab;
  the numbers below were computed locally via the equivalent rank-AUC in PowerShell).

**Method validation:** `delta_psa3d` AUC on the full dataset = **0.51**, reproducing the
known collapsed value (≈0.505) reported for the full 7,297-compound set. The rank-AUC is
correct.

## Results

### Full dataset (n = 7,298; permeable 66%)

| Feature | AUC all | AUC <9 | AUC ≥9 |
|---------|:------:|:-----:|:-----:|
| delta_psa3d | 0.51 | 0.51 | 0.54 |
| psa3d_std | 0.51 | 0.51 | 0.55 |
| psa3d_spread | 0.51 | 0.51 | 0.54 |
| delta_psa3d_per_mw | 0.47 | 0.49 | 0.44 |
| MolLogP | 0.63 | 0.64 | 0.61 |
| MolWt | 0.53 | 0.52 | **0.73** |

*n:* all 7,298 · <9 4,455 · ≥9 2,843.

### Clean subset — Furukawa + Chugai (n = 1,602)

| Feature | AUC all | AUC <9 | AUC ≥9 |
|---------|:------:|:-----:|:-----:|
| delta_psa3d | **0.68** | 0.53 | 0.52 |
| psa3d_std | 0.69 | 0.54 | 0.48 |
| psa3d_spread | 0.68 | 0.53 | 0.52 |
| delta_psa3d_per_mw | 0.59 | 0.55 | **0.66** |
| MolLogP | 0.32 | 0.42 | 0.15 |
| MolWt | 0.62 | 0.43 | 0.25 |

*n:* all 1,602 · <9 737 · ≥9 865. *Permeable rate:* all 75% · <9 60% · **≥9 88%**.

## Findings

1. **Absolute ΔPSA / stretchiness is largely a *size proxy*, not a within-class signal.**
   On the clean subset, `delta_psa3d` scores **0.68 pooled** but collapses to **0.53 (<9)**
   and **0.52 (≥9)** once size is held fixed. The pooled signal is mostly the between-class
   effect that large peptides have both large ΔPSA *and* higher permeability (88% vs 60%).
   **This quantitatively confirms the caveat that absolute ΔPSA scales with size.**

2. **Size-normalization rescues the large-peptide class.** `delta_psa3d_per_mw` is the *only*
   stretchiness feature that **improves** within a size bin — **0.66 in the ≥9 class** vs 0.52
   for absolute ΔPSA. Direct, quantitative support for (a) the Yu et al. 2026 per-size
   normalization and (b) chameleonicity as a **size-gated** phenomenon that emerges in the
   ≥9-residue class once the size confound is removed.

3. **Size-stratification does not rescue the full noisy dataset.** On all 7,298, every
   stretchiness feature stays ≈0.51 in every bin — the pooled-source label noise dominates
   regardless of size. (MolWt ≥9 = 0.73 on full data is the inverse of the clean subset's
   0.25, another fingerprint of source heterogeneity flipping relationships.)

4. **The ≥9 clean bin is 88% permeable** — highly imbalanced, so there is a low ceiling on
   how much *any* single descriptor can add there; most large clean-set peptides are permeable
   regardless.

## Interpretation & decision

- **The "vacuum stretchiness alone" hypothesis is largely negative for the *absolute*
  descriptor** — it carries no within-size-class permeability signal. It works pooled only
  because it encodes size.
- **But the *normalized* capacity feature has real (if modest) within-class signal for the
  large chameleonic peptides (≥9 → 0.66).** So the vacuum layer is not useless; it needs
  (i) size normalization and (ii) the implicit-solvent **propensity** layer to push past ~0.66.
- This *sharpens*, not kills, the multi-fidelity plan: use **normalized** capacity features
  from vacuum as the cheap prior, restricted to the **≥9-residue / clean-source** regime,
  and let the CREST/CPCM-X **propensity** descriptors provide the lift. Vacuum capacity is a
  screen; implicit propensity is the classifier.

## Next

1. Add `delta_psa3d_per_sasa`, `delta_psa3d_per_residue` and re-run this table (per-MW may not
   be the best normalizer; per-SASA is the Yu 2026 form).
2. On the ≥9 clean subset, test the implicit-solvent descriptors (Boltzmann-weighted ΔPSA,
   ΔG_transfer) head-to-head against `delta_psa3d_per_mw`'s 0.66 — the first real
   capacity-vs-propensity comparison.
3. Multivariate model (gradient boosting / TabPFN) on the ≥9 clean subset with the normalized
   feature set, SHAP for which conformational DOF matters — vs. single-feature AUC here.
