# Experiment — Dataset Coverage: Labels vs Ensembles by Residue Count

**Date:** 2026-06-09
**Script:** `scripts/analyze_dataset_coverage.py`
**Output:** `results/dataset_coverage.csv`
**Data:** `results/archive/feature_matrix.csv` (CycPeptMPDB) · `results/archive/cremp_deltapsa.csv` (CREMP)

---

## Question

The hypothesis says chameleonicity matters most above ~9 residues. CREST is too slow for large peptides (CsA CHCl3 took 11+ days). **What data do we actually have for each residue-count partition** — permeability labels, CREMP CHCl3 ensembles, and fast vacuum descriptors — and where is the gap?

---

## Results

| residues | CycPeptMPDB total | PAMPA labels | vacuum descriptors | CREMP CHCl3 ensembles |
|---|---|---|---|---|
| 4 | 55 | 55 | 55 | 0 |
| 5 | 69 | 69 | 69 | 0 |
| 6 | 2077 | 2077 | 2077 | 1153 |
| 7 | 2070 | 2070 | 2070 | 1292 |
| 8 | 111 | 111 | 111 | 0 |
| 9 | 387 | 387 | 387 | 0 |
| 10 | 1711 | 1711 | 1711 | 12 |
| 11 | 219 | 219 | 219 | 0 |
| 12 | 213 | 213 | 213 | 0 |
| 13 | 167 | 167 | 167 | 0 |
| 14 | 120 | 120 | 119 | 0 |
| 15 | 26 | 26 | 26 | 0 |

- **≥9-mer (chameleonic regime):** 2,843 labeled compounds, **12** CREMP ensembles
- **≥11-mer (CsA regime):** 745 labeled compounds, **0** CREMP ensembles

---

## Findings

### 1. Permeability labels are abundant at every size
All 7,298 CycPeptMPDB compounds have PAMPA. The chameleonic regime is well-populated with **labels**: 2,843 compounds ≥9-mer, 745 ≥11-mer. The label side is not the bottleneck.

### 2. Vacuum RDKit dynamic descriptors already cover everything
`aq_psa3d` / `mem_psa3d` / `delta_psa3d` (Chem 269 vacuum ETKDGv3) exist for **all 7,298 compounds, all sizes** — including every large peptide. We are **not blocked on dynamic descriptors for the chameleonic regime** — a fast (if approximate) version already exists.

### 3. CREMP ensembles cover only the small regime
CREMP CHCl3 ensembles: 6-mer (1,153), 7-mer (1,292), 10-mer (12), **nothing else**. The chameleonic regime (≥9-mer) has essentially **zero** solvent-grounded ensembles. The exact regime where chameleonicity matters is where CREMP is silent.

---

## Strategic Implication — the project is not CREST-blocked

The data, by regime:

| Regime | Labels | CREMP CHCl3 | Vacuum desc | Status |
|---|---|---|---|---|
| Small (4–7mer) | ✓ (4,271) | ✓ (2,445) | ✓ | **Fully equipped** — solvent-grounded + fast |
| Large (≥9mer) | ✓ (2,843) | ✗ (12) | ✓ | **Vacuum-only** for dynamics |

**Key realization:** we have been treating CREST as the descriptor engine and getting stuck on large peptides. But fast **vacuum descriptors already exist for the entire dataset, all sizes.** A Phase-1 model can train *today* on vacuum dynamic + 2D static descriptors across all 7,298 compounds, and test the partition hypothesis immediately. CREST was never the path to high-throughput; it is the validation microscope (confirmed by the 2026-06-05/09 CsA experiments).

**How good are the vacuum descriptors?** The three-way comparison (`2026-06-09_csa_threeway_method_comparison.md`) quantified it: vacuum captures the **closed** state well (PSA 95 ≈ experimental 96) but **over-opens** the open state (180 vs 137), inflating ΔPSA ~2×. If that bias is *systematic* across compounds, the descriptors still rank — good enough for a first model, with the open state as the known weak point.

---

## Revised Plan

1. **Phase 1 (now, no new compute):** train on vacuum dynamic + 2D static descriptors, all 7,298 compounds. Run the partition analysis (6/7/8/9/10/11/12-mer) — which descriptor family predicts in each bin. This is the hypothesis test, unblocked.
2. **Phase 1.5 (small regime accuracy):** for 6–7mers, swap in CREMP CHCl3 descriptors (solvent-grounded) and check whether AUC improves vs vacuum — measures how much solvent grounding buys.
3. **Phase 2 (large regime accuracy, if needed):** the chameleonic regime depends on vacuum descriptors whose open state is over-extended. If the partition analysis shows the large bins need better dynamics, generate them with a **fast** method (OpenMM GBSA dual-dielectric, or MACE-OFF GPU) — **not** CREST. Targeted fix: correct only the open state (vacuum's weak point); the closed state is already accurate.

---

## Limitations

1. CREMP residue counts inferred from the period-separated `compound_id` (residue tokens) — may miscount exotic monomers, but the 6/7-mer dominance is unambiguous.
2. Vacuum descriptors use the max/min-PSA conformer proxy (single-ensemble extremes), which inflates ΔPSA (see three-way experiment). Acceptable for ranking if systematic; not for absolute physical magnitudes.
3. "Labels at every size" does not mean *balanced* — the 8-mer (111) and ≥13-mer bins are thin; partition tests there will have low power.
