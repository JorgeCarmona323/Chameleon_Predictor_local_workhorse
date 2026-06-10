# Experiment — CsA 3D Descriptors: RDKit Vacuum vs CREST Implicit vs Experimental

**Date:** 2026-06-09
**Script:** `scripts/compare_methods_csa.py`
**Output:** `results/csa_threeway_descriptors.csv`
**Data:** RDKit vacuum (`results/archive/feature_matrix.csv` row 1) · CREST V1 (`data/CREST_CsA_20260512/`) · experimental crystal (`data/experimental_structure_references_CsA/`)

---

## Question

How do the three ways we can generate CsA conformers — vacuum RDKit, implicit-solvent CREST, and experimental crystal structures — differ in their 3D descriptors? Which method's bias lands where?

---

## Method

Each method's **open** and **closed** states:
- **RDKit vacuum / CREST V1:** max-PSA conformer (open) and min-PSA conformer (closed) of the ensemble.
- **Experimental:** A1 crystal (open, aqueous, CCDC 2149649) and C1 crystal (closed, DEKSAN, CCDC 1138505).

PSA = 3D polar SASA (rdFreeSASA). Note the open/closed definitions differ by design: the computational methods use the ensemble's PSA *extremes*; experimental uses the *actual metastable states*. That difference is the point of the comparison.

---

## Results

| descriptor | RDKit vacuum | CREST V1 implicit | Experimental |
|---|---|---|---|
| open PSA (Å²) | **179.9** | 146.2 | 137.5 |
| closed PSA (Å²) | 95.0 | **51.2** | 95.9 |
| ΔPSA (Å²) | 84.9 | 94.9 | **41.6** |
| Boltzmann PSA | — | 84.1 | — |
| open HB | 1 | 1 | 2 |
| closed HB | 4 | 2 | 4 |
| open Rg | 7.09 | 6.49 | 6.15 |
| closed Rg | 5.86 | 6.37 | 6.42 |

---

## Findings

### 1. The two computational methods have opposite biases
- **RDKit vacuum over-opens the open state but nails the closed state.** Open PSA 179.9 vs experimental 137.5 (Rg 7.09 vs 6.15 — hyperextended); but closed PSA 95.0 ≈ experimental 95.9, and closed HB 4 = experimental 4. With no solvent, the open conformer hyperextends (nothing holds it together), while the closed IMHB-stabilized conformer is the vacuum global minimum and matches the real closed crystal almost exactly.
- **CREST implicit nails the open state but over-collapses the closed state.** Open PSA 146.2 ≈ experimental 137.5; but closed PSA 51.2 ≪ experimental 95.9. Implicit water keeps the open conformer realistic, but the closed conformer collapses further than the real crystal — and with fewer H-bonds (2 vs 4), suggesting a *hydrophobic* collapse rather than the H-bond-stabilized closed fold.

### 2. Both methods overestimate ΔPSA ~2× vs experimental
Computational ΔPSA (84.9 vacuum, 94.9 CREST) is roughly **double** the experimental A1↔C1 difference (41.6). Because each method takes the PSA *extremes* of the ensemble — which may be rare, transient conformers — rather than the real metastable A1/C1 states. **The max-PSA/min-PSA proxy systematically inflates ΔPSA.**

### 3. Conformational range
RDKit vacuum spans a wide Rg (5.86–7.09); CREST is narrow (6.37–6.49); experimental A1/C1 sit at 6.15–6.42. Vacuum explores more extreme geometries in both directions.

---

## Implications

- **For the descriptor pipeline:** the absolute ΔPSA from *either* computational method is inflated ~2× relative to the real metastable-state difference. This is the strongest argument yet for **macrostate-based descriptors** (cluster the ensemble to find real basins) over the max/min-PSA proxy — the proxy picks outliers, not states.
- **For ML:** if the inflation is *systematic* across compounds, relative rankings may still hold — but absolute ΔPSA should not be read as a physical switching magnitude. Prefer normalized/relative forms.
- **Complementary biases suggest a hybrid view:** vacuum captures the closed (IMHB) state well; implicit captures the open (solvated) state well. Neither alone reproduces both endpoints. Fully reproducing the experimental two-state picture likely needs explicit-solvent sampling (Tier-2 OpenMM) — consistent with the CsA exp-vs-CREST finding (2026-06-05).

---

## Limitations

1. **Open/closed defined differently across methods** (PSA extremes vs real crystal states) — intentional, but means the comparison is "each method's most-open/closed conformer vs the real open/closed forms," not a like-for-like state match.
2. **RDKit vacuum** = Chem 269 ETKDGv3 vacuum run with the max/min-PSA proxy (`feature_matrix.csv`).
3. **CREST V1** carries the known issues (no `-notopo` → all-trans, single-start, over-collapse). CsA_v2 rerun pending.
4. **Experimental** = single static crystal structures (one open, one closed), not ensembles; HB counts use the geometric cutoff, not crystallographic assignment.
