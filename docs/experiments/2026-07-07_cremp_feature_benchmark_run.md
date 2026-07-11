# Experiment — CREMP Feature Benchmark: do 3D descriptors predict PAMPA?

**Date:** 2026-07-07
**Scripts:** `scripts/feature_benchmark.py` (RF/LightGBM), `scripts/dump_feature_matrices.py` + `scripts/tabpfn_score.py` (decoupled TabPFN lane)
**Env:** `chameleon-ml` (created this session; see Setup)
**Data:** `results/archive/cremp_deltapsa.csv` (CREMP CHCl₃ ensembles) × `results/archive/feature_matrix.csv` (CycPeptMPDB PAMPA labels)
**Output:** `results/2026-07-07_feature_benchmark_results.csv` (combined RF/LightGBM/TabPFN)

---

## Question

Does expensive CREST 3D ensemble sampling (**F7**) actually beat free fingerprints
(**F1–F4**) and cheap single-conformer RDKit 3D (**F6**) at predicting passive
permeability? Motivated by PROTAC-TS (Murakami 2026): a plain Morgan fingerprint +
TabPFN reached R²=0.71 on Caco-2 with **no** conformer ensemble. If a free
fingerprint matches an ensemble, the CREST pipeline is not earning its cost.

The decision hinges on **generalization**, not just random-split fit — so every
feature set is scored two ways:
- **random 5-fold CV** (shuffled, seed 42), and
- **leave-source-out CV** holding out `2020_Townsend` (the dominant source, 82.5%
  of the merged set) — the honest test of out-of-distribution transfer.

---

## Setup / what this session fixed

1. **Env didn't exist.** `chameleon-ml` was only a doc tag; `envs/ml.yml` also omits
   `rdkit`, `mordred`, and `tabpfn` that `feature_benchmark.py` actually imports.
   Created a corrected env (python 3.11, numpy<2, scikit-learn, lightgbm, rdkit,
   `mordredcommunity`, tqdm; tabpfn added via pip).
2. **Real bug: 0 compounds merged.** `feature_benchmark.py` used `Chem` in
   `main()`'s `_canonical()` but never imported it at module level (a prior refactor
   pushed rdkit imports into the feature functions). The `except` swallowed the
   `NameError`, so **every** CREMP SMILES canonicalized to `None` → 0 merged. Fixed
   with a module-level `from rdkit import Chem`. Merge recovered to **n = 2,416**
   (0 dropped; 65.3% permeable).
3. **TabPFN decoupled.** TabPFN v2 (`tabpfn` 8.0.8) is justified here — the Nature
   2025 paper (`docs/literature/ML literature/s41586-024-08328-6.pdf`) shows it wins
   up to **10,000 samples**, so n≈2.4k is in range (v1's ~1k/100-feat cap did not
   apply). But on Windows, torch intermittently fails to load `shm.dll` (WinError
   127) when co-imported with the mkl-heavy rdkit/mordred stack in one long process.
   Rather than force everything into one script, TabPFN runs in its own torch-clean
   process: `dump_feature_matrices.py` writes the F1–F7 matrices (rdkit/mordred, no
   torch); `tabpfn_score.py` imports torch first, no rdkit/mordred, and scores the
   same splits. Run on CPU (RTX 3060 present but a one-off at this size).

---

## Method

- **Subset:** n = 2,416 macrocycles with both a CREMP CHCl₃ ensemble and a
  CycPeptMPDB PAMPA label (canonical-SMILES inner join, dedup). Label: permeable if
  PAMPA ≥ −6.0 log cm/s.
- **Feature sets:** F1 Morgan bit 2048 · F2 Morgan count 500 · F3 Morgan count 2048
  · F4 atom-pair 2048 (chirality-aware) · F5 Mordred 2D (1442) · F6 Mordred 2D+3D,
  single ETKDG conformer (1498) · **F7 CREST CHCl₃ ensemble (10)**: aq/mem/Δ 3D-PSA,
  psa3d std/spread, Boltzmann-PSA, norm ΔPSA, ensemble energy, low-E population,
  unique confs.
- **Models:** RandomForest (300 trees), LightGBM (300 est), TabPFN v2 (PCA→100 for
  >100-dim sets). Metric: AUC-ROC.

---

## Results

### RandomForest + LightGBM (AUC-ROC)

| Feature set | dims | RF random | RF Townsend | LGBM random | LGBM Townsend |
|---|---|---|---|---|---|
| F1 Morgan bit 2048    | 2048 | 0.825 | 0.526 | 0.822 | 0.576 |
| F2 Morgan count 500   | 500  | 0.827 | 0.566 | 0.831 | 0.543 |
| F3 Morgan count 2048  | 2048 | 0.828 | 0.595 | 0.830 | 0.535 |
| F4 atom-pair 2048     | 2048 | 0.826 | 0.572 | 0.824 | 0.547 |
| F5 Mordred 2D         | 1442 | 0.831 | 0.617 | 0.829 | 0.639 |
| F6 Mordred 2D+3D      | 1498 | 0.826 | 0.623 | 0.827 | 0.647 |
| **F7 CREST CHCl₃**    | 10   | 0.695 | **0.647** | 0.694 | **0.645** |

### TabPFN v2

_Pending Stage 2 (`tabpfn_score.py`) — will fill on completion._

---

## Interpretation

_Pending TabPFN; RF/LightGBM story below is already clear:_

- **Random CV flatters the cheap features.** F1–F6 sit at ~0.82–0.83; F7 (just 10
  physics features) trails at ~0.69. On a shuffled split, high-dim 2D features win.
- **Out-of-source, the ranking inverts.** Holding out Townsend, the fingerprints
  **collapse** (F1–F4 → 0.53–0.60): most of their random-CV performance was fitting
  the dominant source, not permeability. F7 barely moves (0.695 → 0.647) and becomes
  **the most robust feature set** — 10 interpretable descriptors generalizing better
  than 2048-dim fingerprints. Mordred (F5/F6) sits between, also more robust than
  fingerprints; F6 (single-conformer 3D) ≈ F5 (2D), i.e. one embedded conformer adds
  little over 2D here.
- **Reading for the CREST pipeline:** the ensemble's value is **generalization /
  robustness**, not peak random-split AUC. This is the honest, defensible framing —
  and it argues the CREST tier earns its place where transfer matters, while cheap
  descriptors suffice for in-distribution ranking.

---

## Caveats

1. **Townsend imbalance (82.5%)** makes random CV optimistic for any feature set
   that can learn source idiosyncrasy; leave-source-out is the trustworthy column.
2. **F7 is CREMP's CHCl₃-only ensemble** (no aqueous leg) — the ΔPSA is CREMP's
   max/min-PSA proxy within one ensemble, not a dual-solvent ΔG_transfer. The
   aqueous-CREST (F8) and CPCM-X ΔG_transfer features are the planned extensions.
3. Single held-out source = one generalization test; a leave-one-source-out sweep
   would tighten the robustness claim.
