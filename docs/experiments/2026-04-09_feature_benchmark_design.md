# 2026-04-09 — Feature Generation Benchmark Design

## Question

Does expensive CREST-derived 3D conformer sampling provide meaningful predictive signal over cheap RDKit-embedded 3D descriptors? And does MAPC outperform standard Morgan fingerprints for permeability prediction on macrocyclic peptides?

This benchmark is motivated by PROTAC-TS (Murakami et al., *JACS Au* 2026) which achieved R²=0.710 predicting Caco-2 permeability on 43 PROTACs using only a 500-dimensional count-based Morgan fingerprint + TabPFN — no conformer ensemble, no 3D simulation. That result raises a direct challenge to the CREST-pipeline approach: **if a cheap fingerprint nearly matches or beats an expensive conformer ensemble for permeability, we need to know before committing more compute.**

---

## Literature Basis

### AGDIFF (Wu & Zou, *JCIM* 2026, `accurate-3d-structure-prediction...pdf`)

All-atom diffusion model for 3D conformer ensemble generation of cyclic peptides, trained on
CREMP (36,198 macrocyclic peptides). Uses 2D molecular graph representation — inherently
handles NCAAs, D-amino acids, N-methylation, and arbitrary cyclization chemistries. Achieves
RMSD 0.79 Å on CREMP internal holdout, 1.04 Å / rTFD 15.27° on CREMP-CycPeptMPDB (our
dataset, n=40 external test). Higher rTFD for NCAAs attributed to lactone bonds and β-amino
acids absent from training.

**Critical caveat — single solvent training:** AGDIFF was trained exclusively on CREST CHCl3
conformer ensembles. It has no exposure to aqueous conformational behavior. Consequences:

- AGDIFF can approximate membrane conformers well (CHCl3-trained ≈ nonpolar ≈ membrane-mimetic)
- AGDIFF cannot generate valid aqueous conformers — polar group exposure in water is absent
  from its learned distribution
- ΔPSA computed from AGDIFF ensembles systematically underestimates the aqueous PSA,
  compressing the ΔPSA signal
- AGDIFF is therefore **not a substitute for dual-solvent CREST** for the chameleonic ΔPSA
  framework

**Potential role:** Fast membrane-conformer generator for 800K DEL compounds (inference is
orders of magnitude faster than CREST). Could generate F_new as an additional feature set
if AGDIFF-derived ΔPSA is benchmarked alongside CREST CHCl3 and aqueous.

**Scalability caveat:** Trained on 4–6-mers. DEL library goes to 13-mers. Generalization
to larger rings is untested and expected to degrade.

**Thesis implication:** The dual-CREST (CHCl3 + aqueous) approach is a genuine
methodological contribution precisely because existing models (including AGDIFF) are
single-solvent and cannot capture the full chameleonic ΔPSA. Aqueous CREST runs are
not optional — they are the only path to valid aqueous conformers at this scale.

---

### PROTAC-TS (Murakami et al., *JACS Au* 2026, `data-driven-design-of-protac-linkers...pdf`)

Prediction model for Caco-2 permeability of PROTACs (n=43 from PROTAC-DB 3.0). Best results from Table 1:

| Feature | Model | R² |
|---|---|---|
| Morgan count-based, 500-dim | TabPFN | **0.710** |
| Morgan count-based, 2048-dim | TabPFN | 0.661 |
| Mordred 2D+3D | TabPFN | 0.636 |
| Mordred 2D only | TabPFN | 0.617 |
| Morgan bit-based, 2048-dim | TabPFN | 0.646 |

Key observations:
- **Count-based Morgan (r=2, 500-dim) was the single best feature.** Smaller dimensionality beat larger — likely because 500-dim count captures repetitive structural motifs (PEG linkers, alkyl chains) better without curse of dimensionality on n=43.
- **3D Mordred consistently added ~2-3% R² over 2D only**, even on a tiny dataset. Geometric signal matters even when n is small.
- **Count > bit**: count-based Morgan captures multiplicity of substructure occurrences; bit-based loses this.
- **TabPFN won across all feature sets** — a Bayesian tabular prior-data fitted network designed for small datasets. Relevant because our CREMP permeability subset (n=3,258) is still small-data for high-dimensional features.
- Their 3D Mordred came from a **single RDKit ETKDG-embedded conformer**, not a conformer ensemble. This is our potential edge.

### CREMP (Grambow et al., *Scientific Data* 2024, `s41597-024-03698-y.pdf`)

The dataset underlying our CREMP benchmark. Key facts:
- 36,198 unique macrocyclic peptides; 31.3M conformers from CREST in **chloroform** (ALPB implicit solvent — chosen to approximate nonpolar membrane environment)
- 3,258 molecules with experimental passive permeability from CycPeptMPDB
- Aqueous CREST was explicitly described as computationally prohibitive at this scale — they ran 3.9M CPU hours for CHCl3 alone
- Each conformer carries xTB energy; ensembles provide entropy, free energy, and occupancy metadata

### PROTAC-TS GitHub (ycu-iil/PROTAC-TS)

Open source implementation. Uses ChemTSv2 (RNN + Monte Carlo Tree Search) for linker generation guided by the Caco-2 prediction model. Nine filtering modules. Trained on PROTAC-DB 3.0 linker library.

---

## The Fingerprint Gap: Morgan vs. MAPC

PROTAC-TS uses standard RDKit Morgan fingerprints (circular, ECFP-style). Our pipeline uses **MAPC** (MinHashed Atom-Pair Fingerprint Chiral), a fundamentally different algorithm:

| Property | Morgan (ECFP) | MAPC |
|---|---|---|
| Algorithm | Circular neighborhood hashing | MinHashed atom-pair paths |
| Chirality | Optional flag | Native |
| Output | Fixed-length bit or count vector | Fixed-length (default 2048) |
| Captures | Local chemical environment | Pairwise atom relationships across molecule |
| Strength | Repetitive local motifs | Long-range structural relationships |

For macrocyclic peptides with complex ring topology and critical chirality (D-amino acids, N-methylation), MAPC's chirality-awareness and long-range path encoding may be advantageous. But this is untested — the PROTAC-TS result shows Morgan alone is very competitive.

---

## Proposed Benchmark

### Dataset
CREMP-CycPeptMPDB permeability subset: n=3,258 molecules with experimental passive permeability.

### Feature Sets to Compare

| ID | Feature | Source | Compute cost |
|---|---|---|---|
| F1 | Morgan bit-based (r=2, 2048-dim) | RDKit | ~free |
| F2 | Morgan count-based (r=2, 500-dim) | RDKit | ~free |
| F3 | Morgan count-based (r=2, 2048-dim) | RDKit | ~free |
| F4 | MAPC (2048-dim) | mhfp library | ~free |
| F5 | Mordred 2D only | mordred | low |
| F6 | Mordred 2D+3D (single ETKDG conformer) | RDKit + mordred | low |
| F7 | CREST ensemble descriptors, CHCl3 | CREMP pipeline | **done** |
| F8 | CREST ensemble descriptors, aqueous | CREMP pipeline | **pending — compute acquired** |
| F9 | AGDIFF-generated conformer descriptors | AGDIFF (Wu & Zou 2026) | pending — membrane-only |

F7 and F8 may include: ensemble-averaged Mordred 3D, min/max-PSA conformer descriptors, ensemble entropy, Rg distribution, ΔPSA metrics.

### Models
- **TabPFN** (PROTAC-TS winner; designed for small tabular datasets)
- **LightGBM** (their second model)
- Current pipeline model (for direct comparison)

### Evaluation
- Same stratified CV splits across all feature sets (shuffle seed fixed)
- Metrics: R², Pearson r, Spearman ρ, RMSE (continuous); AUC-ROC (binary at -6.0 log cm/s threshold)
- Report mean ± SD across folds

### Three Questions This Answers

1. **MAPC vs. Morgan** — Does the chirality-aware atom-pair fingerprint outperform standard circular FP for macrocyclic permeability?
2. **RDKit 3D vs. CREST CHCl3** — Is expensive conformer sampling (39M CPU-hrs for the full CREMP set) worth it over a free ETKDG embedding for predictive performance?
3. **CREST CHCl3 vs. CREST aqueous** — Does running CREST in aqueous implicit solvent add signal over the chloroform ensemble? This becomes a free comparison once aqueous runs are done.
4. **AGDIFF vs. CREST CHCl3** — Does a fast ML-generated membrane-conformer ensemble match physics-based CREST? If yes, AGDIFF becomes a scalable conformer engine for 800K DEL compounds. If no, CREST remains necessary.
5. **AGDIFF limitation** — AGDIFF cannot generate aqueous conformers (single-solvent CHCl3 training). F9 therefore represents a membrane-only signal upper bound, not a ΔPSA estimate. This is a fundamental constraint, not a benchmark failure.

---

## Why This Is Decision-Critical

If CREST 3D >> RDKit 3D: the ensemble pipeline is justified and aqueous CREST compute is worthwhile.

If RDKit 3D ≈ CREST 3D: the expensive conformer generation adds little predictive value; future effort should go to model architecture and feature engineering, not simulation.

If Morgan >> MAPC: we should switch fingerprints for the ML pipeline.

If MAPC >> Morgan: we have evidence that chirality-aware long-range features matter for macrocyclic permeability — a publishable finding.

The answer also directly informs whether the CREST aqueous runs (now computationally accessible) are a high-priority experiment or a diminishing-returns one.

---

## Implementation Plan

1. Compute F1–F6 for all 3,258 CREMP permeability compounds (one script, ~minutes on CPU)
2. Load F7 from existing CREMP pipeline outputs
3. Train TabPFN + LightGBM on each feature set with same CV splits
4. F8 added as final comparison once aqueous CREST runs complete
5. Plot: R² heatmap (features × models) + parity plots for top-3 feature sets

---

## References

- Murakami et al., *JACS Au* 2026, 6, 1400–1410. DOI: 10.1021/jacsau.6c00033
- Grambow et al., *Scientific Data* 2024, 11, 859. DOI: 10.1038/s41597-024-03698-y
- PROTAC-TS GitHub: https://github.com/ycu-iil/PROTAC-TS
- Yu et al. (Delta PSA), bioRxiv 2026. DOI: 10.64898/2026.01.06.697862
- Hung, Venkatesan & Chang (ICoN-v1), bioRxiv 2026. DOI: 10.64898/2026.03.12.711417
