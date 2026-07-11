# Experiment — CsA v2 (`-notopo`) validation: does it reproduce the A1 NMR conformer?

**Date:** 2026-07-07
**Data:** v1 `data/CREST_CsA_20260512/` (23 conf) · v2 `results/conformers/CSA_v2_water/` (1,019 conf, CREST 2.12 `--noreftopo -notopo`, ALPB water)
**Scripts:** `build_ensemble_from_crest.py` (raw xyz → sdf/json), `validate_csa_water.py` (A1 fingerprint), consistent `phys_descriptors_v3` recompute of both ensembles.
**Reference:** A1 aqueous conformer, Limbach/Bhatt et al., *JACS* 2022, 144, 12602 (solution NMR + X-ray/neutron).

---

## Question

The June 5 exp-vs-CREST experiment showed CREST **v1** misses the A1 conformer two ways
(100% trans — no cis MeVal11–MeBmt1; over-collapsed water PSA). It proposed **v2
(`--noreftopo -notopo`)** to let the sampling cross the cis/trans barrier. This tests: does
v2 reproduce A1, and **which ensemble is more accurate** vs the NMR A1 conformer?

## Method

The v2 CREST run was stopped during the final tight-optimization pass, but the metadynamics
completed (all 8 MDs) and a CREGEN wrote `crest_conformers.xyz` (1,019 conformers) +
`crest.energies`. `build_ensemble_from_crest.py` rebuilt `ensemble.sdf` + `ensemble.json`
(bonds from the v1 template — atom order verified identical; Boltzmann weights from the GFN2
energies at 298.15 K). v1 and v2 were then scored with the **identical** `phys_descriptors_v3`
method so the comparison is apples-to-apples (the June-5 v1 numbers used a different PSA/Rg
method and are not directly comparable in absolute terms).

## Results

**A1 fingerprint (`validate_csa_water.py`, Boltzmann-weighted over 1,019 conf):**

| A1 criterion | v2 population |
|---|---|
| 1. Exactly one cis amide (MeVal11–MeBmt1) | **0.0%** (0 / 1019 conformers) |
| 2. Abu2 NH H-bonded | 5.2% |
| 3. Ala7 NH H-bonded | 1.1% |
| 4. Val5 NH solvent-exposed | 0.3% |
| **Full A1 match (all 4)** | **0.0%** |

**Continuous descriptors, identical method:**

| ensemble | 3D-PSA (Å²) | max PSA | Rg (Å) | IMHB | cis | |err vs A1| PSA / Rg / IMHB |
|---|---|---|---|---|---|---|
| **A1 aqueous (target)** | 137.5 | — | 6.15 | ~2 | cis | — |
| CREST **v1** (23) | 102.0 | 171.6 | 5.85 | 2.31 | trans 0% | 35.5 / 0.30 / 0.31 |
| CREST **v2** `-notopo` (1019) | 90.5 | 129.7 | 5.58 | 3.93 | trans 0% | 47.0 / 0.57 / 1.93 |

## Findings

1. **v2 does not reproduce A1.** 0 of 1,019 conformers have the cis MeVal11–MeBmt1 amide;
   full A1 match 0.0%. Still 100% trans, like v1.
2. **`-notopo` did not help accuracy — it made the match worse.** v2 is MORE collapsed than
   v1 on every axis (PSA 90.5 < 102.0, Rg 5.58 < 5.85, IMHB 3.93 > 2.31). **v1 is closer to
   the NMR A1 conformer on all three continuous descriptors.**
3. **But `-notopo` did help *sampling*** (23 → 1,019 conformers, far more converged). The
   thorough sampling revealed that CREST+ALPB genuinely *prefers* collapsed, trans, heavily
   internally-H-bonded folds. So v1's better match is best read as an **under-sampling
   artifact** (a small ensemble that stayed near the open-ish start), while v2 is the more
   truthful picture of the implicit-solvent free-energy landscape — and that landscape is wrong.

## Interpretation — the barrier is implicit solvation, not sampling

- **`-notopo` is the wrong lever for cis/trans.** It relaxes reference *topology/connectivity*
  (the log shows it acted on the 173 H/C atoms only); the amide cis/trans barrier is a
  high-barrier rotation about a partial-double-bond C–N, not a connectivity change.
- **ALPB has no cavity waters.** The two crystal waters that bridge A1's H-bond network
  (CCDC 2149649) are absent in implicit solvent, so the open/cis A1 fold is energetically
  disfavored — more sampling simply finds more of the collapsed basins ALPB prefers.
- **Conclusion:** neither v1 nor v2 reproduces A1; the dominant limitation is **implicit
  solvation**, confirming (and promoting) the "second-order limitation" flagged on June 5 to
  the primary one. Reproducing A1 needs the **explicit-water tier (OpenMM / TIP3P)**, not a
  CREST sampling flag. This retires the CsA_v2 `-notopo` hypothesis.

## Caveats

1. v2 stopped during final tight optimization; metadynamics was complete and optimization
   cannot convert trans→cis, so the all-trans / over-collapsed conclusion is robust.
2. Absolute PSA gap to A1 is method-approximate (A1 PSA computed via a different route on
   June 5); the **v1-vs-v2 direction** (same method) and the **cis result** (geometric,
   definitive) are the robust claims.
3. Implicit vs explicit solvent is a known, expected divergence for the aqueous CsA conformer.
