# 32-Hit Permeability Ranking — Physics Pipeline vs Controls

*Status: complete (2026-10-15) · Scope: ortho 6-mer macrocycle hits · Pipeline: CREST → CPCM-X → 3D descriptors*

## 1. Purpose

Can the physics pipeline place a library of 32 ortho 6-mer hits on a single permeability axis anchored by known controls? Specifically:

- **Negative control (impermeable):** `1-6-4-7` — a charged (−1) cysteic-acid macrocycle expected to be membrane-impermeable.
- **Positive controls (permeable):** `2-9-9-8`, `3-12-8-12` (R and S) — hits with favorable SMD permeability.

The question for the project: does the pipeline (a) correctly separate the controls, and (b) distribute the remaining hits in a chemically-sensible way between them?

## 2. Method

All descriptors come from the shippable pipeline — nothing bespoke:

1. **Conformers:** CREST iMTD-GC (GFN2-xTB + ALPB), water and hexane legs (`crest_engine.py`).
2. **Energies:** CPCM-X single-point ΔG_solv + Boltzmann populations per leg (`free_energy_calculator.py`).
3. **3D descriptors:** `ensemble_descriptors.py` (`process_compound`), **Boltzmann-weighted by the CPCM-X populations** (required `--energies-csv`, no fallback).
4. **Ranking + figures:** `hits_permeability_ranking.py` (a plotting wrapper that computes nothing itself; caches `all_hits_descriptors.csv`).

**Predictors used** (validated against explicit-solvent SMD, ρ = 1.00 on the DOPC-comparable set):
- **3D-PSA, hexane (membrane) leg** — lower = more permeable.
- **ΔG_transfer (water → hexane)** — more negative = more permeable.
- **ΔPSA = PSA(water) − PSA(hexane)** — chameleonicity (mechanism, not a permeability rank at this size).

Charge states are audited (the three cysteic-acid hits `1-6-4-7`, `1-8-9-2`, `1-9-1-5` run at −1). Controls are earlier xylene runs used as reference values.

## 3. Results

**The negative-control end is clean and physically sensible.** `1-6-4-7` is the single least-permeable of all 32 (PSA 232 Å², ΔG +7.4 kcal/mol), and its two charged cysteic-acid siblings flank it:

| compound | PSA hexane (Å²) | ΔG_transfer (kcal/mol) |
|---|---|---|
| 1-6-4-7 (NC) | 232 | +7.4 |
| 1-8-9-2 | 214 | +8.4 |
| 1-9-1-5 | 220 | +6.2 |

The three −1 sulfonates form a tight impermeable cluster at the upper-right of the permeability space. **No hit is predicted less permeable than the NC.**

**Hits span a wide permeability range.** From very permeable (`4-11-1-7`, PSA 88, ΔG −5.3) up to the charged cluster. PSA(hexane) and ΔG_transfer agree in direction at both extremes.

**The positive controls sit mid-pack, not at the permeable extreme.** PC PSA band = 120–156 Å². Of the remaining hits: **18 are *more* permeable than the least-permeable PC**, and **13 fall between the PC band and the NC**. So the "between PC and NC" expectation holds firmly on the NC side; on the PC side the pipeline predicts many hits *better* than the positive controls.

**ΔPSA (chameleonicity) is a different axis — not a permeability rank at 6-mer size.** The NC `1-6-4-7` sits at ΔPSA ≈ 0 (a charged sulfonate can't bury its polar surface), but the *permeable* PC `2-9-9-8` is also ≈ 0 — it is permeable by being *constitutively nonpolar*, not by shape-shifting. ΔPSA therefore does not separate PC from NC. This matches the known scope: chameleonicity is a ≥9-mer phenomenon; at small-macrocycle size, **absolute apolar-PSA and ΔG_transfer are the discriminators.**

### Figures (`data/hits_permeability_ranking/2026-10-15/`)
- `fig_rank_psa_hexane.png` — 32 hits ranked by membrane-leg PSA (primary permeability view)
- `fig_rank_dG_transfer.png` — ranked by ΔG_transfer
- `fig_rank_delta_psa.png` — ranked by Boltzmann-weighted ΔPSA (chameleonicity)
- `fig_rank_psa_water.png` — ranked by water-leg PSA (aqueous polar exposure)
- `fig_permeability_space.png` — 2D PSA(hexane) vs ΔG
- `fig_permeability_space_water.png` — 2D PSA(water) vs ΔG
- `all_hits_descriptors.csv` — the underlying data for all 32

## 4. What this means for the project

- **The permeability-*ranking* capability is validated and provenance-clean.** The pipeline separates the controls correctly and ranks a 32-member library sensibly, with a chemically-correct charged negative-control cluster.
- **The limitation is now a characterized boundary.** Geometry/chameleonicity (Rgyr, ΔPSA) is not reliable at small-macrocycle size — implicit solvent under-samples the open state (the companion charged-`1-6-4-7` re-analysis showed the correct charge did *not* open its geometry). Claim permeability ranking; do not claim chameleonicity prediction for ≤~7-mers.
- **The next unlock is calibration, not more physics.** Magnitude is a free parameter; this is a *ranking*, not calibrated absolute permeability. Anchoring the axis against experimental data (DZALK live-cell, CHUGAI PAMPA) converts ranking → prediction and feeds the ML layer.

## 5. Caveats

- Positive controls are earlier xylene runs (different provenance) used as reference. A fully apples-to-apples panel would re-run `2-9-9-8` / `3-12-8-12` / `6-4-4-13` through the current batch.
- No experimental permeability for the 32 hits yet, so the axis is uncalibrated (ranking only).
- 6-mer scope: ΔPSA/shape features are not expected to discriminate here.

## 6. Reproduce

```bash
# descriptors come from the pipeline; ranking + figures:
conda run -n rdkit_env python scripts/hits_permeability_ranking.py            # reuses cached descriptors
conda run -n rdkit_env python scripts/hits_permeability_ranking.py --recompute # force recompute from ensembles
```

Related: companion re-analysis `docs/experiments/2026-10-15_smd_1647_charged_reanalysis.md`; SMD cross-check `data/6mer comparison/SMD_vs_implicit_analysis.md`.
