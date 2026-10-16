# 2026-10-15 — 32-hit permeability ranking vs PC/NC controls

**Context.** With all 32 ortho 6-mer hits sampled (CREST) and scored (CPCM-X), rank the whole library on one permeability axis anchored by controls, and check where the hits fall. Controls: **NC = 1-6-4-7** (charged −1 cysteic-acid, expected impermeable); **PC = 2-9-9-8, 3-12-8-12 R/S** (favorable SMD permeability). Companion to [the charged-1-6-4-7 SMD re-analysis](2026-10-15_smd_1647_charged_reanalysis.md).

**Method.** Descriptors from the pipeline script `ensemble_descriptors.py` (`process_compound`), Boltzmann-weighted by CPCM-X populations (required `--energies-csv`). Ranking + figures by `scripts/hits_permeability_ranking.py` (plotting wrapper only — computes no descriptors; caches `all_hits_descriptors.csv`, `--recompute` to regenerate). Inputs from `results/conformers/hits_transfer/` (ensembles + energy CSVs). PCs are earlier xylene runs used as reference values. Figure style from `scripts/plot_style.py`.

Predictors (validated vs explicit-solvent SMD, ρ = 1.00 on the DOPC-comparable set): **PSA(hexane)** and **ΔG_transfer(water→hexane)** (lower / more-negative = more permeable); **ΔPSA = PSA(water)−PSA(hexane)** (chameleonicity — mechanism, not a rank at this size).

## Result

**NC end is clean and physically sensible.** 1-6-4-7 is the least permeable of all 32 (PSA 232, ΔG +7.4); its two charged siblings flank it — a tight impermeable cluster of the three −1 sulfonates. **No hit is predicted less permeable than the NC.**

| compound | PSA hexane (Å²) | ΔG_transfer |
|---|---|---|
| 1-6-4-7 (NC) | 232 | +7.4 |
| 1-8-9-2 | 214 | +8.4 |
| 1-9-1-5 | 220 | +6.2 |

**Hits span a wide range** (most permeable: 4-11-1-7, PSA 88, ΔG −5.3). PSA and ΔG agree in direction at both extremes.

**PCs sit mid-pack, not at the permeable extreme.** PC PSA band = 120–156. Of the rest: **18 hits more permeable than the least PC**, **13 between PC band and NC**. So "between PC and NC" holds firmly on the NC side; on the PC side the pipeline predicts many hits *better* than the positive controls.

**ΔPSA is a mechanism axis, not a permeability rank here.** NC 1-6-4-7 ≈ 0 (charged, can't bury), but permeable PC 2-9-9-8 ≈ 0 too (constitutively nonpolar, not chameleonic). ΔPSA doesn't separate PC/NC — consistent with the 6-mer scope (chameleonicity is a ≥9-mer story; absolute apolar-PSA dominates at small size). See [[project_dpsa_mechanism_scope]].

## Figures / data — `data/hits_permeability_ranking/2026-10-15/`

- `fig_rank_psa_hexane.png` — primary permeability ranking
- `fig_rank_dG_transfer.png` — ranked by ΔG_transfer
- `fig_rank_delta_psa.png` — ranked by Boltzmann-weighted ΔPSA (chameleonicity)
- `fig_rank_psa_water.png` — ranked by water-leg PSA (aqueous polar exposure)
- `fig_permeability_space.png` / `fig_permeability_space_water.png` — 2D PSA vs ΔG
- `all_hits_descriptors.csv` — underlying data for all 32

## Interpretation (project state)

- Permeability **ranking** is validated and provenance-clean; the pipeline separates controls correctly and distributes a 32-member library sensibly.
- Chameleonicity/geometry is **out of scope at 6-mer size** (implicit under-opens; ΔPSA doesn't discriminate). Claim ranking, not chameleonicity, below ~7-mers.
- Next unlock = **calibration**, not more physics: anchor the PSA/ΔG axis against experimental permeability (DZALK / CHUGAI PAMPA) → ranking becomes prediction and feeds the ML layer.

## Caveats

- PCs are earlier xylene runs (different provenance); a fully consistent panel re-runs 2-9-9-8 / 3-12-8-12 / 6-4-4-13 through the current batch.
- No experimental permeability for the 32 yet → axis uncalibrated (ranking only).

## Reproduce

```bash
conda run -n rdkit_env python scripts/hits_permeability_ranking.py             # reuse cached descriptors
conda run -n rdkit_env python scripts/hits_permeability_ranking.py --recompute # recompute from ensembles
```

Notion report: `docs/notion/Hits32_Permeability_Ranking.md`.
