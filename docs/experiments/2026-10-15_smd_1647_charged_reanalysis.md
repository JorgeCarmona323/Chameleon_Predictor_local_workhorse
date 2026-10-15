# 2026-10-15 — SMD vs implicit re-analysis with charged (−1) 1-6-4-7

**Context.** The 6-mer SMD-vs-implicit cross-check (colleague's explicit-DOPC SMD on Cmp1–4 = 2-9-9-8 / 3-12-8-12 / 6-4-4-13 / 1-6-4-7 vs our implicit CREST→CPCM-X descriptors) was originally run with **1-6-4-7 protonated by mistake** (its cysteic-acid sulfonate should be −1). 1-6-4-7 has now been re-run **deprotonated (charge −1)** — CREST GFN2 + CPCM-X, both legs, ensemble **2026-10-13** (from the 32-hit batch `scripts/crest_fe_hits_array_slurm.sh`; scored via `free_energy_hits_ortho_slurm.sh`). This entry records the re-analysis.

**Method.** Pulled the charged 1-6-4-7 ensembles + energy CSV to `results/conformers/hits_transfer/`, computed descriptors with `ensemble_descriptors.py` (`--apolar hexane`, CPCM-X-population-weighted), swapped the values into `scripts/smd_vs_implicit_stats.py` (`BASE["1-6-4-7"]`), and regenerated the package.

## Charged vs neutral 1-6-4-7

| | neutral (old) | charged −1 (new) |
|---|---|---|
| ΔG_transfer (w→hex) | +2.57 | **+7.41** |
| PSA water / hexane | 239 / 199 | 231 / **232** |
| Rgyr water / hexane | 4.74 / 4.57 | 4.70 / **4.64** |
| IMHB water / hexane | 3.7 / 5.4 | 4.9 / 7.0 |

## Results (figures + CSVs: `../../data/6mer comparison/2026-10-15_1647_charged_reanalysis/`)

- **PSA ranks permeability — unchanged, now un-asterisked.** ρ = 0.80 (n=4) / **1.00** (n=3 DOPC-comparable). 1-6-4-7 sits at the top-right of the PSA panel (least permeable, highest barrier 61). The result no longer carries the "1-6-4-7 was mischarged" caveat.
- **ΔG_transfer → +7.41** (was +2.57): the −1 sulfonate strongly disfavors the apolar phase, correctly the most-unfavorable of the set. n=4 ρ stays 0.40 **only** because of 6-4-4-13 (the non-DOPC confound); excluding it, ρ = **1.00**.
- **Key reframe (Rgyr).** Charging it did **not** enlarge Rgyr — it stayed compact (~4.7) even though explicit-water SMD has it splayed (core Rg 9.1; water n_eff ≈ 3.9, a few folded states dominate). So the Rgyr disagreement with SMD is **not a charge artifact** but a **sampling/population limit of implicit solvent** — the same under-opening seen at larger peptide size. This strengthens the "implicit under-samples the open state" thesis rather than being a one-off bug.

## Conclusion

The correct charge **strengthens** the headline (PSA/ΔG rank the DOPC-comparable permeability order perfectly, ρ=1.00, with no mischarge caveat) and **clarifies** the limitation (the Rgyr/geometry gap is implicit-solvent under-opening, not protonation). Chameleonicity/shape via implicit remains unreliable; PSA/ΔG-based permeability *ranking* holds.

## Files

- Script: `scripts/smd_vs_implicit_stats.py` (charged `BASE["1-6-4-7"]`, captions updated)
- Report: `data/6mer comparison/SMD_vs_implicit_analysis.md` (caveat → resolved; panel + next-steps updated)
- Dated snapshot: `data/6mer comparison/2026-10-15_1647_charged_reanalysis/` (7 figures + 4 CSVs + report copy)

## All-hits permeability ranking (done 2026-10-15)

`scripts/hits_permeability_ranking.py` → `data/hits_permeability_ranking/2026-10-15/` (`all_hits_descriptors.csv` + 3 figures). Computed CPCM-X-weighted descriptors for all 32 ortho hits; ranked by PSA(hexane) and ΔG_transfer against the controls (**NC = 1-6-4-7**, **PC = 2-9-9-8 / 3-12-8-12 R,S** from the earlier xylene SMD runs, reference values).

**Result — the NC anchor is excellent; the PCs sit mid-pack.**
- **Negative-control end is clean and physically sensible.** 1-6-4-7 is the least permeable of all 32 (PSA 232, ΔG +7.4), and its two charged cysteic-acid siblings **1-8-9-2 (214, +8.4)** and **1-9-1-5 (220, +6.2)** flank it — the three −1 sulfonates form a tight impermeable cluster at the upper-right of the 2D space. **0 hits are less permeable than the NC.**
- **Hits span a wide permeability range**, from very permeable (4-11-1-7 PSA 88, ΔG −5.3) up to the charged cluster.
- **The PCs are not at the permeable extreme.** PC PSA band is 120–156; **18 of the remaining hits have *lower* PSA than the least-permeable PC**, and 13 fall between the PC band and the NC. So "all hits between PC and NC" holds for the NC side but not the PC side — the pipeline predicts many hits *more* permeable than 2-9-9-8 / 3-12-8-12.
- PSA and ΔG agree in direction (both put the charged cluster at the impermeable corner), with the usual rank scatter in the middle.

**Caveats:** PCs are earlier xylene runs (different provenance) used as reference; a fully consistent panel needs 2-9-9-8 / 3-12-8-12 / 6-4-4-13 re-run at scale. The NC cluster is really "the charged −1 sulfonates," so it doubles as a charge sanity check (all three correctly extreme).

## Next

- Re-run 2-9-9-8 / 3-12-8-12 / 6-4-4-13 at scale for a fully consistent all-hits panel (currently the earlier xylene runs).
- With experimental permeability for any of these hits, calibrate the PSA/ΔG axis (magnitude is a free parameter; the ranking is the signal).
