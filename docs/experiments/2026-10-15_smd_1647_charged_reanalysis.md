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

## Next

- **All-hits trend run (started 2026-10-15):** rank all 32 ortho hits by PSA(hexane)/ΔG_transfer with **1-6-4-7 = negative control** (least permeable) and **2-9-9-8 + 3-12-8-12 = positive controls** (most permeable); check the remaining hits fall between the controls.
- Re-run 2-9-9-8 / 3-12-8-12 / 6-4-4-13 at scale for a fully consistent panel (currently the earlier xylene runs).
