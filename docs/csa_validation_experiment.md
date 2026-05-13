# CsA Validation Experiment

## Goal

Validate that the CREST water ensemble captures the experimentally known A1 aqueous conformer of CsA, as characterized by Bhatt et al. JACS 2022 and Limbach et al. JACS 2022.

## A1 Structural Fingerprint to Verify

- Cis-amide at MeVal11−MeBmt1 (ω ≈ 0°, confirmed by Hα−Hα ROE in NMR)
- H-bond: Abu2(NH)···MeLeu10(C=O)
- H-bond: Val5(NH)···Ala7(C=O)
- Two water molecules inside the macrocycle cavity (explicit solvent only)

## Current Status

- Water ensemble: 23 conformers, dominant C1 at 46.3% Boltzmann weight
  - Location on cluster: `results/runs/run_20260503_150449_1_CsA/water/`
  - Files: `ensemble.xyz`, `ensemble.json`, `ensemble.sdf`
- CHCl3 (mem) job 259118: running as of 2026-05-12 (5+ days elapsed)
- Validation scripts written: `scripts/validate_csa_water.py`, `scripts/visualize_csa_vmd.tcl`

## Scripts

- `scripts/validate_csa_water.py` — checks cis/trans omega at MeVal11-MeBmt1, H-bond distances, reports Boltzmann-weighted fractions vs A1 fingerprint
- `scripts/visualize_csa_vmd.tcl` — VMD visualization script; load ensemble.xyz as mol 0, source this file from Tk Console; highlights NH donors (blue), C=O acceptors (red), backbone (gray)

## Next Steps

1. Check CHCl3 job status: `squeue -u j3carmona`
2. Once CHCl3 job completes, transfer both ensembles and run cross-solvent analysis
3. Run validate_csa_water.py and confirm C1 matches A1 fingerprint
4. Generate VMD visualization for PI presentation
