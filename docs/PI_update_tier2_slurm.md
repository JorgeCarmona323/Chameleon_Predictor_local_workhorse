# Tier-2 Validation — SLURM Run Update
**Date:** April 29, 2026
**Project:** Chameleon Predictor — Cyclic Peptide Membrane Permeability
**Status:** 🟡 Running on cluster

---

## What We Are Running

We submitted **5 independent SLURM jobs** on the cluster, one per reference compound. Each job runs the full **Tier-2 CREST pipeline** (`crest_v3.1.py`) — a physics-based conformer sampling workflow that computes Boltzmann-weighted ΔPSA in two dielectric environments (water ε=80, membrane ε=4.8).

This directly implements the proposal goal: *"MMFF94 minimization at ε=78 and ε=4"*, using GFN2-xTB + ALPB solvation (a more rigorous level of theory than MMFF94).

---

## Reference Compounds

These 5 compounds were selected to span the full size and permeability landscape:

| # | Compound | Size | PAMPA | Expected |
|---|---|---|---|---|
| 0 | Hexapeptide (Rezai & Lokey 2006) | 6-mer | −6.20 | Impermeable — small, polar, no switching |
| 1 | Cyclosporin A (Witek JCTC 2016) | 11-mer | −5.90 | Permeable — gold-standard chameleonic, ΔPSA ~75 Å² |
| 2 | c\*[PSLYF] (Hickey JMedChem 2016) | 11-mer | −9.10 | Impermeable — large but no switching |
| 3 | DP-955 (CHUGAI 2013) | 15-mer | −5.20 | Permeable — largest permeable in set |
| 4 | DP-944 (CHUGAI 2013) | 15-mer | −7.00 | Impermeable — largest impermeable |

---

## The Pipeline (Per Compound)

Each job runs the same 3-step workflow **twice** — once in water (ε=80) and once in membrane (ε=4.8):

```
Step 1 — RDKit ETKDGv3
  Embed 5,000 conformers → MMFF94 optimization → RMSD filter → top 50

Step 2 — GFN2-xTB geometry optimization
  50 parallel workers (one per conformer), --alpb {solvent}
  → select lowest-energy conformer as CREST seed

Step 3 — CREST iMTD-GC
  crest {xyz} --gfn2 --alpb {solvent} -T 20 --keepdir
  → generates full conformational ensemble in each environment

Output — Boltzmann-weighted ΔPSA
  ΔPSA = PSA_boltz(water) − PSA_boltz(membrane)
  ΔHB  = HB_boltz(membrane) − HB_boltz(water)
```

This matches the **CREMP protocol** (Atz et al. 2024) exactly.

---

## Technical Notes

- **CREST version:** v2.12 (downgraded from v3.x — reproducible crashes on cluster with v3)
- **Conda environment:** `chameleon_crest212`
- **Resources per job:** 20 CPUs, 16 GB RAM, partition `all`
- **Submitter script:** `scripts/submit_tier2_slurm_updated.py`
- **Output:** `results/crest_runs/run_{timestamp}_{idx}_{short}/`
- **Logs:** `results/slurm_logs/run_{timestamp}/`

---

## What We Hope to See

| Compound | ΔPSA prediction | Reasoning |
|---|---|---|
| HexPep | ~0 Å² | Too small to adopt chameleonic conformations |
| **CsA** | **~75 Å²** | Literature NMR + MD confirm large switching; key validation |
| c\*[PSLYF] | ~0 Å² | Rigid despite size |
| DP-955 | Large positive | Permeable 15-mer — should show switching |
| DP-944 | ~0 Å² | Impermeable 15-mer — no switching |

If CREST ΔPSA correctly separates permeable from impermeable compounds (particularly CsA vs. HexPep and DP-955 vs. DP-944), the pipeline is validated for broader screening.

---

## Next Steps

1. **Collect results** from the 5 running SLURM jobs
2. **Compare CREST ΔPSA** against PAMPA labels — does Boltzmann-weighted ΔPSA predict permeability?
3. **Next week's run:** Add **1NMe3** (Bockus et al., Lokey lab 2015, ID=980 in CycPeptMPDB) as a 6th reference compound alongside hit compounds from the screening library
   - 1NMe3 is a tri-N-methylated cyclic hexapeptide — permeable via backbone N-methylation (not chameleonic switching), PAMPA = −5.52
   - Together with HexPep, it will confirm that the pipeline distinguishes *mechanism* of permeability, not just outcome
4. **Tier-1 → Tier-2 comparison:** Cross-validate CREST ΔPSA against the Tier-1 ETKDG-vacuum ΔPSA already computed for the full 7k dataset
