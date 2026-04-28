# Questions for Ricardo — Tier-2 CREST Experiment Review

**Context:** First CREST run on 5 reference compounds (dual-solvent: water + CHCl3 ALPB).
Script: `scripts/tier2_crest.py` | Submission: `scripts/submit_tier2_slurm.py`
Repo: https://github.com/JorgeCarmona323/Chameleon_Predictor

---

## Sampling Parameters

**1. `--quick` / `--mquick` flags — should we drop them?**
The script currently passes both to CREST, which reduces sampling depth for speed.
CREMP ran full iMTD-GC (3.9M CPU hours for ~36k compounds — no quick flags).
Compound 0 (Hexapeptide) is in CREMP so we can cross-validate our CHCl3 run against
their published ensemble — but only if our sampling depth is comparable.
- Is `--quick`/`--mquick` acceptable for a validation run, or should we match CREMP exactly?
- At minimum, should we drop the flags for Compound 0 (Hexapeptide) where we're cross-checking?

**2. Conformer cap — is 200 enough?**
We cap at 200 conformers before Boltzmann weighting (`--max-confs 200`).
CREMP reports up to 5,743 conformers per compound (mean ~700) for 6-mers in CHCl3.
- Are we truncating meaningfully, or do high-energy conformers contribute negligibly enough that 200 is fine?
- Should we raise the cap for the 15-mers (DP-955, DP-944) which likely generate more?

**3. Energy window**
CREST default is 6 kcal/mol (consistent with Ketzel et al. 2025 for a cyclic octapeptide).
- Is 6 kcal/mol appropriate for 15-mers, or should it be wider to capture the full conformational space?

---

## Cluster Resources

**4. Max wall time per job**
Script currently requests 8h per compound (10 CREST runs total: 5 compounds × 2 solvents).
- What is the hard wall time limit on the Jinich cluster (`all` partition)?
- Is 8h realistic for a 15-mer without `--quick`? CsA (11-mer) and DP-955/DP-944 (15-mer) are the expensive ones.

**5. CPU-hour quota**
- Is there a per-user CPU-hour quota per week/month?
- 5 compounds × 2 solvents × 8 CPUs × up to 8h = up to 640 CPU-hours. Is that within limits?

---

## Parameters to Match CREMP

**6. Full parameter check**
Our setup: CREST iMTD-GC + ALPB + GFN2-xTB + 298.15 K + CHCl3 (membrane) / water (aqueous).
CREMP setup: same but CHCl3 only, full iMTD-GC (no `--quick`).
- Anything else in CREMP's methodology we should match for the CHCl3 runs?
- Charge handling: all 5 reference compounds are neutral — is `--chrg 0` implicit or should we pass it explicitly?

---

## Downstream Use (New)

**7. Top-N conformer output**
We now save the top-10 Boltzmann-weighted conformers per solvent as multi-conformer XYZ
(`results/conformers/{compound}/{solvent}/top10_boltzmann.xyz`) for downstream
docking / QM / MD input. Comment lines encode rank, Boltzmann weight, and GFN2-xTB energy.
- Is XYZ the right format for the docking engine you'd recommend, or should we also write SDF?
- For QM input (xTB single-point or ORCA): is the CREST XYZ geometry ready to use or does it need further optimization first?

---

## Open / Nice to Know

- Is there a CREST version pinned on the cluster? (`crest --version`) — CREMP used CREST 2.x; we target 3.x.
- Any macrocycle-specific CREST flags you'd recommend beyond what we have?
