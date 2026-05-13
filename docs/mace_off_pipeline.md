# MACE-OFF + OpenMM Pipeline Proposal (v4.0 target)

CREST/xTB is CPU-only and does not scale to a library. This document tracks the proposed GPU-compatible replacement and its validation status.

## Why the Switch

- HexPep (6-mer): 11–15h per solvent on 20 CPU threads. Multilevel OPT = 43–49% of wall time.
- 11-mer multi-start projection: weeks per molecule per solvent. Not feasible for 50+ compounds.
- CREST/xTB has no GPU path; the architecture doesn't map to GPU.

## MACE-OFF23

Message-passing neural network (equivariant) from Csanyi group, Cambridge. Trained on SPICE dataset at ωB97M-D3(BJ)/def2-TZVPPD level — near-DFT accuracy.

- Supports H, C, N, O, F, P, S, Cl, Br — covers all macropeptide atoms
- Three sizes: S (fast, poor density), M (recommended, 0.25 kcal/mol torsion MAE), L (no improvement over M)
- Install: `pip install mace-torch`
- Single-point on CsA (196 atoms): ~10–50 ms on GPU vs ~1–2s for xTB

## Benchmark Result (2026-05-13)

Ran `benchmark_mace_vs_xtb.py` on 23 CsA water conformers vs GFN2-xTB+ALPB.

**Result: r = +0.217 — FAILED the benchmark gate (threshold: r > 0.85)**

- C1 (A-like, xTB dominant at 46.3%) ranked LAST by MACE (rank shift +21)
- MACE energy range: 32.9 kcal/mol vs xTB 5.95 kcal/mol
- Cause: gas-phase MACE vs ALPB-solvated xTB. The A conformer is stabilized by explicit water interactions that MACE (gas phase) cannot see.

**Conclusion:** MACE-OFF alone cannot replace xTB+ALPB for conformer ranking. Solvation is required.

## Proposed Pipeline (v4.0)

1. ETKDGv3 + cis/trans enumeration (keep, CPU, seconds)
2. MACE-OFF geometry optimization per cis/trans seed (GPU, ~2–5 min/conformer)
3. Short MD or metadynamics via OpenMM + MACE-OFF forces + explicit TIP3P water
4. Macrostate classification (data-driven RMSD clustering)
5. Descriptor extraction: Boltzmann PSA, HBD, ΔG(W-M) per solvent, bias index

## Key Limitations

1. **No native implicit solvation** — MACE-OFF is gas-phase only. The paper demonstrates explicit TIP3P water via OpenMM (`openmmml`). GBn2 hybrid is an option but unvalidated for this molecule class.

2. **SPICE training data has limited macrocycle coverage** — MACE-OFF23 trained on fragments/short peptides up to ~90 atoms. CsA has 196 atoms (extrapolation). No macrocycles in training set.

3. **Explicit TIP3P is expensive** — feasibility test running (`mace_tip3p_feasibility.py`). At ~50 steps/s on A100: 25M steps (50 ns) = ~140h per replica. RTX 4090 result pending.

4. **T-REMD still slow for 11-mers** — practical compromise: GFN-FF for MTD/MD sampling + MACE-OFF for energy re-ranking of selected conformers.

5. **Descriptor pipeline must be rebuilt** — crest_v3.3.py is coupled to CREST/xTB I/O. v4 needs a new pipeline rather than patching v3.3.

## Next Steps (in order)

1. Get TIP3P feasibility test step rate from `mace_tip3p_slurm.sh` job
2. If step rate < 20 steps/s: benchmark GBn2 implicit solvent as cheaper alternative
3. If GBn2 also fails: use MACE for geometry optimization only, xTB+ALPB for ranking (hybrid)
4. Validation gate before committing to v4: r > 0.85 on CsA water ensemble with chosen solvation approach

## Validation Gate

Run `benchmark_mace_vs_xtb.py` (or GBn2 variant):
- Pearson r > 0.85 → proceed with MACE as energy engine
- r 0.70–0.85 → use with caution, flag in paper
- r < 0.70 → solvation gap dominant; stay with xTB or add explicit solvent
