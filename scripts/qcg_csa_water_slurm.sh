#!/bin/bash
#SBATCH --job-name=qcg_csa_water
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# QCG (Quantum Cluster Growth = explicit water micro-solvation, CREST) of CsA with GFN-FF.
# DIAGNOSTIC: does an EXPLICIT water shell reopen the polar surface that implicit solvent
# over-shields? CsA's chameleonic dynamic range collapsed under implicit (ΔPSA ~0 vs real ~48).
# We seed the wrongly-CLOSED implicit-water conformer, grow an explicit water shell with GFN-FF
# (the only method that scales to CsA), and produce a solvated ensemble. Then compare the
# solute's PSA in the explicit shell vs implicit:
#   opens  -> the fix is the SOLVENT MODEL, and cheap GFN-FF + explicit water is adequate at scale
#   doesn't-> try a GFN2 single-point on the clusters; if THAT opens it, the force field (GFN-FF)
#             is too crude -> need GFN2-explicit / MACE-OFF.
# Flags (reconciled with this build's `crest -qcg` help):
#   -nofix           : let the SOLUTE relax/open during growth. For WATER, CREST fixes the solute
#                      by DEFAULT -- which would lock CsA in its closed shape and defeat the test.
#                      -nofix is REQUIRED here so CsA can open. (This is the key correction.)
#   -ensemble        : generate a solvated conformer ENSEMBLE (QCG Example 2).
#   -enslvl <m>      : method for the ensemble search;  -<m> sets the growth/opt level too.
#                      GFN-FF for CsA affordability.
#   -nsolv <N>       : shell size; NSOLV=0 -> omit it so QCG auto-grows to convergence (maxsolv 150).
#   --alpb water     : implicit model applied ONLY to the final single-point energy RANKING of the
#                      ensemble (per QCG Example 2); the solvent shell is still fully explicit.
#   --mdtime <ps>    : MTD time for the ensemble search. Bigger cluster (CsA + shell) => larger
#                      mdtime for a complete ensemble (docs' tip). Default here is a pilot value.
#   --wscal 1.0      : outer-wall potential scaling (as in Example 2).
#   Tunables via env: NSOLV (0=auto-grow), QMETHOD (gfnff|gfn2), MDTIME (ps).
# See data/6mer comparison/SMD_vs_implicit_analysis.md and the chameleonicity discussion.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/qcg/csa_water
source scripts/env.sh
JOBS="${SLURM_CPUS_PER_TASK:-20}"
NSOLV="${NSOLV:-40}"          # explicit waters in the shell (tunable; 0 = let QCG auto-grow to ~conv)
QMETHOD="${QMETHOD:-gfnff}"   # gfnff (fast, scales to CsA) or gfn2
MDTIME="${MDTIME:-100}"       # MTD time (ps); increase for a more complete ensemble on this big cluster

WORK="results/qcg/csa_water"
# 1) seed: lowest-E CsA conformer from the IMPLICIT-WATER ensemble (the wrongly-closed state).
#    (Use the hexane ensemble instead for the most-shielded start: CSA_SEED=".../hexane/ensemble.xyz")
SEED="${CSA_SEED:-results/conformers/Cyclosporin A GFN_FF/water/ensemble.xyz}"
[ -f "$SEED" ] || { echo "ERROR: CsA seed ensemble not found: $SEED" >&2; exit 1; }
NAT=$(head -1 "$SEED" | tr -d ' \r')
[ "$NAT" -gt 0 ] 2>/dev/null || { echo "ERROR: could not read atom count from $SEED" >&2; exit 1; }
head -n $((NAT + 2)) "$SEED" > "$WORK/solute.xyz"          # frame 0 (lowest-energy conformer)
echo "CsA seed: $SEED  ($NAT atoms)  ->  $WORK/solute.xyz"

# 2) solvent: a single water molecule
cat > "$WORK/water.xyz" <<'XYZ'
3

O   0.00000   0.00000   0.00000
H   0.75700   0.58600   0.00000
H  -0.75700   0.58600   0.00000
XYZ

# 3) QCG: grow the explicit water shell + solvated ensemble.  CsA is neutral (charge 0).
#    --nofix so the solute can OPEN; --alpb water only ranks the ensemble; GFN-FF growth + ensemble.
NSOLV_ARG=""; [ "${NSOLV:-0}" -gt 0 ] 2>/dev/null && NSOLV_ARG="--nsolv $NSOLV"
echo "===== QCG | CsA | water | $QMETHOD | nsolv=${NSOLV:-auto} | mdtime=${MDTIME}ps | --nofix | $(date) ====="
( cd "$WORK" && crest solute.xyz --qcg water.xyz $NSOLV_ARG --nofix --ensemble --alpb water \
      --wscal 1.0 --mdtime "$MDTIME" --"$QMETHOD" --enslvl "$QMETHOD" --chrg 0 --T "$JOBS" > qcg.out 2>&1 )
rc=$?
echo "crest rc=$rc  |  log: $WORK/qcg.out"
tail -20 "$WORK/qcg.out" 2>/dev/null
[ $rc -eq 0 ] || { echo "ERROR: QCG did not finish cleanly -- check $WORK/qcg.out and VERIFY --qcg flags for this CREST build" >&2; exit 1; }

echo "===== Done | QCG CsA water | $(date) ====="
echo "Outputs (QCG Example 2): $WORK/full_ensemble.xyz (solvated ensemble), crest_best.xyz (lowest),"
echo "                         full_population.dat (cluster populations)."
echo "NEXT: strip the waters from full_ensemble.xyz -> CsA-only conformers, population-weight by"
echo "      full_population.dat, compute PSA, and compare to the implicit-water PSA (does it reopen?)."
echo "      Optional: GFN2 single-point on the clusters to separate solvent-model vs GFN-FF cause."
