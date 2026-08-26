#!/bin/bash
#SBATCH --job-name=fe_job             # <-- EDIT: shows in squeue + names the log
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
# NOTE: no #SBATCH --time -- never cap walltime on this cluster.
#
# ============================================================================
#  Reusable STAGE-2 energy job: CPCM-X (or ALPB) dG_transfer + per-conformer
#  populations for an EXISTING CREST ensemble. Copy, edit the EDIT block, sbatch.
#
#  Single-point only -- scores the fixed CREST geometries, no re-optimization.
#  For each leg: GFN2 + CPCM-X(solvent) single-point per conformer -> solvated G
#  -> ensemble G -> dG_transfer(REF -> each leg) + the 'pop' column the
#  descriptor stage needs. Output -> results/free_energy/fe_<NAME>.{csv,summary.csv}
# ============================================================================

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh                        # conda env + CPCM-X-enabled xtb on PATH
export OMP_NUM_THREADS=1                      # 1 thread/xtb -> $JOBS single-points in parallel
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# ─────────────────────────────── EDIT THIS ────────────────────────────────
NAME="MyMolecule"                             # names the output: fe_<NAME>.csv
RUN_DIR="results/runs/<ts>_MyMolecule"        # the stage-1 CREST run dir (has <leg>/ensemble.xyz)
METHOD="cpcmx"                                 # "cpcmx" (partition-accurate) or "alpb" (fallback)
REF="water"                                   # reference leg: dG_transfer(REF -> each other leg)
CHARGE="0"
EWIN="8"                                       # kcal/mol pre-trim before scoring (throughput)
# One entry per solvent leg: "label=<path to that leg's ensemble.xyz>".
# label = the folder name / xtb solvent keyword (water, chloroform, hexane, dmso, ...).
LEGS=(
  "water=$RUN_DIR/water/ensemble.xyz"
  "chloroform=$RUN_DIR/chloroform/ensemble.xyz"
  "hexane=$RUN_DIR/hexane/ensemble.xyz"
)
# ───────────────────────────────────────────────────────────────────────────

# sanity-check every leg exists before spending compute
leg_args=()
for l in "${LEGS[@]}"; do
    path="${l#*=}"
    [ -f "$path" ] || { echo "ERROR: missing ensemble: $path" >&2; exit 1; }
    leg_args+=(--leg "$l")
done

echo "===== FE | $NAME | $METHOD | ref=$REF | ${#LEGS[@]} legs | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python scripts/free_energy_calculator.py \
    --method "$METHOD" --ewin "$EWIN" --ref "$REF" --charge "$CHARGE" --jobs "$JOBS" \
    "${leg_args[@]}" \
    --out "results/free_energy/fe_${NAME}.csv"

echo "===== Done | results/free_energy/fe_${NAME}.summary.csv | $(date) ====="
