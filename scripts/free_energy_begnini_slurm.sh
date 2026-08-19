#!/bin/bash
#SBATCH --job-name=fe_begnini_cpcmx
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --array=0-1
#
# CPCM-X ΔG_transfer for the two Begnini 6-mer NMR-validation compounds, all THREE solvents:
#   water -> chloroform   (literature-standard membrane mimic; Begnini's own paper works in CHCl3)
#   water -> hexane       (for consistency with the project hits)
# Scores the already-generated CREST ensembles (no re-search). --ewin 8 pre-trims each ensemble.
# One array task per molecule. Small rigid 6-mers -> finishes in minutes.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh                    # conda env + CPCM-X-enabled xtb
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

MOLS=(Begnini_1 Begnini_2)
MOL="${MOLS[$SLURM_ARRAY_TASK_ID]}"

# Default: the <MOL> solvent-data folder with the LATEST timestamp (newest wins — simple &
# reproducible). Override a molecule explicitly via <MOL>_DIR, e.g. Begnini_1_DIR=<path>.
ovr="${MOL}_DIR"; BASE="${!ovr:-}"
[ -z "$BASE" ] && BASE=$(ls -dt -d "results/conformers/$MOL" results/runs/run_*_"$MOL" 2>/dev/null | head -1 || true)
[ -n "$BASE" ] || { echo "ERROR: no directory found for $MOL (set ${MOL}_DIR=<dir> to override)" >&2; exit 1; }
echo "$MOL base dir (latest timestamp): $BASE"

# find each leg's ensemble (the chloroform folder may be named 'chloroform' or 'chcl3')
find_leg() { local d; for d in "$@"; do [ -f "$BASE/$d/ensemble.xyz" ] && { echo "$BASE/$d/ensemble.xyz"; return; }; done; }
W=$(find_leg water)
C=$(find_leg chloroform chcl3)
H=$(find_leg hexane)
[ -n "$W" ] && [ -n "$C" ] && [ -n "$H" ] || { echo "ERROR: $MOL missing a leg (water='$W' chcl3='$C' hexane='$H')" >&2; exit 1; }
echo "  water      = $W"
echo "  chloroform = $C"
echo "  hexane     = $H"

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    --leg "water=$W" \
    --leg "chloroform=$C" \
    --leg "hexane=$H" \
    --out "results/free_energy/fe_${MOL}.csv"

echo
echo "Done: $MOL ΔG_transfer(water->chloroform) and (water->hexane) -> results/free_energy/fe_${MOL}.summary.csv"
