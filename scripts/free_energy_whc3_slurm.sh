#!/bin/bash
#SBATCH --job-name=fe_whc3_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --time=04:00:00
#
# CPCM-X ΔG_transfer for WhC3 (White compound 3) across water/chloroform/hexane. Run AFTER
# crest_whc3_slurm.sh finishes generating the fresh 3-solvent set. Gives water->chloroform
# (literature-standard) + water->hexane. Scores the existing CREST ensembles; --ewin 8 pre-trim.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# newest WhC3 run dir (falls back to organized conformers/ if present)
BASE=""
[ -d "results/conformers/WhC3" ] && BASE="results/conformers/WhC3"
[ -z "$BASE" ] && BASE=$(ls -dt results/runs/run_*_WhC3 2>/dev/null | head -1)
[ -n "$BASE" ] || { echo "ERROR: no WhC3 directory found" >&2; exit 1; }
echo "WhC3 base dir: $BASE"

# find each leg (chloroform folder may be 'chloroform', 'chcl3', or the old 'mem')
find_leg() { local d; for d in "$@"; do [ -f "$BASE/$d/ensemble.xyz" ] && { echo "$BASE/$d/ensemble.xyz"; return; }; done; }
W=$(find_leg water aq)
C=$(find_leg chloroform chcl3 mem)
H=$(find_leg hexane)
[ -n "$W" ] && [ -n "$C" ] || { echo "ERROR: WhC3 missing water or chcl3 (water='$W' chcl3='$C')" >&2; exit 1; }
echo "  water      = $W"
echo "  chloroform = $C"
echo "  hexane     = ${H:-<none — water->chcl3 only>}"

LEGS=(--leg "water=$W" --leg "chloroform=$C")
[ -n "$H" ] && LEGS+=(--leg "hexane=$H")

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    "${LEGS[@]}" \
    --out "results/free_energy/fe_WhC3.csv"

echo
echo "Done: WhC3 ΔG_transfer -> results/free_energy/fe_WhC3.summary.csv"
