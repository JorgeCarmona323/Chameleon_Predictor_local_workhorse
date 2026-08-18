#!/bin/bash
#SBATCH --job-name=fe_whc3_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
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

# Default: the WhC3 solvent-data folder with the LATEST timestamp (newest run wins — a simple,
# predictable, reproducible rule). To pin a specific folder instead, the user sets it EXPLICITLY:
#   WHC3_DIR=results/runs/run_..._WhC3 sbatch scripts/free_energy_whc3_slurm.sh
# We deliberately do NOT auto-guess by folder contents — that trades predictability for magic.
BASE="${WHC3_DIR:-}"
[ -z "$BASE" ] && BASE=$(ls -dt -d results/runs/run_*_WhC3 results/conformers/WhC3 2>/dev/null | head -1)
[ -n "$BASE" ] || { echo "ERROR: no WhC3 folder found (run crest_whc3_slurm.sh first, or set WHC3_DIR=<dir>)" >&2; exit 1; }
echo "WhC3 base dir (latest timestamp): $BASE"

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
