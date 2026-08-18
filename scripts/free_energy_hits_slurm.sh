#!/bin/bash
#SBATCH --job-name=fe_hits_cpcmx
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --array=0-13%3
#
# Task 6 — CPCM-X ΔG_transfer(water -> hexane) for the 14 project hits.
# On the HPC the water and hexane legs live in SEPARATE run dirs (generated at different
# times), so we pair them by COMPOUND INDEX = field 4 of run_<date>_<time>_<idx>_<name>,
# taking the LATEST run for each solvent. Both legs are GFN2/ALPB and phase-specific.
# One array task per hit; --ewin 8 pre-trims each ensemble before scoring.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh                    # conda env + CPCM-X-enabled xtb-dist
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# 14 hits, same indices/order as crest_hexane_array_slurm.sh.
# NOTE (2026-07-28): the diazirine integrity check proved 12/13/14 have a distorted N=N
# (1.43/1.318 A vs 1.23) in >=1 leg, AND the v2 reruns (16-19) are ALSO broken (water 1.318)
# -> v2 is NOT a valid substitute. 24 is broken too (both legs ~1.319, ΔG cancelled to a
# plausible -2.80 mirage). Trustworthy diazirines: 15, 21, 23 only. 12/13/14/24 need a
# WORKING constrained regen (stronger force constant) before their ΔG can be believed.
INDICES=(20 21 22 23 7 24 5 6 12 13 10 11 14 15)
idx="${INDICES[$SLURM_ARRAY_TASK_ID]}"

# latest results/runs/run_*_<idx>_*/<solvent>/ensemble.xyz  (idx = field 4 of the run dir)
latest_leg() {
    local want="$1" solv="$2" d rn
    for d in $(ls -dt results/runs/run_*/"$solv" 2>/dev/null); do
        rn=$(basename "$(dirname "$d")")
        if [ "$(echo "$rn" | cut -d_ -f4)" = "$want" ] && [ -f "$d/ensemble.xyz" ]; then
            echo "$d/ensemble.xyz"; return
        fi
    done
}

HEX=$(latest_leg "$idx" hexane)
WAT=$(latest_leg "$idx" water)
[ -n "$HEX" ] && [ -n "$WAT" ] || { echo "idx $idx: missing leg (hex='$HEX' wat='$WAT')" >&2; exit 1; }

# charge from the hexane run's metadata (default 0)
charge=$(grep -oE '"charge"[^0-9-]*(-?[0-9]+)' "$(dirname "$HEX")/metadata.json" 2>/dev/null \
         | grep -oE '\-?[0-9]+$' || true); charge="${charge:-0}"

label=$(basename "$(dirname "$(dirname "$HEX")")" | cut -d_ -f4-)   # run dir -> e.g. 20_1-6-4-7_xylene
echo "idx=$idx  label=$label  charge=$charge"
echo "  water  = $WAT"
echo "  hexane = $HEX"

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge "$charge" --jobs "$JOBS" \
    --leg "water=$WAT" \
    --leg "hexane=$HEX" \
    --out "results/free_energy/fe_hit_${label}.csv"

echo "Done: ΔG_transfer(water->hexane) -> results/free_energy/fe_hit_${label}.summary.csv"
