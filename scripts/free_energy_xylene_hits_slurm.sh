#!/bin/bash
#SBATCH --job-name=fe_xylene
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=24G
#SBATCH --partition=all
#SBATCH --array=0-4
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# CPCM-X dG_transfer (water->chloroform + water->hexane) for the 5 XYLENE-linked hits
# (incl. 3-12-8-12 R/S) = the PRE-refinement energies/populations the pre/post CENSO
# comparison needs. Finds ensembles on the HPC by COMPOUND INDEX = field 4 of
# run_<date>_<time>_<idx>_<name>, latest run per solvent (same convention as
# free_energy_hits_slurm.sh), and scores chloroform too (the old hits job did water+hexane only).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# idx per crest_hexane_array_slurm.sh:
#   20=1-6-4-7   22=2-9-9-8   5=3-12-8-12 R   6=3-12-8-12 S   7=6-4-4-13   (all xylene)
INDICES=(20 22 5 6 7)
idx="${INDICES[${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-4}]}"

# find latest results/runs/run_*_<idx>_*/<solv>/ensemble.xyz ; $2.. = candidate folder names
latest_leg() {
    local want="$1"; shift
    local solv d rn
    for solv in "$@"; do
        for d in $(ls -dt results/runs/run_*/"$solv" 2>/dev/null); do
            rn=$(basename "$(dirname "$d")")
            if [ "$(echo "$rn" | cut -d_ -f4)" = "$want" ] && [ -f "$d/ensemble.xyz" ]; then
                echo "$d/ensemble.xyz"; return
            fi
        done
    done
}

WAT=$(latest_leg "$idx" water)
HEX=$(latest_leg "$idx" hexane)
CHL=$(latest_leg "$idx" chcl3 chloroform)
[ -n "$WAT" ] && [ -n "$HEX" ] || { echo "idx $idx: missing water/hexane (wat='$WAT' hex='$HEX')" >&2; exit 1; }

label=$(basename "$(dirname "$(dirname "$WAT")")" | cut -d_ -f4-)   # -> e.g. 20_1-6-4-7_xylene
charge=$(grep -oE '"charge"[^0-9-]*(-?[0-9]+)' "$(dirname "$WAT")/metadata.json" 2>/dev/null \
         | grep -oE '\-?[0-9]+$' || true); charge="${charge:-0}"

legs=( --leg "water=$WAT" --leg "hexane=$HEX" )
if [ -n "$CHL" ]; then legs+=( --leg "chloroform=$CHL" ); else echo "note: no chloroform leg for idx $idx -> water->hexane only"; fi

echo "===== fe(CPCM-X) idx=$idx label=$label charge=$charge | $(date) ====="
echo "  water  = $WAT"
echo "  hexane = $HEX"
echo "  chcl3  = ${CHL:-<none>}"

python scripts/free_energy_calculator.py --method cpcmx --ewin 8 --ref water --charge "$charge" --jobs "$JOBS" \
    "${legs[@]}" \
    --out "results/free_energy/fe_xylene_${label}.csv"

echo "===== Done | results/free_energy/fe_xylene_${label}.summary.csv | $(date) ====="
