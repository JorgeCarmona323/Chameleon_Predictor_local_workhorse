#!/bin/bash
#SBATCH --job-name=fe_hits_ortho
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --array=0-31%4
#
# SCORING-ONLY: CPCM-X ΔG_transfer(water -> hexane) for the 32 ortho 6-mer hits whose CREST
# ensembles ALREADY EXIST on disk.  (The combined crest_fe_hits_array wrapper sampled all 32
# correctly but its post-CREST heredoc crashed writing a PosixPath -> every task false-failed
# BEFORE scoring, so no fe CSV was written.  Sampling is intact; this just scores it.)
#
# Layout written by crest_engine.generate_conformers(outdir="results/runs"):
#     results/runs/<YYYYMMDD_HHMMSS>_<short>/<label>/ensemble.xyz   (label = water | hexane)
#     results/runs/<YYYYMMDD_HHMMSS>_<short>/<label>/metadata.json  (name, charge, ...)
# We index by the SAME run==1 order as the generation wrapper (so task 31 = 1-6-4-7, etc.),
# resolve each molecule's LATEST run dir by matching "name" in metadata.json (latest-timestamp
# convention), read its charge from metadata, and score both legs. Standard output path:
#     results/free_energy/fe_hit_<name>.csv
# Re-run safe: re-scoring overwrites the CSV; CREST is never re-run here.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh                     # conda env + CPCM-X-enabled xtb-dist
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"
CSV="data/hits_ortho_batch.csv"
[ -f "$CSV" ] || { echo "ERROR: $CSV not found" >&2; exit 1; }

# to-run list (run==1), CRLF-safe; index by array task id -> name (same order as generation).
mapfile -t NAMES < <(awk -F',' 'NR>1 { for(i=1;i<=NF;i++) gsub(/\r/,"",$i); if($4==1) print $2 }' "$CSV")
idx="${SLURM_ARRAY_TASK_ID:?run as: sbatch scripts/free_energy_hits_ortho_slurm.sh}"
[ "$idx" -lt "${#NAMES[@]}" ] || { echo "ERROR: task $idx >= ${#NAMES[@]} rows" >&2; exit 1; }
NAME="${NAMES[$idx]}"

# latest results/runs/<ts>_<short>/ whose water metadata names this molecule AND has both legs.
find_base() {
    local name="$1" d md
    for d in $(ls -dt results/runs/*/ 2>/dev/null); do
        md="${d}water/metadata.json"
        [ -f "$md" ] || continue
        if grep -q "\"name\": \"$name\"" "$md" \
           && [ -f "${d}water/ensemble.xyz" ] && [ -f "${d}hexane/ensemble.xyz" ]; then
            echo "${d%/}"; return 0
        fi
    done
    return 1
}

BASE="$(find_base "$NAME")" || { echo "ERROR: no complete water+hexane run dir for $NAME" >&2; exit 1; }

# charge from the water leg's metadata (default 0 if absent)
charge=$(grep -oE '"charge"[^0-9-]*(-?[0-9]+)' "$BASE/water/metadata.json" 2>/dev/null \
         | grep -oE '\-?[0-9]+$' | head -1 || true); charge="${charge:-0}"

echo "===== fe_hits_ortho[$idx] $NAME | charge=$charge | $(date) ====="
echo "  base   = $BASE"
echo "  water  = $BASE/water/ensemble.xyz"
echo "  hexane = $BASE/hexane/ensemble.xyz"

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge "$charge" --jobs "$JOBS" \
    --leg "water=$BASE/water/ensemble.xyz" \
    --leg "hexane=$BASE/hexane/ensemble.xyz" \
    --out "results/free_energy/fe_hit_${NAME}.csv"

echo "===== Done fe_hits_ortho[$idx] $NAME -> results/free_energy/fe_hit_${NAME}.csv | $(date) ====="
