#!/bin/bash
#SBATCH --job-name=crest_fe_hits
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-31%3
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# BATCH: the 32 not-yet-run ortho 6-mer hits (34 total - 2-9-9-8 & 3-12-8-12 already done).
# SMILES source: data/hits_ortho_batch.csv, built from the AUTHORITATIVE ortho set = the
# Highlight_ID-labeled rows of 34hits FINAL dwar (verified ortho = matches our done runs; the
# 34_Hit_values_extracted.csv was a bad META export and is NOT used). The three cysteic-acid
# hits (1-8-9-2, 1-9-1-5, 1-6-4-7) are DEPROTONATED -> charge -1 in the CSV (see protonation QC).
# Per molecule: CREST GFN2 (water + hexane) -> CPCM-X dG_transfer (both legs). Standard outputs.
#
# Throttled to 3 concurrent (%3). If the CsA QCG job runs at the same time that is 4 total ->
# either run QCG first, or resubmit this as --array=0-31%2 to stay <= 3 concurrent.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs results/free_energy
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"
CSV="data/hits_ortho_batch.csv"
[ -f "$CSV" ] || { echo "ERROR: $CSV not found (is data/ on the HPC? see note in the header)" >&2; exit 1; }

# to-run list (run==1), CRLF-safe; index by array task id.  fields: name|charge|smiles
mapfile -t RUN < <(awk -F',' 'NR>1 { for(i=1;i<=NF;i++) gsub(/\r/,"",$i); if($4==1) print $2"|"$3"|"$5 }' "$CSV")
echo "to-run rows in CSV: ${#RUN[@]}"
idx="${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-31%3 scripts/crest_fe_hits_array_slurm.sh}"
[ "$idx" -lt "${#RUN[@]}" ] || { echo "ERROR: task $idx >= ${#RUN[@]} run rows" >&2; exit 1; }
IFS='|' read -r NAME CHG SMILES <<< "${RUN[$idx]}"
[ -n "$NAME" ] && [ -n "$SMILES" ] || { echo "ERROR: bad CSV row for task $idx" >&2; exit 1; }
echo "===== hit[$idx] $NAME | charge=$CHG | GFN2 | water/hexane | $(date) ====="
echo "  SMILES: $SMILES"

# 1) CREST GFN2, water + hexane
WORKDIR_FILE="$(mktemp)"
python - "$JOBS" "$NAME" "$CHG" "$SMILES" "$WORKDIR_FILE" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce
ce.GFN_METHOD = "2"
jobs, name, chg, smiles, wf = sys.argv[1:6]
res = ce.generate_conformers(smiles, name=name, outdir="results/runs",
        solvent_pairs=[("water", "water"), ("hexane", "hexane")],
        charge=int(chg), n_threads=int(jobs))
open(wf, "w").write(res.get("work_dir", "") or "")
print("ok:", res.get("ok"), " charge:", res.get("charge"), " work_dir:", res.get("work_dir"))
PY
BASE="$(cat "$WORKDIR_FILE")"; rm -f "$WORKDIR_FILE"
[ -n "$BASE" ] && [ -f "$BASE/water/ensemble.xyz" ] && [ -f "$BASE/hexane/ensemble.xyz" ] \
    || { echo "ERROR: CREST produced no water+hexane ensemble for $NAME" >&2; exit 1; }
echo "  CREST work dir: $BASE"

# 2) CPCM-X dG_transfer (both legs), passing the molecule's charge
python scripts/free_energy_calculator.py --method cpcmx --ewin 8 --ref water --charge "$CHG" --jobs "$JOBS" \
    --leg "water=$BASE/water/ensemble.xyz" \
    --leg "hexane=$BASE/hexane/ensemble.xyz" \
    --out "results/free_energy/fe_hit_${NAME}.csv"

echo "===== Done hit[$idx] $NAME | $(date) ====="
