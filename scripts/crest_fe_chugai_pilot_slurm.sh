#!/bin/bash
#SBATCH --job-name=crest_fe_chugai_pilot
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-4%3
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# PILOT: the 5 Chugai (2013_CHUGAI) decapeptides that overlap CREMP (IDs 743,744,745,783,792).
# All are N-methylated, NEUTRAL (charge 0 -- charge-audited, no ionizable groups).
#
# Purpose (NOT correlation -- PAMPA spans only -4.4..-5.7, n=5):
#   1. shake out the large-peptide GFN-FF pipeline end-to-end before scaling Chugai;
#   2. cross-check our GFN-FF chloroform ensemble vs CREMP's GFN2 chloroform ensemble
#      (same 5 molecules) -> validates GFN-FF for the Chugai scale-up.
#
# Method: CREST GFN-FF sampling (large peptides; CsA precedent) in CHLOROFORM ONLY -- the only
#         phase CREMP has a GFN2 ensemble for, so the only apples-to-apples benchmark. No water
#         leg => NO dG_transfer: this pilot is a pure SAMPLING benchmark (GFN-FF chloroform PSA/Rg
#         vs CREMP GFN2 chloroform PSA/Rg). CPCM-X still runs on the chloroform ensemble (self-
#         referenced) so the Boltzmann-weighted descriptors use CPCM-X populations (project rule),
#         matching CREMP's bw_psa3d. 3D descriptors are computed LOCALLY after transfer, not here.
#
# Layout written by crest_engine.generate_conformers(outdir="results/runs"):
#     results/runs/<YYYYMMDD_HHMMSS>_<short>/chloroform/ensemble.xyz (+ metadata.json)
# Output: results/free_energy/fe_chugai_<name>.csv  (chloroform-leg CPCM-X populations)
# Re-run safe: skips a molecule whose fe CSV already exists; reuses a prior run dir with a
#              chloroform ensemble instead of re-sampling.  Throttled %3 (<=3 concurrent).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs results/free_energy
source scripts/env.sh                      # conda env + CPCM-X-enabled xtb-dist
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"
CSV="data/chugai_pilot/chugai_pilot_input.csv"
[ -f "$CSV" ] || { echo "ERROR: $CSV not found (is data/ on the HPC?)" >&2; exit 1; }

idx="${SLURM_ARRAY_TASK_ID:?run as: sbatch scripts/crest_fe_chugai_pilot_slurm.sh}"

# pull name<TAB>charge<TAB>smiles for this row (header-parsed, robust to column order)
ROW="$(python - "$CSV" "$idx" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline='')))
i = int(sys.argv[2])
if i >= len(rows):
    sys.exit(f"task {i} >= {len(rows)} rows")
r = rows[i]
print(f"{r['name']}\t{r['charge']}\t{r['smiles']}")
PY
)" || { echo "ERROR: row extract failed for task $idx" >&2; exit 1; }
IFS=$'\t' read -r NAME CHG SMILES <<< "$ROW"
[ -n "$NAME" ] && [ -n "$SMILES" ] || { echo "ERROR: bad row for task $idx" >&2; exit 1; }

FE_OUT="results/free_energy/fe_chugai_${NAME}.csv"
if [ -s "$FE_OUT" ]; then
    echo "SKIP $NAME: $FE_OUT already exists"; exit 0
fi

echo "===== chugai_pilot[$idx] $NAME | charge=$CHG | GFN-FF | chloroform only | $(date) ====="
echo "  SMILES: $SMILES"

# reuse a prior run dir with a chloroform ensemble for this molecule, if any (latest-timestamp)
find_base() {
    local name="$1" d
    for d in $(ls -dt results/runs/*/ 2>/dev/null); do
        [ -f "${d}chloroform/metadata.json" ] || continue
        grep -q "\"name\": \"$name\"" "${d}chloroform/metadata.json" || continue
        if [ -f "${d}chloroform/ensemble.xyz" ]; then echo "${d%/}"; return 0; fi
    done
    return 1
}

BASE="$(find_base "$NAME" || true)"
if [ -n "$BASE" ]; then
    echo "  reusing existing run dir: $BASE (skipping CREST)"
else
    # 1) CREST GFN-FF, chloroform leg only
    WORKDIR_FILE="$(mktemp)"
    python - "$JOBS" "$NAME" "$CHG" "$SMILES" "$WORKDIR_FILE" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce
ce.GFN_METHOD = "ff"                         # GFN-FF sampling (large peptides)
jobs, name, chg, smiles, wf = sys.argv[1:6]
res = ce.generate_conformers(smiles, name=name, outdir="results/runs",
        charge=int(chg), n_threads=int(jobs),
        solvent_pairs=[("chcl3", "chloroform")])   # chloroform leg only
open(wf, "w").write(str(res.get("work_dir", "") or ""))
print("ok:", res.get("ok"), " charge:", res.get("charge"), " work_dir:", res.get("work_dir"))
PY
    BASE="$(cat "$WORKDIR_FILE")"; rm -f "$WORKDIR_FILE"
fi

[ -n "$BASE" ] && [ -f "$BASE/chloroform/ensemble.xyz" ] \
    || { echo "ERROR: missing chloroform ensemble for $NAME" >&2; exit 1; }
echo "  CREST work dir: $BASE"

# 2) CPCM-X on the chloroform ensemble, self-referenced (GFN2 single points).
#    No transfer (single phase) -- this gives the per-conformer CPCM-X Boltzmann
#    populations used to weight the 3D descriptors (project rule), matching CREMP.
python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref chcl3 --charge "$CHG" --jobs "$JOBS" \
    --leg "chcl3=$BASE/chloroform/ensemble.xyz" \
    --out "$FE_OUT"

echo "===== Done chugai_pilot[$idx] $NAME -> $FE_OUT | $(date) ====="
