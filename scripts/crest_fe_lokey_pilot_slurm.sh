#!/bin/bash
#SBATCH --job-name=crest_fe_lokey_pilot
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-5%3
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# PILOT: 6 Lokey (2020_Townsend) 6-mer CONTROLS, all neutral, all in CREMP.
#   C titration (same 2D constitution, 2D-blinded -- only 3D stereochemistry differs):
#     lokey_3161 PAMPA -10.00 (floor) | lokey_2958 -6.49 | lokey_2868 -5.73 | lokey_3168 -4.81 (top)
#   B cross-scaffold anchors (different constitution, MW/TPSA-matched):
#     lokey_3252 PAMPA -4.00 (PC) | lokey_3418 -10.00 (NC)
#
# Purpose (n=6, extremes/titration -- a SEPARATION + monotonicity check, not a powered correlation):
#   1. shake out the full two-leg Lokey pipeline end-to-end before the ~2,435 signal-band campaign;
#   2. C: do our 3D descriptors TRACK PAMPA within one identical 2D graph (3D beats 2D)?
#   3. B: does the pipeline also order UNRELATED scaffolds correctly (generalization)?
#   4. free cross-check: our GFN2 chloroform PSA vs CREMP's cached GFN2 chloroform PSA for these 6.
#   Also settles the solvent-leg path with data: we score BOTH water->chcl3 (Path A, Rezai-aligned)
#   and water->hexane (Path B, hits/SMD-aligned) and see which separates the controls better.
#
# Method: CREST GFN2 sampling (small 6-mers -- validated hexapeptide path, matches CREMP's GFN2)
#         in water + chloroform + hexane (crest_engine default solvent_pairs). CPCM-X GFN2 single
#         points, water reference -> dG_transfer(water->chcl3) AND (water->hexane). 3D descriptors
#         are computed LOCALLY after transfer (same flow as the hits), not here.
#
# Layout written by crest_engine.generate_conformers(outdir="results/runs"):
#     results/runs/<YYYYMMDD_HHMMSS>_<short>/{water,chloroform,hexane}/ensemble.xyz (+ metadata.json)
# Output: results/free_energy/fe_lokey_<name>.csv
# Re-run safe: skips a molecule whose fe CSV already exists; reuses a complete prior run dir
#              (all 3 legs) instead of re-sampling.  Throttled %3 (<=3 concurrent).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs results/free_energy
source scripts/env.sh                      # conda env + CPCM-X-enabled xtb-dist
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"
CSV="data/lokey_pilot/lokey_pilot_input.csv"
[ -f "$CSV" ] || { echo "ERROR: $CSV not found (is data/ on the HPC?)" >&2; exit 1; }

idx="${SLURM_ARRAY_TASK_ID:?run as: sbatch scripts/crest_fe_lokey_pilot_slurm.sh}"

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

FE_OUT="results/free_energy/fe_lokey_${NAME}.csv"
if [ -s "$FE_OUT" ]; then
    echo "SKIP $NAME: $FE_OUT already exists"; exit 0
fi

echo "===== lokey_pilot[$idx] $NAME | charge=$CHG | GFN2 | water/chloroform/hexane | $(date) ====="
echo "  SMILES: $SMILES"

# reuse a complete prior run dir (all 3 legs) for this molecule, if any (latest-timestamp)
find_base() {
    local name="$1" d
    for d in $(ls -dt results/runs/*/ 2>/dev/null); do
        [ -f "${d}water/metadata.json" ] || continue
        grep -q "\"name\": \"$name\"" "${d}water/metadata.json" || continue
        if [ -f "${d}water/ensemble.xyz" ] && [ -f "${d}chloroform/ensemble.xyz" ] \
           && [ -f "${d}hexane/ensemble.xyz" ]; then echo "${d%/}"; return 0; fi
    done
    return 1
}

BASE="$(find_base "$NAME" || true)"
if [ -n "$BASE" ]; then
    echo "  reusing existing run dir: $BASE (skipping CREST)"
else
    # 1) CREST GFN2 (default), default solvent_pairs = water + chloroform + hexane
    WORKDIR_FILE="$(mktemp)"
    python - "$JOBS" "$NAME" "$CHG" "$SMILES" "$WORKDIR_FILE" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce                        # GFN_METHOD default "2" (GFN2) -- do NOT override
jobs, name, chg, smiles, wf = sys.argv[1:6]
res = ce.generate_conformers(smiles, name=name, outdir="results/runs",
        charge=int(chg), n_threads=int(jobs))     # default legs: water/chloroform/hexane
open(wf, "w").write(str(res.get("work_dir", "") or ""))
print("ok:", res.get("ok"), " charge:", res.get("charge"), " work_dir:", res.get("work_dir"))
PY
    BASE="$(cat "$WORKDIR_FILE")"; rm -f "$WORKDIR_FILE"
fi

[ -n "$BASE" ] && [ -f "$BASE/water/ensemble.xyz" ] \
    && [ -f "$BASE/chloroform/ensemble.xyz" ] && [ -f "$BASE/hexane/ensemble.xyz" ] \
    || { echo "ERROR: missing water+chloroform+hexane ensembles for $NAME" >&2; exit 1; }
echo "  CREST work dir: $BASE"

# 2) CPCM-X dG_transfer, water reference -> chcl3 and hexane legs (GFN2 single points)
python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge "$CHG" --jobs "$JOBS" \
    --leg "water=$BASE/water/ensemble.xyz" \
    --leg "chcl3=$BASE/chloroform/ensemble.xyz" \
    --leg "hexane=$BASE/hexane/ensemble.xyz" \
    --out "$FE_OUT"

echo "===== Done lokey_pilot[$idx] $NAME -> $FE_OUT | $(date) ====="
