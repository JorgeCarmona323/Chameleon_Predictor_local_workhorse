#!/bin/bash
#SBATCH --job-name=goat
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-7
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# GOAT (ORCA 6.1.1) conformer search @ GFN2-xTB + ALPB -- the GOAT-vs-CREST comparison layer.
# 8 jobs: 3 validation molecules (cmpd4/cmpd10/CsA) + 5 xylene hits. Seeds = frame 0 (lowest-E)
# of the existing CREST ensembles. Output -> results/goat/<name>/<solvent>/ensemble.xyz,
# parallel to (never overwriting) the CREST results. Feeds refine_engine.py + the descriptor scripts.
#
# FIRST RUN: smoke-test one task -- `sbatch --array=0 scripts/goat_slurm.sh` -- and VERIFY:
#   (a) the ORCA GOAT + ALPB input syntax below is accepted by this ORCA build,
#   (b) GOAT wrote *.finalensemble.xyz (goat_to_ensemble.py looks for that).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/goat results/slurm_logs
JOBS="${SLURM_CPUS_PER_TASK:-20}"
ORCA="$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"

# env: CENSO/ORCA python + OpenMPI (orca env) + ORCA binary. ORCA's GFN2-xTB is bundled (no external xtb).
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate orca
export PATH="$(dirname "$ORCA"):$PATH"
[ -x "$ORCA" ] || { echo "ERROR: ORCA not executable at $ORCA" >&2; exit 1; }

# "seedspec|solvent|charge|name"   seedspec = glob:<path>  OR  idx:<n> (hit run-dir index, field 4)
JOBS_SPEC=(
  "glob:results/*_Fairlie_4/chloroform/ensemble.xyz|chloroform|0|cmpd4"
  "glob:results/*_Fairlie_10/dmso/ensemble.xyz|dmso|0|cmpd10"
  "glob:results/*[Cc]yclospor*/*/ensemble.xyz|chloroform|0|CsA"
  "idx:20|chloroform|0|1-6-4-7_xylene"
  "idx:22|chloroform|0|2-9-9-8_xylene"
  "idx:5|chloroform|0|3-12-8-12_R_xylene"
  "idx:6|chloroform|0|3-12-8-12_S_xylene"
  "idx:7|chloroform|0|6-4-4-13_xylene"
)
IFS='|' read -r seedspec solvent charge name <<< "${JOBS_SPEC[${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-7}]}"

resolve_seed() {
    local spec="$1"
    if [[ "$spec" == glob:* ]]; then
        ls -t ${spec#glob:} 2>/dev/null | head -1
    elif [[ "$spec" == idx:* ]]; then
        local want="${spec#idx:}" s d rn
        for s in chcl3 chloroform; do
            for d in $(ls -dt results/runs/run_*/"$s" 2>/dev/null); do
                rn=$(basename "$(dirname "$d")")
                [ "$(echo "$rn" | cut -d_ -f4)" = "$want" ] && [ -f "$d/ensemble.xyz" ] && { echo "$d/ensemble.xyz"; return; }
            done
        done
    fi
}
SEED_ENS=$(resolve_seed "$seedspec")
[ -n "$SEED_ENS" ] && [ -f "$SEED_ENS" ] || { echo "ERROR: seed not found for '$name' (spec=$seedspec)" >&2; exit 1; }

work="results/goat/${name}/${solvent}"; mkdir -p "$work"
nat=$(head -1 "$SEED_ENS" | tr -d ' \r')
head -n $((nat + 2)) "$SEED_ENS" > "$work/seed.xyz"       # frame 0 = lowest-E CREST conformer

echo "===== GOAT $name / $solvent | charge=$charge | $(date) ====="
echo "  seed = $SEED_ENS  ($nat atoms)"

# ---- ORCA GOAT input (VERIFY syntax on this ORCA build; see header) ----
cat > "$work/goat.inp" <<INP
! GOAT GFN2-xTB ALPB($solvent)
%pal nprocs $JOBS end
%maxcore 3000
* xyzfile $charge 1 seed.xyz
INP

( cd "$work" && "$ORCA" goat.inp > goat.out 2>&1 )
grep -qiE 'ORCA TERMINATED NORMALLY' "$work/goat.out" \
    || { echo "ERROR: ORCA GOAT did not finish (tail of goat.out):" >&2; tail -25 "$work/goat.out" >&2; exit 1; }

python scripts/goat_to_ensemble.py --goat-dir "$work" --out-dir "$work" --solvent "$solvent" --charge "$charge"
echo "===== Done $name / $solvent -> $work/ensemble.xyz | $(date) ====="
