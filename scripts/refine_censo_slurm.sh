#!/bin/bash
#SBATCH --job-name=refine_xylene
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --partition=all
#SBATCH --array=0-4
# (no #SBATCH --time -- r2SCAN-3c refinement is CPU-heavy, ~days/compound.)
#
# STAGE 1.5: r2SCAN-3c + CPCM (CENSO 3.0.8, ORCA 6.1.1) refinement of the 5 xylene hits'
# <SOLVENT> leg. PRE = GFN2 ensemble; POST = this refined ensemble. NON-DESTRUCTIVE
# (writes <leg>/refined/). Set REFINE_SOLVENT to pick the leg (default chloroform).
#
# FIRST TIME: run `conda run -n orca python scripts/refine_engine.py --probe` once to confirm
# the CENSO 3.0.8 config keys / output filenames, then adjust refine_engine.py if needed.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs
JOBS="${SLURM_CPUS_PER_TASK:-20}"
ORCA="$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"
SOLVENT="${REFINE_SOLVENT:-chloroform}"

# --- env: xtb (RRHO/GFN parts) + CENSO/OpenMPI (orca env) + ORCA binary ---
source scripts/env.sh                                    # xtb-dist + XTBPATH on PATH
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate orca                                      # CENSO 3.0.8 + OpenMPI 4.1.8
export PATH="$(dirname "$ORCA"):$PATH"
[ -x "$ORCA" ] || { echo "ERROR: ORCA not executable at $ORCA" >&2; exit 1; }
command -v censo >/dev/null || { echo "ERROR: censo not on PATH after 'conda activate orca'" >&2; exit 1; }

# idx per crest_hexane_array_slurm.sh: 20=1-6-4-7 22=2-9-9-8 5=3-12-8-12R 6=3-12-8-12S 7=6-4-4-13
INDICES=(20 22 5 6 7)
idx="${INDICES[${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-4}]}"

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

case "$SOLVENT" in
  chloroform|chcl3) ENS=$(latest_leg "$idx" chcl3 chloroform) ;;
  *)                ENS=$(latest_leg "$idx" "$SOLVENT") ;;
esac
[ -n "$ENS" ] || { echo "idx $idx: no '$SOLVENT' ensemble found under results/runs/" >&2; exit 1; }
label=$(basename "$(dirname "$(dirname "$ENS")")" | cut -d_ -f4-)
chg=$(grep -oE '"charge"[^0-9-]*(-?[0-9]+)' "$(dirname "$ENS")/metadata.json" 2>/dev/null \
      | grep -oE '\-?[0-9]+$' || true); chg="${chg:-0}"

echo "===== refine idx=$idx label=$label solvent=$SOLVENT charge=$chg | $(date) ====="
echo "  ensemble = $ENS"
python scripts/refine_engine.py --ensemble "$ENS" --solvent "$SOLVENT" \
    --orca "$ORCA" --maxcores "$JOBS" --charge "$chg" --ewin 6.0
echo "===== Done idx=$idx -> $(dirname "$ENS")/refined/ | $(date) ====="
