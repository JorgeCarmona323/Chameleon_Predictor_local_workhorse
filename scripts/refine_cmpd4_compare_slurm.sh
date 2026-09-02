#!/bin/bash
#SBATCH --job-name=refine_cmpd4
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=all
#SBATCH --array=0-1
# (no #SBATCH --time -- r2SCAN-3c refinement is CPU-heavy.)
#
# Pre/post CENSO r2SCAN-3c + CPCM refinement of cmpd4's CHLOROFORM ensemble, from BOTH samplers,
# to compare GOAT-vs-CREST *after* refinement:
#   task 0 = CREST ensemble  -> <crest leg>/refined/
#   task 1 = GOAT  ensemble  -> results/goat/cmpd4/chloroform/refined/
# Serial ORCA (omp=1) + parallel across conformers, matching smd_hexpep_slurm.sh (no MPI).
#
# PREREQ: verify CENSO 3.0.8 config first -- `conda run -n orca python scripts/refine_engine.py --probe`.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs
JOBS="${SLURM_NTASKS:-20}"   # CENSO (Dask) counts SLURM TASKS as cores -> maxcores must be <= ntasks
ORCA="$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"

source scripts/env.sh
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate orca
export PATH="$(dirname "$ORCA"):$PATH"
export OMPI_MCA_rmaps_base_oversubscribe=1               # CENSO runs mpirun -np (omp-min) ORCA jobs under 1 SLURM task
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$(dirname "$ORCA"):${LD_LIBRARY_PATH:-}"
[ -x "$ORCA" ] || { echo "ERROR: ORCA not executable at $ORCA" >&2; exit 1; }
command -v censo >/dev/null || { echo "ERROR: censo not on PATH" >&2; exit 1; }
# CENSO's RRHO/GFN parts want an xtb; prefer the orca env's own, else the crest env's.
XTBEXE="$(command -v xtb 2>/dev/null || true)"; [ -x "$XTBEXE" ] || XTBEXE="$HOME/miniconda3/envs/chameleon_crest212/bin/xtb"
export XTBEXE; export PATH="$(dirname "$XTBEXE"):$PATH"

# locate cmpd4 CREST chloroform ensemble (idx 20 pattern OR Fairlie_4 run dir)
crest_ens() {
    for d in results/runs/run_*_Fairlie_4/chloroform results/runs/run_*_Fairlie_4/chcl3; do
        [ -f "$d/ensemble.xyz" ] && { echo "$d/ensemble.xyz"; return; }
    done
    ls -t results/runs/run_*/chcl3/ensemble.xyz results/runs/run_*/chloroform/ensemble.xyz 2>/dev/null \
        | while read f; do rn=$(basename "$(dirname "$(dirname "$f")")"); [[ "$rn" == *Fairlie_4* ]] && { echo "$f"; return; }; done
}

case "${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-1}" in
  0) SRC="crest"; ENS="$(crest_ens)" ;;
  1) SRC="goat";  ENS="results/goat/cmpd4/chloroform/ensemble.xyz" ;;
esac
[ -n "${ENS:-}" ] && [ -f "$ENS" ] || { echo "ERROR: cmpd4 $SRC ensemble not found (ENS='$ENS')" >&2; exit 1; }

echo "===== refine cmpd4 [$SRC] | $(date) ====="
echo "  ensemble = $ENS"
echo "  XTBEXE   = $XTBEXE"
python scripts/refine_engine.py --ensemble "$ENS" --solvent chloroform \
    --orca "$ORCA" --maxcores "$JOBS" --charge 0 --ewin 6.0
echo "===== Done cmpd4 [$SRC] -> $(dirname "$ENS")/refined/ | $(date) ====="
