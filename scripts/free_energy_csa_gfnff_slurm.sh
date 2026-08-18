#!/bin/bash
#SBATCH --job-name=fe_csa_gfnff_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=24G
#SBATCH --partition=all
#
# Task 3 — CPCM-X ΔG_transfer for Cyclosporin A (GFN-FF ensembles): water + hexane.
# CsA is ~196 atoms, so CPCM-X single-points are much slower than HexPep; --ewin 8 pre-trims
# each ensemble to its low-energy window before scoring. Geometry is GFN-FF; the CPCM-X
# single-points are GFN2-level solvation on those geometries (no re-optimization).
#
# Auto-discovers the CsA GFN-FF ensembles by their metadata (source "GFN-FF" + name
# "Cyclosporin"), so it works whether they sit in results/runs/run_*_1_CsA/ or a renamed
# results/conformers/ folder. Override with:  CSA_DIR="<dir>" sbatch scripts/free_energy_csa_gfnff_slurm.sh

set -euo pipefail
REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh                    # conda env + CPCM-X-enabled xtb-dist
export OMP_NUM_THREADS=1                  # 1 thread/xtb -> $JOBS single-points in parallel
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# --- locate the CsA GFN-FF ensembles ---------------------------------------------------
CSA_DIR="${CSA_DIR:-}"
if [ -z "$CSA_DIR" ]; then
    while IFS= read -r meta; do
        if grep -q "GFN-FF" "$meta" 2>/dev/null && grep -qi "cyclosporin" "$meta" 2>/dev/null; then
            CSA_DIR="$(dirname "$(dirname "$meta")")"; break
        fi
    done < <(find results -path "*/water/metadata.json" 2>/dev/null | sort -r)
fi
[ -n "$CSA_DIR" ] || { echo "ERROR: no CsA GFN-FF ensembles found. Set CSA_DIR=<dir> and re-submit." >&2; exit 1; }
echo "CsA GFN-FF dir: $CSA_DIR"

WATER="$CSA_DIR/water/ensemble.xyz"
HEXANE="$CSA_DIR/hexane/ensemble.xyz"
CHCL3="$CSA_DIR/chcl3/ensemble.xyz"
for p in "$WATER" "$HEXANE" "$CHCL3"; do
    [ -f "$p" ] || { echo "MISSING ensemble: $p" >&2; exit 1; }
done
echo "water      = $WATER"
echo "chloroform = $CHCL3"
echo "hexane     = $HEXANE"

# both partition legs: water->chloroform AND water->hexane (CPCM-X single-points are cheap)
python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    --leg "water=$WATER" \
    --leg "chloroform=$CHCL3" \
    --leg "hexane=$HEXANE" \
    --out "results/free_energy/csa_gfnff_cpcmx.csv"

echo
echo "Done. ΔG_transfer(water->chloroform) and (water->hexane) in results/free_energy/csa_gfnff_cpcmx.summary.csv"
