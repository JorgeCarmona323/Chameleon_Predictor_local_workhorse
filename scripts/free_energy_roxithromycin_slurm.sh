#!/bin/bash
#SBATCH --job-name=fe_roxi_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=24G
#SBATCH --partition=all
#
# CPCM-X ΔG_transfer for Roxithromycin across water/chloroform/hexane. Run AFTER
# crest_roxithromycin_gfnff_slurm.sh finishes. Geometry is GFN-FF; the CPCM-X single-points are
# GFN2-level solvation on those geometries (uniform scoring level with the rest of the set).
# Both partition legs: water->chloroform AND water->hexane.
#
# Default: newest Roxithromycin run (latest timestamp). Override with ROXI_DIR=<dir>.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

BASE="${ROXI_DIR:-}"
[ -z "$BASE" ] && BASE=$(ls -dt -d results/runs/run_*_Roxithromycin results/conformers/Roxithromycin 2>/dev/null | head -1 || true)
[ -n "$BASE" ] || { echo "ERROR: no Roxithromycin folder found (run crest_roxithromycin_gfnff_slurm.sh first, or set ROXI_DIR=<dir>)" >&2; exit 1; }
echo "Roxithromycin base dir (latest timestamp): $BASE"

find_leg() { local d; for d in "$@"; do [ -f "$BASE/$d/ensemble.xyz" ] && { echo "$BASE/$d/ensemble.xyz"; return; }; done; }
W=$(find_leg water)
C=$(find_leg chloroform chcl3 mem)
H=$(find_leg hexane)
[ -n "$W" ] && [ -n "$C" ] && [ -n "$H" ] || { echo "ERROR: Roxithromycin missing a leg (water='$W' chcl3='$C' hexane='$H')" >&2; exit 1; }
echo "  water      = $W"
echo "  chloroform = $C"
echo "  hexane     = $H"

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    --leg "water=$W" \
    --leg "chloroform=$C" \
    --leg "hexane=$H" \
    --out "results/free_energy/fe_Roxithromycin.csv"

echo
echo "Done: Roxithromycin ΔG_transfer(water->chloroform) and (water->hexane) -> results/free_energy/fe_Roxithromycin.summary.csv"
