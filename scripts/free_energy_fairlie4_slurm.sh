#!/bin/bash
#SBATCH --job-name=fe_fairlie4_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# CPCM-X dG_transfer for Fairlie compound 4 (7L96), all 3 legs: water->chloroform (NMR solvent,
# membrane mimic) + water->hexane. Scores the GFN2 CREST ensembles (stage 1) and also writes the
# per-conformer 'pop' CSV the descriptor stage requires. Override the run dir with FAIRLIE4_DIR=<dir>.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# find the Fairlie_4 ensemble dir (safe: never dies on a non-matching glob)
BASE="${FAIRLIE4_DIR:-}"
if [ -z "$BASE" ]; then
  for c in results/runs/*_Fairlie_4 results/*_Fairlie_4 results/conformers/fairlie_6mer_cmpd4 results/conformers/*cmpd4*; do
    [ -d "$c" ] && ls "$c"/water/ensemble.xyz >/dev/null 2>&1 && { BASE="$c"; break; }
  done
fi
[ -n "$BASE" ] || { echo "ERROR: no Fairlie_4 ensemble dir found -- set FAIRLIE4_DIR=<dir>" >&2; exit 1; }
echo "Fairlie_4 base dir: $BASE"

python scripts/free_energy_calculator.py --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    --leg "water=$BASE/water/ensemble.xyz" \
    --leg "chloroform=$BASE/chloroform/ensemble.xyz" \
    --leg "hexane=$BASE/hexane/ensemble.xyz" \
    --out results/free_energy/fe_Fairlie_4.csv

echo "Done: results/free_energy/fe_Fairlie_4.summary.csv"
