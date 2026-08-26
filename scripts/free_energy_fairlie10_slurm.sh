#!/bin/bash
#SBATCH --job-name=fe_fairlie10_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# CPCM-X dG_transfer for Fairlie compound 10 (7L98, sanguinamide N-Me-Phe), all 3 legs:
# water->dmso (the NMR solvent, needed for the descriptor comparison to the 7L98 structure) +
# water->hexane. DMSO is a valid CPCM-X solvent. Scores the GFN2 CREST ensembles and writes the
# per-conformer 'pop' CSV the descriptor stage needs. Override the run dir with FAIRLIE10_DIR=<dir>.
#
# Note: DMSO is not a membrane mimic, so the water->dmso dG is NOT a permeability proxy -- the
# permeability number for cmpd 10 comes from the paper's PAMPA (11.0), not this dG. If CPCM-X
# rejects 'dmso' on this xtb build, drop that --leg and it still scores water->hexane.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

BASE="${FAIRLIE10_DIR:-}"
if [ -z "$BASE" ]; then
  for c in results/runs/*_Fairlie_10 results/*_Fairlie_10 results/conformers/fairlie_6mer_cmpd10 results/conformers/*cmpd10*; do
    [ -d "$c" ] && ls "$c"/water/ensemble.xyz >/dev/null 2>&1 && { BASE="$c"; break; }
  done
fi
[ -n "$BASE" ] || { echo "ERROR: no Fairlie_10 ensemble dir found -- set FAIRLIE10_DIR=<dir>" >&2; exit 1; }
echo "Fairlie_10 base dir: $BASE"

python scripts/free_energy_calculator.py --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    --leg "water=$BASE/water/ensemble.xyz" \
    --leg "dmso=$BASE/dmso/ensemble.xyz" \
    --leg "hexane=$BASE/hexane/ensemble.xyz" \
    --out results/free_energy/fe_Fairlie_10.csv

echo "Done: results/free_energy/fe_Fairlie_10.summary.csv"
