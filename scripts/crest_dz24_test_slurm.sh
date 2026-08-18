#!/bin/bash
#SBATCH --job-name=dz24_plain
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#
# One-off test: regenerate compound 24 (Brain_6-4-4-13_diazirine) WATER leg with the current
# pipeline, to check whether the plain N=N distance constraint holds through CREST metadynamics
# now that --cinp reaches the CREST call. DIAZIRINE_FC defaults to 0.25 (the ORIGINAL force
# constant) so a PASS proves the fix was purely the --cinp wiring — nothing else needed.
# Escalate by submitting with:  DIAZIRINE_FC=1.0 sbatch scripts/crest_dz24_test_slurm.sh

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs
source scripts/env.sh
export DIAZIRINE_FC="${DIAZIRINE_FC:-0.25}"
echo "DIAZIRINE_FC=$DIAZIRINE_FC | $(date)"

python scripts/crest_v3.2.py --compound 24 --threads 20 --outdir results --solvents water

echo "Done | $(date). Verify with: python scripts/verify_diazirine_integrity.py"
