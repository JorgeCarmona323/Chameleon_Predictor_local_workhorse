#!/bin/bash
#SBATCH --job-name=cso_crest
#SBATCH --output=logs/cso_crest_%j.log
#SBATCH --error=logs/cso_crest_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p logs results

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon

echo "===== CsO CREST run — $(date) ====="
echo "Python: $(which python)"

# CsO is index 2 in REFERENCE_COMPOUNDS (0-indexed, after Hexapeptide and CsA)
python scripts/crest_v3.2.py --compound 2 --threads 20 --outdir results

echo "===== Done — $(date) ====="
