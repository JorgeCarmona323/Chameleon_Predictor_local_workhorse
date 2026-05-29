#!/bin/bash
#SBATCH --job-name=dopc_6mer_crest
#SBATCH --output=logs/dopc_6mer_crest_%j.log
#SBATCH --error=logs/dopc_6mer_crest_%j.err
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

COMPOUND=$1
echo "===== 6-mer CREST run compound=$COMPOUND — $(date) ====="
echo "Python: $(which python)"

python scripts/crest_v3.2.py --compound "$COMPOUND" --threads 20 --outdir results

echo "===== Done compound=$COMPOUND — $(date) ====="
