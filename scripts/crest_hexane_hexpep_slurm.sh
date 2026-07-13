#!/bin/bash
#SBATCH --job-name=crest_hexane_hexpep
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#
# HexPep (Rezai hexamer, --compound 0) hexane leg — the one phase we're missing for the
# solvent-model comparison. Same engine/settings as the 14-compound hexane array; hexane
# is ALPB-valid so no surrogate. Writes results/runs/run_<ts>_0_HexPep/hexane/.

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon_crest212

echo "===== CREST hexane | HexPep (compound 0) | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python scripts/crest_v3.2.py \
    --compound 0 \
    --threads 20 \
    --outdir results \
    --solvents hexane

echo "===== Done HexPep hexane | $(date) ====="
