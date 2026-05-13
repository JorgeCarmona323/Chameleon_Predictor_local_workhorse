#!/bin/bash
#SBATCH --job-name=mace_tip3p
#SBATCH --output=logs/mace_tip3p_%j.log
#SBATCH --error=logs/mace_tip3p_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --gres=gpu:rtx4090:1

set -euo pipefail

cd "$HOME/Chameleon_Predictor"
mkdir -p logs results

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate MACE

echo "=== Environment ==="
echo "Python:  $(which python)"
echo "PyTorch: $(python -c 'import torch; print(torch.__version__)')"
echo "GPU:     $(python -c 'import torch; print(torch.cuda.get_device_name(0))')"
echo "=== Starting feasibility test ==="

python scripts/mace_tip3p_feasibility.py \
    --model models/MACE-OFF23_medium.model \
    --steps 50000 \
    --equil-steps 10000 \
    --padding 12

echo "=== Done ==="
