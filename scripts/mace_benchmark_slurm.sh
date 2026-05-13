#!/bin/bash
#SBATCH --job-name=mace_bench
#SBATCH --output=logs/mace_benchmark_%j.log
#SBATCH --error=logs/mace_benchmark_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --gres=gpu:rtx4090:1

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"

mkdir -p logs results

source "$HOME/miniconda3/etc/profile.d/conda.sh"

conda activate MACE

echo "=== Environment ==="
echo "Python:  $(which python)"
echo "PyTorch: $(python -c 'import torch; print(torch.__version__)')"
echo "CUDA:    $(python -c 'import torch; print(torch.cuda.is_available())')"
echo "GPU:     $(python -c 'import torch; print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none")')"
echo "mace:    $(python -c 'import mace; print(mace.__version__)' 2>/dev/null || echo 'not installed')"
echo "=== Starting benchmark ==="

python scripts/benchmark_mace_vs_xtb.py --model models/MACE-OFF23_medium.model

echo "=== Done ==="
echo "Results: results/mace_vs_xtb_CsA_water.csv"
echo "Plot:    results/mace_vs_xtb_CsA_water.png"
