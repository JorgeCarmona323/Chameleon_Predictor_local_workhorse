#!/bin/bash
#SBATCH --job-name=feat_bench
#SBATCH --output=results/slurm_logs/benchmark_%j.out
#SBATCH --error=results/slurm_logs/benchmark_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=all

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"

mkdir -p results/slurm_logs results/runs

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon

echo "Python: $(which python)"
echo "Numpy: $(python -c 'import numpy; print(numpy.__version__)')"
echo "Torch: $(python -c 'import torch; print(torch.__version__)' 2>/dev/null || echo 'not installed')"
echo "TabPFN: $(python -c 'import tabpfn; print(tabpfn.__version__)' 2>/dev/null || echo 'not installed')"

python scripts/feature_benchmark.py

echo "Done. Results written to results/feature_benchmark_results.csv"
