#!/bin/bash
#SBATCH --job-name=crest_cyclohexane
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-13
# To throttle how many array tasks run at once, append %N to the range above,
# e.g. --array=0-13%7 runs at most 7 concurrently.

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon_crest212

# 14 jobs: 5 scaffolds × {xylene, diazirine}; the 3-12 series also × {R, S}.
# Each value is a --compound index into REFERENCE_COMPOUNDS (scripts/crest_v3.2.py).
# Diazirine entries auto-get the N=N distance constraint from crest_engine.py.
COMPOUNDS=(
  20   #  0  1-6-4-7     xylene
  21   #  1  1-6-4-7     diazirine
  22   #  2  2-9-9-8     xylene
  23   #  3  2-9-9-8     diazirine
  7    #  4  6-4-4-13    xylene      (Brain_6-4-4-13; xylene linker, naphthalene residue)
  24   #  5  6-4-4-13    diazirine
  5    #  6  3-12-8-12   R xylene
  6    #  7  3-12-8-12   S xylene
  12   #  8  3-12-8-12   R diazirine
  13   #  9  3-12-8-12   S diazirine
  10   # 10  3-12-10-12  R xylene
  11   # 11  3-12-10-12  S xylene
  14   # 12  3-12-10-12  R diazirine
  15   # 13  3-12-10-12  S diazirine
)

IDX="${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}"

echo "===== CREST cyclohexane | array task=$SLURM_ARRAY_TASK_ID compound=$IDX | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

# Cyclohexane leg only (adds a cyclohexane/ ensemble alongside any existing water/chloroform).
python scripts/crest_v3.2.py \
    --compound "$IDX" \
    --threads 20 \
    --outdir results \
    --solvents cyclohexane=cyclohexane

echo "===== Done array task=$SLURM_ARRAY_TASK_ID compound=$IDX | $(date) ====="
