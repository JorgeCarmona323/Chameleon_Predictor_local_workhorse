#!/bin/bash
#SBATCH --job-name=crest_hex_dz_v2
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-3%3
#
# Constrained-hexane regeneration for the DOPC-series diazirines.
# The ORIGINAL diazirines (idx 12-15) were generated BEFORE the auto N=N constraint existed,
# so their ensembles let GFN2 stretch the diazirine N=N to ~1.43 A -> the +35/+9/+7 kcal/mol
# ΔG_transfer artifacts. The v2 entries (idx 16-19) are the N=N-constrained reruns; their
# water/mem ensembles already exist (constrained, Jun 16) but they lack a hexane leg.
# This generates ONLY that missing hexane leg; the current pipeline auto-applies the N=N
# distance constraint (SMARTS-detected in crest_engine.py) for these diazirine compounds.
# 21/23/24 were generated with the auto-constraint already and are unaffected.

set -euo pipefail
REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh

# v2 N=N-constrained DOPC diazirines (reruns of 12/13/14/15 respectively)
COMPOUNDS=(16 17 18 19)   # DOPCdz_R_v2, DOPCdz_S_v2, DOPCsardz_R_v2, DOPCsardz_S_v2
IDX="${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}"

echo "===== CREST hexane (N=N-constrained v2) | task=$SLURM_ARRAY_TASK_ID compound=$IDX | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python scripts/crest_v3.2.py \
    --compound "$IDX" \
    --threads 20 \
    --outdir results \
    --solvents hexane

echo "===== Done task=$SLURM_ARRAY_TASK_ID compound=$IDX | $(date) ====="
