#!/bin/bash
#SBATCH --job-name=crest_csa_gfnff
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#
# Cyclosporin A (--compound 1) in water + chloroform + hexane using GFN-FF (the GFN force
# field) instead of GFN2-xTB. GFN-FF is ~100-1000x faster, so this large 11-mer samples far
# more cheaply — a throughput experiment vs. the (slow, cis-failing) GFN2 implicit-solvent runs.
# Writes results/runs/run_<ts>_1_CsA/{water,chcl3,hexane}/ (source tag = "CREST GFN-FF ALPB").
#
# NOTE: this changes the LEVEL OF THEORY, not the solvation model — CREST still uses ALPB
# implicit solvent. The prior CsA validation showed A1 (cis-amide) fails under implicit solvent
# regardless of method, so GFN-FF here is expected to help throughput, not the cis population.

set -euo pipefail
REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs

source scripts/env.sh                     # conda env + CPCM-X-enabled xtb (single source of truth)

echo "===== CREST GFN-FF | CsA (compound 1) | water/chcl3/hexane | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)   crest: $(which crest)"

python scripts/crest_v3.2.py \
    --compound 1 \
    --threads 20 \
    --outdir results \
    --solvents water,chcl3,hexane \
    --gfn ff

echo "===== Done CsA GFN-FF | $(date) ====="
