#!/bin/bash
#SBATCH --job-name=crest_wat_chcl3_hexpep
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#
# HexPep (Rezai hexamer, --compound 0) water + chloroform legs, RE-RUN under the CURRENT
# engine so all three phases (water/chcl3/hexane) come from the same crest_engine version
# with full metadata. The existing results/conformers/HexPep/{aq,mem}/full_ensemble.xyz are
# from an OLD (Apr 30) pipeline: bare ensembles, no metadata.json / sdf / crest tree, and
# unverifiable CREST settings. This job makes the set uniform for any publication rerun.
#
# Both legs share ONE cached RDKit embedding (results/embeddings/), so water and chcl3 start
# from the identical 3D structure — same convention as the hexane leg. Non-destructive:
# writes results/runs/run_<ts>_0_HexPep/{water,chcl3}/, leaving the old aq/mem tree untouched.
#
# Folder labels match the canonical tree: water -> water/ (xtb --alpb water),
# chcl3 -> chcl3/ (xtb --alpb chcl3). Pairs with crest_hexane_hexpep_slurm.sh.

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs

source scripts/env.sh                     # conda env + CPCM-X-enabled xtb (single source of truth)

echo "===== CREST water+chcl3 | HexPep (compound 0) | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python scripts/crest_v3.2.py \
    --compound 0 \
    --threads 20 \
    --outdir results \
    --solvents water,chcl3

echo "===== Done HexPep water+chcl3 | $(date) ====="
