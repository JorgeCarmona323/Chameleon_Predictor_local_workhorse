#!/bin/bash
#SBATCH --job-name=crest_roxi_gfnff
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#
# Roxithromycin (10-mer / 14-membered macrolide NMR validation anchor) in water + chloroform +
# hexane using GFN-FF, exactly like the CsA run that worked (crest_csa_gfnff_slurm.sh). The prior
# GFN2 validation run only finished the WATER leg before timing out — roxithromycin is a big,
# floppy macrolide, so GFN2 CREST is too slow to fit all 3 legs. GFN-FF samples far more cheaply.
# Policy: use GFN-FF for CsA-scale (and here roxithromycin, which hit the same wall).
# No --time cap (inherits the partition default), matching the multi-day CsA GFN-FF run.

set -euo pipefail
REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh

# resolve Roxithromycin's --compound index by name (robust to registry reordering)
IDX=$(python - <<'PY'
import importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("c32", Path("scripts/crest_v3.2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for i, c in enumerate(m.REFERENCE_COMPOUNDS):
    if c["name"] == "Roxithromycin":
        print(i); break
else:
    raise SystemExit("Roxithromycin not found in REFERENCE_COMPOUNDS")
PY
)

echo "===== CREST GFN-FF | Roxithromycin (compound $IDX) | water/chcl3/hexane | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)   crest: $(which crest)"

python scripts/crest_v3.2.py \
    --compound "$IDX" \
    --threads 20 \
    --outdir results \
    --solvents water,chcl3,hexane \
    --gfn ff

echo "===== Done Roxithromycin GFN-FF | $(date) ====="
