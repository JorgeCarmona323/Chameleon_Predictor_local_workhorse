#!/bin/bash
#SBATCH --job-name=crest_valid
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-2%2
#
# NMR pipeline-validation compounds: generate water/chloroform/hexane CREST ensembles for the
# 6-mer (Begnini_1, Begnini_2) and 10-mer (Roxithromycin) size anchors. Default solvent legs are
# water,chloroform=chcl3,hexane, exactly the three we need. Compounds are resolved BY NAME at
# runtime (not a hardcoded --compound index) so a registry reorder can't point at the wrong molecule.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh

NAMES=(Begnini_1 Begnini_2 Roxithromycin)
NAME="${NAMES[$SLURM_ARRAY_TASK_ID]}"

# name -> 0-based index into REFERENCE_COMPOUNDS
IDX=$(python - "$NAME" <<'PY'
import sys, importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("c32", Path("scripts/crest_v3.2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
want = sys.argv[1]
for i, c in enumerate(m.REFERENCE_COMPOUNDS):
    if c["name"] == want:
        print(i); break
else:
    sys.exit(f"name {want!r} not found in REFERENCE_COMPOUNDS")
PY
)

echo "===== CREST validation | $NAME = compound $IDX | 3 solvents | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python scripts/crest_v3.2.py --compound "$IDX" --threads 20 --outdir results
# (no --solvents -> default water,chloroform=chcl3,hexane)

echo "===== Done | $NAME (compound $IDX) | $(date) ====="
