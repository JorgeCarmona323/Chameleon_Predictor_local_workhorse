#!/bin/bash
#SBATCH --job-name=crest_whc3
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --time=24:00:00
#
# Remake White compound 3 (WhC3) — Lokey N-methylated cyclic hexapeptide, the rigid IMHB-locked
# permeable 6-mer (White, Nat Chem Biol 2011). The old WhC3 run (May, "mem" naming) lacked a
# hexane leg; this regenerates a clean, consistent water/chloroform/hexane set with GFN2 (small
# 6-mer, so GFN2 is fine — same footing as HexPep and Begnini). Resolves the index by name.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh

IDX=$(python - <<'PY'
import importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("c32", Path("scripts/crest_v3.2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for i, c in enumerate(m.REFERENCE_COMPOUNDS):
    if c["short"] == "WhC3" or c["name"] == "White_compd3":
        print(i); break
else:
    raise SystemExit("WhC3 not found in REFERENCE_COMPOUNDS")
PY
)

echo "===== CREST GFN2 | WhC3 (compound $IDX) | water/chcl3/hexane | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)   crest: $(which crest)"

python scripts/crest_v3.2.py \
    --compound "$IDX" \
    --threads 20 \
    --outdir results \
    --solvents water,chcl3,hexane

echo "===== Done WhC3 | $(date) ====="
