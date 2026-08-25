#!/bin/bash
#SBATCH --job-name=crest_job          # <-- EDIT: shows in squeue + names the log
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
# NOTE: no #SBATCH --time -- never cap walltime on this cluster (the partition default
# already allows long jobs; a cap has killed CREST runs mid-search before).
#
# ============================================================================
#  Reusable CREST conformer-generation job.
#  Copy this file, edit the four variables in the EDIT block, then `sbatch` it.
#  Registry-free: runs straight from a SMILES (no REFERENCE_COMPOUNDS entry needed).
#
#  Output -> results/runs/<timestamp>_<NAME>/<leg>/{ensemble.xyz,ensemble.sdf,metadata.json}
#  This is stage 1 only (conformers). Score energies + descriptors separately.
# ============================================================================

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh                       # conda env + CPCM-X xtb on PATH
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# ─────────────────────────────── EDIT THIS ────────────────────────────────
NAME="MyMolecule"                            # run label -> results/runs/<ts>_<NAME>/
SMILES="O=C1N[C@H](...)...C1=O"              # the molecule (VERIFY it first!)
GFN="2"                                       # "2" = GFN2-xTB (6-8mers); "ff" = GFN-FF (large/flexible)
SOLVENTS="water,chloroform=chcl3,hexane"      # comma-sep. bare = folder==keyword;
                                              # LABEL=KEYWORD to name the folder != the xtb --alpb keyword
                                              # (e.g. dmso  ,  chloroform=chcl3). Match the NMR/assay solvent.
# ───────────────────────────────────────────────────────────────────────────

echo "===== CREST | $NAME | GFN$GFN | $SOLVENTS | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python - "$JOBS" "$NAME" "$SMILES" "$GFN" "$SOLVENTS" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce

jobs, name, smiles, gfn, solv = sys.argv[1:6]
ce.GFN_METHOD = gfn                          # switches CREST + xtb pre-opt level of theory

# parse "water,chloroform=chcl3,hexane" -> [(xtb_keyword, folder_label), ...]
pairs = []
for tok in (t.strip() for t in solv.split(",")):
    if not tok:
        continue
    label, keyword = tok.split("=", 1) if "=" in tok else (tok, tok)
    pairs.append((keyword.strip(), label.strip()))   # generate_conformers wants (solvent, label)

res = ce.generate_conformers(smiles, name=name, outdir="results/runs",
                             solvent_pairs=pairs, n_threads=int(jobs))
print("ok:", res.get("ok"), " work_dir:", res.get("work_dir"))
if not res.get("ok"):
    sys.exit("one or more solvent legs failed -- see above")
PY

echo "===== Done | $NAME | $(date) ====="
