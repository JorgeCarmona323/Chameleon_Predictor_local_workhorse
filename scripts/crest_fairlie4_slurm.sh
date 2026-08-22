#!/bin/bash
#SBATCH --job-name=crest_fairlie4
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# Fairlie compound 4 (PDB 7L96): Leu cyclic hexapeptide with N-Me-Tyr, cyclo(Leu-D-Leu-Leu-
# N-Me-Tyr-D-Pro-Leu). Solution NMR in CDCl3. Generate GFN2 CREST ensembles across
# water / chloroform (= the NMR solvent, used for the RMSD-vs-NMR comparison) / hexane.
# SMILES verified visually (matches the deposited structure). Registry-free: runs from SMILES.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh
JOBS="${SLURM_CPUS_PER_TASK:-20}"

echo "===== CREST | Fairlie_4 (7L96) | GFN2 | water/chloroform/hexane | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python - "$JOBS" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce
ce.GFN_METHOD = "2"
SMILES = "CC(C)C[C@@H]1NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccc(O)cc2)N(C)C(=O)[C@H]2CCCN2C(=O)[C@H](CC(C)C)NC1=O"
res = ce.generate_conformers(
    SMILES, name="Fairlie_4", outdir="results",
    solvent_pairs=[("water", "water"), ("chcl3", "chloroform"), ("hexane", "hexane")],
    n_threads=int(sys.argv[1]))
print("ok:", res.get("ok"), " work_dir:", res.get("work_dir"))
PY

echo "===== Done | Fairlie_4 | $(date) ====="
