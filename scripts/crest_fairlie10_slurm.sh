#!/bin/bash
#SBATCH --job-name=crest_fairlie10
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --partition=all
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# Fairlie compound 10 (PDB 7L98): sanguinamide A with N-Me-Phe, the most permeable analog
# (PAMPA 11.0), cyclo(Ala-N-Me-Phe-Pro-Ile-Pro-Ile-thiazole). Solution NMR in DMSO-d6.
# Generate GFN2 CREST ensembles across water / DMSO (= the NMR solvent, for the RMSD-vs-NMR
# comparison) / hexane. SMILES verified against PubChem sanguinamide A + N-methyl (connectivity
# match confirmed). Registry-free: runs from SMILES.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/slurm_logs results/runs
source scripts/env.sh
JOBS="${SLURM_CPUS_PER_TASK:-20}"

echo "===== CREST | Fairlie_10 (7L98) | GFN2 | water/dmso/hexane | $(date) ====="
echo "Node: $(hostname)   Python: $(which python)"

python - "$JOBS" <<'PY'
import sys
sys.path.insert(0, "scripts")
import crest_engine as ce
ce.GFN_METHOD = "2"
SMILES = "CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)N(C)C(=O)[C@H](C)NC(=O)c2csc(n2)[C@H]([C@@H](C)CC)NC(=O)[C@@H]2CCCN2C1=O"
res = ce.generate_conformers(
    SMILES, name="Fairlie_10", outdir="results/runs",
    solvent_pairs=[("water", "water"), ("dmso", "dmso"), ("hexane", "hexane")],
    n_threads=int(sys.argv[1]))
print("ok:", res.get("ok"), " work_dir:", res.get("work_dir"))
PY

echo "===== Done | Fairlie_10 | $(date) ====="
