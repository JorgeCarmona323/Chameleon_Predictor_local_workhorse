#!/bin/bash
#SBATCH --job-name=fe_solvcmp_hexpep
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --time=04:00:00
#
# Solvent-MODEL comparison on HexPep (Rezai hexamer). Score the SAME native per-phase
# ensembles (water + hexane) with each xtb solvation model and compare the resulting
# ΔG_transfer(water -> hexane). Goal: gauge how CPCM-X energies behave for our macrocycle
# system relative to the other models (later anchored to HexPep's Rezai NMR).
#
# Models run here (native xtb single-points only): ALPB, CPCM-X.
#   COSMO(-RS) and SMD move to a separate ORCA arm (more synergy with ORCA; DFT-level),
#   prepped once ORCA is installed. GBSA (legacy GB) is available in xtb too if wanted.
#
# Only water + hexane are compared: their xtb keyword is identical across all models,
# whereas chloroform is "chcl3" for ALPB/GBSA/ddCOSMO but "chloroform" for CPCM-X — and
# hexane is the transfer phase of interest regardless.
#
# Submit AFTER crest_hexane_hexpep_slurm.sh has finished (needs HexPep's hexane ensemble).

set -euo pipefail
REPO_DIR="$HOME/Chameleon_Predictor"
cd "$REPO_DIR"
mkdir -p results/free_energy results/slurm_logs
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon_crest212
export OMP_NUM_THREADS=1                 # 1 thread/xtb → run $JOBS single-points in parallel
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# --- HexPep per-phase ensembles (EDIT if the HPC paths differ) -------------------------
# hexane: auto-pick the most recent HexPep hexane run just generated
HEXANE=$(ls -dt results/runs/run_*_0_HexPep/hexane/ensemble.xyz 2>/dev/null | head -1 || true)
# water: HexPep's existing ensemble (adjust path/filename — it may be aq/full_ensemble.xyz)
WATER="results/conformers/HexPep/aq/full_ensemble.xyz"

for p in "$HEXANE" "$WATER"; do
    [ -n "${p:-}" ] && [ -f "$p" ] || { echo "MISSING ensemble: '${p:-<empty>}' — fix the path in this script"; exit 1; }
done
echo "water  = $WATER"
echo "hexane = $HEXANE"

# --- score water + hexane with each xtb solvation model --------------------------------
for METHOD in alpb cpcmx; do
    echo "=== HexPep | model=$METHOD | $(date) ==="
    python scripts/free_energy_calculator.py \
        --method "$METHOD" --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
        --leg "water=$WATER" \
        --leg "hexane=$HEXANE" \
        --out "results/free_energy/hexpep_solvcmp_${METHOD}.csv"
done

echo
echo "Done. Compare ΔG_transfer across models:"
echo "  results/free_energy/hexpep_solvcmp_{alpb,cosmo,cpcmx}.summary.csv"
echo "NOTE (CPCM-X): if its ΔG_transfer looks off vs the others, the JSON total is likely"
echo "  gas-phase — re-run cpcmx with --cpcmx-add-gsolv and compare. That IS the open question"
echo "  this experiment is meant to answer."
