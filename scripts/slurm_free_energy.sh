#!/bin/bash
#SBATCH --job-name=fe_calc
#SBATCH --output=fe_%A_%a.out
#SBATCH --array=0-0            # set to 0-(N-1) for N molecules
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#
# Free-energy / ΔG_transfer scoring of native per-phase CREST ensembles.
# One array task per molecule. Edit MOLS to point at each molecule's phase dirs;
# each dir must contain a native ensemble.xyz (water + cyclohexane, CREST-generated).
#
# xTB is Linux-only -> this runs on the HPC. Confirm the solvent keywords
# ("cyclohexane" and CPCM-X support) with `xtb --version` / `xtb --help` first.

set -euo pipefail
module load xtb || true            # or: conda activate <env-with-xtb>
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_STACKSIZE=4G

# --- one line per molecule: "<label> <water_dir> <apolar_dir>" -----------------
MOLS=(
  "DOPC_3-12-8-12_S  results/conformers/.../DOPC_3-12-8-12_S/water  results/conformers/.../DOPC_3-12-8-12_S/cyclohexane"
  # "DOPC_3-12-8-12_R  .../water  .../cyclohexane"
)

read -r LABEL WATER_DIR APOLAR_DIR <<< "${MOLS[$SLURM_ARRAY_TASK_ID]}"
echo "scoring $LABEL"

python scripts/free_energy_calculator.py \
    --method cpcmx \
    --leg "water=${WATER_DIR}/ensemble.xyz" \
    --leg "cyclohexane=${APOLAR_DIR}/ensemble.xyz" \
    --ref water \
    --charge 0 \
    --jobs "${SLURM_CPUS_PER_TASK:-8}" \
    --out "fe_${LABEL}.csv"
