#!/bin/bash
#SBATCH --job-name=goat
#SBATCH --output=results/slurm_logs/%x_%A_%a.out
#SBATCH --error=results/slurm_logs/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=all
#SBATCH --array=0-7%4     # %4 = at most 4 tasks running at once (throttle)
# (no #SBATCH --time -- never cap walltime on this cluster.)
#
# GOAT (ORCA 6.1.1) conformer search @ GFN2-xTB + ALPB -- the GOAT-vs-CREST comparison layer.
# 8 jobs: 3 validation molecules (cmpd4/cmpd10/CsA) + 5 xylene hits.
#
# SEEDING (FlexiSol-style): each GOAT run starts from a SMILES -> ETKDG 3D structure
# (make_goat_seed.py, same REFERENCE_COMPOUNDS SMILES CREST uses) -- NOT a CREST conformer.
# That keeps the GOAT-vs-CREST comparison independent and matches FlexiSol's protocol.
# Output -> results/goat/<name>/<solvent>/ensemble.xyz (parallel to CREST, non-destructive).
#
# FIRST RUN: smoke-test one task -- `sbatch --array=0 scripts/goat_slurm.sh` -- and VERIFY:
#   (a) ORCA accepts the `! GOAT GFN2-xTB ALPB(<solvent>)` input below,
#   (b) GOAT wrote *.finalensemble.xyz (goat_to_ensemble.py reads that).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/goat results/slurm_logs
JOBS="${SLURM_NTASKS:-20}"        # ORCA nprocs = MPI ranks = SLURM tasks (not cpus-per-task)
ORCA="$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca"
SEED_ENV="chameleon_crest212"        # env with RDKit + the registry (for make_goat_seed.py)

# "name|solvent"  (charge is derived from the SMILES by make_goat_seed.py)
JOBS_SPEC=(
  "cmpd4|chloroform"
  "cmpd10|dmso"
  "CsA|chloroform"
  "1-6-4-7_xylene|chloroform"
  "2-9-9-8_xylene|chloroform"
  "3-12-8-12_R_xylene|chloroform"
  "3-12-8-12_S_xylene|chloroform"
  "6-4-4-13_xylene|chloroform"
)
IFS='|' read -r name solvent <<< "${JOBS_SPEC[${SLURM_ARRAY_TASK_ID:?run as: sbatch --array=0-7}]}"
work="results/goat/${name}/${solvent}"; mkdir -p "$work"

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# ---- 1) SMILES -> ETKDG 3D seed (RDKit env) ; capture the formal charge ----
echo "===== GOAT $name / $solvent | $(date) ====="
seedlog=$(conda run -n "$SEED_ENV" python scripts/make_goat_seed.py --name "$name" --out "$work/seed.xyz")
echo "$seedlog"
[ -f "$work/seed.xyz" ] || { echo "ERROR: seed generation failed for $name" >&2; exit 1; }
charge=$(echo "$seedlog" | grep -oE 'SEED_CHARGE=-?[0-9]+' | cut -d= -f2); charge="${charge:-0}"

# ---- 2) GOAT (ORCA env: censo/openmpi python + ORCA binary; ORCA's GFN2-xTB is bundled) ----
conda activate orca
export PATH="$(dirname "$ORCA"):$PATH"
export OMPI_MCA_rmaps_base_oversubscribe=1     # belt-and-suspenders vs OpenMPI slot accounting under SLURM
[ -x "$ORCA" ] || { echo "ERROR: ORCA not executable at $ORCA" >&2; exit 1; }

# ORCA GOAT input (VERIFY syntax on this build; see header)
cat > "$work/goat.inp" <<INP
! GOAT GFN2-xTB ALPB($solvent)
%pal nprocs $JOBS end
%maxcore 3000
* xyzfile $charge 1 seed.xyz
INP

echo "  charge=$charge  ORCA=$ORCA"
( cd "$work" && "$ORCA" goat.inp > goat.out 2>&1 )
grep -qiE 'ORCA TERMINATED NORMALLY' "$work/goat.out" \
    || { echo "ERROR: ORCA GOAT did not finish (tail of goat.out):" >&2; tail -25 "$work/goat.out" >&2; exit 1; }

# ---- 3) convert ORCA final ensemble -> our schema ----
python scripts/goat_to_ensemble.py --goat-dir "$work" --out-dir "$work" --solvent "$solvent" --charge "$charge"
echo "===== Done $name / $solvent -> $work/ensemble.xyz | $(date) ====="
