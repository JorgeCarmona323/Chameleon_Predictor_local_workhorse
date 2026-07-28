#!/bin/bash
#SBATCH --job-name=smd_hexpep
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --time=08:00:00
#
# Task 9 (SMD arm) — ORCA DFT + SMD ΔG_transfer for the HexPep 6-mer, to benchmark how far
# off xTB CPCM-X (-7.74 kcal/mol water->hexane) is. Scores the SAME CREST water/hexane/chcl3
# ensembles CPCM-X used; DFT single-points at r2SCAN-3c + SMD, --top 12 conformers/phase.
# Runs ORCA SERIAL per conformer (%pal nprocs 1) and parallelizes ACROSS conformers via
# --jobs, so no MPI is needed. Output lands next to the CPCM-X results for direct comparison.
#
# COSMO-RS arm is NOT here — ORCA only builds the sigma-surface; full COSMO-RS needs
# openCOSMO-RS or COSMOtherm (see memory project_cosmo_rs_backend_dependency).

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs

# --- ORCA (shared build) on PATH + lib path, INSIDE this job's shell only -----------------
ORCA_DIR="$HOME/orca_6.1.1/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg"
ORCA="$ORCA_DIR/orca"
export PATH="$ORCA_DIR:$PATH"
export LD_LIBRARY_PATH="$ORCA_DIR:${LD_LIBRARY_PATH:-}"
[ -x "$ORCA" ] || { echo "ERROR: ORCA not found/executable at $ORCA" >&2; exit 1; }

source scripts/env.sh                    # conda python (harmlessly also sets xtb-dist PATH)
JOBS="${SLURM_CPUS_PER_TASK:-12}"

# --- locate the HexPep ensembles (water/chcl3/hexane live in separate run dirs) -----------
# latest run_*_HexPep/<solvent>/ensemble.xyz per phase
leg_path() { ls -dt results/runs/run_*HexPep/"$1" 2>/dev/null | head -1; }
W="$(leg_path water)";  W="${W:+$W/ensemble.xyz}"
H="$(leg_path hexane)"; H="${H:+$H/ensemble.xyz}"
C="$(leg_path chcl3)";  C="${C:+$C/ensemble.xyz}"

[ -n "${W:-}" ] && [ -f "$W" ] || { echo "ERROR: no HexPep water ensemble found"  >&2; exit 1; }
[ -n "${H:-}" ] && [ -f "$H" ] || { echo "ERROR: no HexPep hexane ensemble found" >&2; exit 1; }
echo "water  = $W"
echo "hexane = $H"
LEGS=(--leg "water=$W" --leg "hexane=$H")
if [ -n "${C:-}" ] && [ -f "$C" ]; then echo "chcl3  = $C"; LEGS+=(--leg "chcl3=$C"); fi

python scripts/free_energy_orca_smd.py \
    --orca "$ORCA" --level "r2SCAN-3c" \
    --ewin 5 --top 12 --charge 0 --mult 1 --jobs "$JOBS" \
    --workdir results/free_energy/hexpep_smd_orca \
    "${LEGS[@]}" --ref water \
    --out results/free_energy/hexpep_smd.csv

echo
echo "Done. SMD ΔG_transfer -> results/free_energy/hexpep_smd.summary.csv"
echo "Compare against CPCM-X (HexPep water->hexane was -7.74 kcal/mol)."
echo "Sanity-check SMD engaged: grep -i 'CDS' results/free_energy/hexpep_smd_orca/conf*_water/conf.out | head"
