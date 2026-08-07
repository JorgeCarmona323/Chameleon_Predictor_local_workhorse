#!/bin/bash
#SBATCH --job-name=fe_hexpep_cpcmx3
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --time=06:00:00
#
# CPCM-X ΔG_transfer for HexPep across ALL THREE solvents (water, chloroform, hexane), to complete
# the CPCM-X-vs-SMD methods comparison. We already have CPCM-X water/hexane (-7.74); this adds the
# CHLOROFORM leg — the literature-standard membrane mimic (Rezai/EnsembleCycPerm) and the leg where
# SMD gave a sane -6.75. Scores the same CREST ensembles SMD used (no re-search). --ewin 8 pre-trim.

set -euo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/free_energy results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# latest HexPep ensemble per solvent (water+chcl3 and hexane may live in different run dirs)
leg() { ls -dt results/runs/run_*HexPep/"$1" 2>/dev/null | head -1; }
W="$(leg water)";  W="${W:+$W/ensemble.xyz}"
C="$(leg chcl3)";  C="${C:+$C/ensemble.xyz}"
H="$(leg hexane)"; H="${H:+$H/ensemble.xyz}"
[ -n "${W:-}" ] && [ -f "$W" ] || { echo "ERROR: no HexPep water ensemble"  >&2; exit 1; }
[ -n "${C:-}" ] && [ -f "$C" ] || { echo "ERROR: no HexPep chcl3 ensemble"  >&2; exit 1; }
echo "water  = $W"; echo "chcl3  = $C"; echo "hexane = ${H:-<none>}"

LEGS=(--leg "water=$W" --leg "chloroform=$C")
[ -n "${H:-}" ] && [ -f "$H" ] && LEGS+=(--leg "hexane=$H")

python scripts/free_energy_calculator.py \
    --method cpcmx --ewin 8 --ref water --charge 0 --jobs "$JOBS" \
    "${LEGS[@]}" \
    --out results/free_energy/hexpep_cpcmx_3solv.csv

echo
echo "Done. CPCM-X ΔG_transfer -> results/free_energy/hexpep_cpcmx_3solv.summary.csv"
echo "Compare water->chloroform to SMD's -6.75; water->hexane should reconfirm ~-7.74."
