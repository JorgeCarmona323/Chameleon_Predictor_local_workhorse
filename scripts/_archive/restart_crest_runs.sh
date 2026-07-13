#!/bin/bash
#SBATCH --job-name=crest_restart
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Restart CREST for CsA — all 14 MTDs completed but post-MTD processing was
# killed by the 24h script timeout. crest_rotamers_0.xyz (853 MB) is present
# for both aq and mem; --restart resumes multilevel optimization from structure
# 4271/69986. DP955 and DP944 MTDs were incomplete — submit those fresh.

set -e

REPO="$HOME/Chameleon_Predictor"
THREADS=20

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate chameleon_crest212
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

RUNS=(
    "run_20260429_173210_1_CsA:CsA_aq_start.xyz:water:water"
    "run_20260429_173210_1_CsA:CsA_mem_start.xyz:chcl3:mem"
)

for entry in "${RUNS[@]}"; do
    IFS=':' read -r run_id start_xyz solvent label <<< "$entry"

    crest_dir="$REPO/results/crest_runs/$run_id/$label/crest"
    ensemble_out="$REPO/results/crest_runs/$run_id/$label/ensemble.xyz"

    echo "======================================================="
    echo "  Run:     $run_id"
    echo "  Solvent: $solvent ($label)"
    echo "  Dir:     $crest_dir"
    echo "  Started: $(date)"
    echo "======================================================="

    if [ ! -d "$crest_dir" ]; then
        echo "  ERROR: directory not found — skipping"
        continue
    fi

    cd "$crest_dir"

    crest "$start_xyz" \
        -T $THREADS \
        --gfn2 \
        --chrg 0 \
        --alpb "$solvent" \
        --keepdir \
        --restart \
        > crest_restart.out 2>&1

    EXIT=$?
    echo "  CREST exit=$EXIT at $(date)"

    if [ -f "crest_conformers.xyz" ] && [ -s "crest_conformers.xyz" ]; then
        cp crest_conformers.xyz "$ensemble_out"
        n_confs=$(grep -c "^[0-9]" crest_conformers.xyz 2>/dev/null || echo "?")
        echo "  Ensemble saved → $ensemble_out"
        echo "  Conformers: $n_confs"
    else
        echo "  WARNING: no crest_conformers.xyz produced"
    fi

    echo ""
done

echo "All restarts complete: $(date)"
