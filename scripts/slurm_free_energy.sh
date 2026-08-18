#!/bin/bash
#SBATCH --job-name=fe_cpcmx
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --partition=all
#
# ΔG_transfer scoring — GFN2-xTB + CPCM-X single-points on the native per-phase
# CREST ensembles (no re-search, no re-optimization).
#
# Auto-discovers every compound under CONFORMERS_ROOT that has a hexane/ensemble.xyz
# and pairs it with the sibling water/ and chcl3/ ensembles. One free_energy_calculator
# run per compound scores all its phases together and reports:
#     ΔG_transfer(water -> chloroform)  and  ΔG_transfer(water -> hexane)
# CPCM-X solvent keywords used: water, chloroform, hexane (all in the CPCM-X DB).
#
# xTB is Linux-only, so this runs on the HPC. The conformer tree must be present here:
# the curated results/conformers/ was organized locally, so rsync/scp it to the HPC
# (under $REPO_DIR/results/conformers) before submitting.
#
# Cost: single-points are seconds each and run --jobs-parallel, and --ewin trims each
# ensemble to its low-energy window first — so all 14 finish in ~tens of minutes, not
# the hours CREST took. Submit with:  sbatch scripts/slurm_free_energy.sh

set -euo pipefail

REPO_DIR="$HOME/Chameleon_Predictor"
CONFORMERS_ROOT="$REPO_DIR/results/conformers"
OUTDIR="$REPO_DIR/results/free_energy"
EWIN=8                                   # kcal/mol GFN2 pre-filter before CPCM-X
JOBS="${SLURM_CPUS_PER_TASK:-20}"

cd "$REPO_DIR"
mkdir -p "$OUTDIR" results/slurm_logs

source scripts/env.sh                    # conda env + CPCM-X-enabled xtb (single source of truth)
export OMP_NUM_THREADS=1                 # 1 thread/xtb → run $JOBS single-points in parallel

echo "Discovering compounds under: $CONFORMERS_ROOT"
mapfile -t HEXDIRS < <(find "$CONFORMERS_ROOT" -type d -name hexane | sort)
echo "Found ${#HEXDIRS[@]} compound(s) with a hexane ensemble:"
printf '  %s\n' "${HEXDIRS[@]}"
echo

for hexdir in "${HEXDIRS[@]}"; do
    cdir=$(dirname "$hexdir")            # the compound directory

    # skip compounds scored elsewhere (CsA has its own GFN-FF wrapper; HexPep/validation
    # compounds use their own aq/-named legs) so this batch stays to the project hits.
    case "$hexdir" in
        *[Cc]yclosporin*|*_CsA*|*HexPep*|*Roxithromycin*|*Begnini*)
            echo "skip (handled separately): $hexdir"; continue ;;
    esac

    # label from the hexane manifest filename (e.g. DOPCdz_R_manifest.json -> DOPCdz_R)
    manifest=$(find "$hexdir" -maxdepth 1 -name '*_manifest.json' | head -1)
    if [ -n "$manifest" ]; then
        label=$(basename "$manifest" _manifest.json)
    else
        label=$(echo "${cdir#"$CONFORMERS_ROOT"/}" | tr ' /' '__')
    fi

    # match the charge used at generation (read from the hexane metadata; default 0)
    charge=$(grep -oE '"charge"[^0-9-]*(-?[0-9]+)' "$hexdir/metadata.json" 2>/dev/null \
             | grep -oE '\-?[0-9]+$' || true)
    charge=${charge:-0}

    # build legs from whichever phases exist (water + chloroform + hexane)
    legs=()
    [ -f "$cdir/water/ensemble.xyz" ] && legs+=(--leg "water=$cdir/water/ensemble.xyz")
    [ -f "$cdir/chcl3/ensemble.xyz" ] && legs+=(--leg "chloroform=$cdir/chcl3/ensemble.xyz")
    legs+=(--leg "hexane=$hexdir/ensemble.xyz")

    echo "=== $label  |  $(( ${#legs[@]} / 2 )) legs  |  charge $charge  |  $(date) ==="
    python scripts/free_energy_calculator.py \
        --method cpcmx --ewin "$EWIN" --ref water --charge "$charge" --jobs "$JOBS" \
        "${legs[@]}" \
        --out "$OUTDIR/fe_${label}.csv"
done

echo
echo "All done. Per-compound results (+ .summary.csv with ΔG_transfer) in: $OUTDIR"
