#!/bin/bash
#SBATCH --job-name=run_validation
#SBATCH --output=results/slurm_logs/%x_%j.out
#SBATCH --error=results/slurm_logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=24G
#SBATCH --partition=all
# (NO #SBATCH --time -- project convention: never cap walltime on this cluster.)
#
# ONE driver for the whole validation set. Replaces the five bespoke fe_* wrappers
# (fe_begnini / csa_gfnff / hexpep_* / fe_roxi / fe_whc3), each of which had its own
# output name, override var, and copy of the same bugs. For every molecule it:
#   1. locates the newest CREST ensemble dir (override with <MOL>_DIR=...),
#   2. NORMALIZES the (inconsistently-named) solvent legs into
#        results/validation/<Mol>/{water,chcl3,hexane}   (symlinks -- no copy, no regen),
#   3. runs the SINGLE core pipeline  `pipeline.py --run-dir`  on it, which scores BOTH
#      legs (water->chloroform + water->hexane) AND computes the 3D descriptors.
# Uniform outputs per molecule, same names for all:
#   results/validation/<Mol>/free_energy_cpcmx.summary.csv
#   results/validation/<Mol>/descriptors_chcl3.csv  +  descriptors_hexane.csv
#   results/validation/<Mol>/pipeline_manifest.json
#
# All six validation compounds are neutral -> --charge 0. Scoring is uniform GFN2/CPCM-X
# regardless of how the ensembles were generated (GFN2 for the 6-mers, GFN-FF for CsA/roxi).

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
mkdir -p results/validation results/slurm_logs
source scripts/env.sh
export OMP_NUM_THREADS=1                  # 1 thread/xtb -> $JOBS single-points in parallel
JOBS="${SLURM_CPUS_PER_TASK:-20}"

# --- newest candidate dir (from an array of literal paths/globs) that has a water leg -----
# Safe: never dies on a non-matching glob or a missing sibling (the bug that killed roxi).
resolve_base() {
    local mol="$1"; shift
    local ovr_name="${mol}_DIR" ovr; ovr="${!ovr_name:-}"
    if [ -n "$ovr" ]; then echo "$ovr"; return; fi
    local d best="" bt=0 t
    for d in "$@"; do
        [ -d "$d" ] || continue
        ls "$d"/water*/ensemble.xyz >/dev/null 2>&1 || continue   # must have a water leg
        t=$(stat -c %Y "$d" 2>/dev/null || echo 0)
        if [ -z "$best" ] || [ "$t" -gt "$bt" ]; then best="$d"; bt="$t"; fi
    done
    echo "$best"
}

# --- symlink base/<solvent-ish>/ -> results/validation/<Mol>/{water,chcl3,hexane} --------
stage_legs() {
    local base="$1" dest="$2" w a h
    rm -rf "$dest"; mkdir -p "$dest"
    w=$(ls -d "$base"/water*  2>/dev/null | head -1)
    a=$(ls -d "$base"/chcl3* "$base"/chloroform* "$base"/mem* 2>/dev/null | head -1)
    h=$(ls -d "$base"/hexane* 2>/dev/null | head -1)
    [ -n "$w" ] && ln -sfn "$(readlink -f "$w")" "$dest/water"
    [ -n "$a" ] && ln -sfn "$(readlink -f "$a")" "$dest/chcl3"
    [ -n "$h" ] && ln -sfn "$(readlink -f "$h")" "$dest/hexane"
    echo "  water  = ${w:-<none>}"
    echo "  chcl3  = ${a:-<none>}"
    echo "  hexane = ${h:-<none>}"
    [ -f "$dest/water/ensemble.xyz" ] && [ -f "$dest/chcl3/ensemble.xyz" ]   # need >=2 legs
}

run_one() {
    local mol="$1" base="$2"
    echo "=================================================================="
    if [ -z "$base" ]; then
        echo "SKIP $mol -- no ensemble dir found (still generating, or set ${mol}_DIR=<dir>)"
        return
    fi
    echo "$mol  <-  $base"
    local dest="results/validation/$mol"
    if ! stage_legs "$base" "$dest"; then
        echo "SKIP $mol -- fewer than 2 legs available at $base"
        return
    fi
    python pipeline.py --run-dir "$dest" --name "$mol" \
        --method cpcmx --ewin 8 --charge 0 --threads "$JOBS" \
        || echo "WARN $mol -- pipeline returned non-zero (see above)"
}

# --- the six validation compounds. Candidate dirs are LITERAL (arrays handle the space) ---
run_one Begnini_1     "$(resolve_base Begnini_1     $(ls -d results/runs/run_*_Begnini_1 2>/dev/null) results/conformers/Begnini_1)"
run_one Begnini_2     "$(resolve_base Begnini_2     $(ls -d results/runs/run_*_Begnini_2 2>/dev/null) results/conformers/Begnini_2)"
run_one Roxithromycin "$(resolve_base Roxithromycin $(ls -d results/runs/run_*_Roxithromycin 2>/dev/null) results/conformers/Roxithromycin)"
run_one HexPep        "$(resolve_base HexPep        $(ls -d results/runs/run_*_HexPep 2>/dev/null) results/conformers/HexPep)"
run_one WhC3          "$(resolve_base WhC3          $(ls -d results/runs/run_*_WhC3 2>/dev/null) results/conformers/WhC3)"
# CsA: its ensemble folder name contains a space -> quote it explicitly
run_one CsA           "$(resolve_base CsA $(ls -d results/runs/run_*_CsA 2>/dev/null) "results/conformers/Cyclosporin A GFN_FF")"

echo "=================================================================="
echo "DONE. Uniform per-molecule outputs under results/validation/<Mol>/:"
ls -1 results/validation/*/free_energy_cpcmx.summary.csv 2>/dev/null
