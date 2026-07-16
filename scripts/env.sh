#!/bin/bash
# env.sh — single source of truth for the Chameleon pipeline runtime.
# ---------------------------------------------------------------------------------------
# SOURCE this at the top of every SLURM wrapper (generation AND scoring) so conformer
# generation and free-energy scoring run on the SAME xtb binary. That binary is the
# official grimme-lab xtb release ("xtb-dist"), which is built WITH the CPCM-X solvation
# library. The conda xtb bundled with CREST 2.12 was compiled WITHOUT CPCM-X
# ("CPCM-X library was not included in this version of xTB"), so we prepend xtb-dist to
# PATH to prefer it. It is the SAME xtb version/commit the conda env pairs with CREST, so
# this is a drop-in for generation (no re-validation) and it unlocks --cpcmx for scoring.
#
#   Usage in a wrapper:   source scripts/env.sh        (after `cd` to the repo root)
#   One-time setup:       bash   scripts/setup_xtb.sh  (installs xtb-dist)
#
# Overridable via environment: CHAMELEON_CONDA_ENV, XTB_DIST.
# NOTE: this file intentionally does NOT set OMP_NUM_THREADS — threading is job-specific
# (scoring wants 1 thread/xtb x N parallel confs; CREST manages its own). Set it per-job.
# ---------------------------------------------------------------------------------------

# --- conda: CREST 2.12 + python deps ---------------------------------------------------
CONDA_ENV="${CHAMELEON_CONDA_ENV:-chameleon_crest212}"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# --- xtb: prefer the CPCM-X-enabled release binary over the conda one ------------------
XTB_DIST="${XTB_DIST:-$HOME/xtb-dist}"
if [ -x "$XTB_DIST/bin/xtb" ]; then
    export PATH="$XTB_DIST/bin:$PATH"
    export XTBPATH="$XTB_DIST/share/xtb"
else
    echo "[env.sh] WARNING: $XTB_DIST/bin/xtb not found — run 'bash scripts/setup_xtb.sh' first." >&2
    echo "[env.sh]          Falling back to '$(command -v xtb || echo none)' (may lack CPCM-X)." >&2
fi

# --- report which binaries are live ----------------------------------------------------
echo "[env.sh] conda=$CONDA_ENV  xtb=$(command -v xtb || echo NONE)  ($(xtb --version 2>&1 | grep -i version | head -1 | sed 's/^ *//'))"
