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
#   Check everything:     bash   scripts/tests/run_all.sh
#
# PORTABLE BY DESIGN — works on an HPC (conda + SLURM) or on Google Colab (no conda; deps
# come from pip). Nothing here is fatal: if conda or xtb is missing we warn and continue,
# so the stages that don't need them still run. tests/check_env.sh is what tells you
# whether the stage you actually want is ready.
#
# Overridable via environment: CHAMELEON_CONDA_ENV, XTB_DIST.
# NOTE: this file intentionally does NOT set OMP_NUM_THREADS — threading is job-specific
# (scoring wants 1 thread/xtb x N parallel confs; CREST manages its own). Set it per-job.
# ---------------------------------------------------------------------------------------

# --- conda: activate the pipeline env, if conda exists ---------------------------------
# Looks for conda via $CONDA_EXE first, then the usual install locations (HPC miniconda,
# Colab/docker /opt/conda). On Colab there is typically no conda and deps are pip-installed
# into the system python — that is fine, we just skip activation.
CONDA_ENV="${CHAMELEON_CONDA_ENV:-chameleon_crest212}"
_conda_sh=""
if [ -n "${CONDA_EXE:-}" ] && [ -x "${CONDA_EXE:-}" ]; then
    _conda_sh="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
fi
if [ ! -f "${_conda_sh:-/nonexistent}" ]; then
    for _base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge" "/opt/conda" "/usr/local"; do
        if [ -f "$_base/etc/profile.d/conda.sh" ]; then _conda_sh="$_base/etc/profile.d/conda.sh"; break; fi
    done
fi
if [ -f "${_conda_sh:-/nonexistent}" ]; then
    # shellcheck disable=SC1090
    source "$_conda_sh"
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
        conda activate "$CONDA_ENV"
    else
        echo "[env.sh] NOTE: conda env '$CONDA_ENV' not found — using the current python." >&2
    fi
else
    echo "[env.sh] NOTE: no conda found — using the current python (normal on Colab/pip setups)." >&2
fi

# --- xtb: prefer the CPCM-X-enabled release binary over any conda/system one -----------
XTB_DIST="${XTB_DIST:-$HOME/xtb-dist}"
if [ -x "$XTB_DIST/bin/xtb" ]; then
    export PATH="$XTB_DIST/bin:$PATH"
    export XTBPATH="$XTB_DIST/share/xtb"
else
    echo "[env.sh] WARNING: $XTB_DIST/bin/xtb not found — run 'bash scripts/setup_xtb.sh' first." >&2
    echo "[env.sh]          Falling back to '$(command -v xtb || echo none)' (may lack CPCM-X)." >&2
fi

# --- report which binaries are live ----------------------------------------------------
_xtb_bin="$(command -v xtb || true)"
if [ -n "$_xtb_bin" ]; then
    _xtb_ver="$("$_xtb_bin" --version 2>&1 | grep -i version | head -1 | sed 's/^ *//')"
else
    _xtb_bin="NONE"; _xtb_ver="not installed"
fi
echo "[env.sh] python=$(command -v python || echo NONE)"
echo "[env.sh] xtb=$_xtb_bin  ($_xtb_ver)"
