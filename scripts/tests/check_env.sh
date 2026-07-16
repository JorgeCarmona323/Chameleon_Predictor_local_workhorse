#!/bin/bash
# check_env.sh — preflight "doctor" for the Chameleon pipeline.
# ---------------------------------------------------------------------------------------
# Verifies a fresh checkout has everything needed to GENERATE conformers and SCORE energies.
# Sources scripts/env.sh (conda + CPCM-X-enabled xtb), then checks each dependency and prints
# PASS / FAIL / WARN. Exits 0 iff every REQUIRED check passes (optional ones only WARN).
# Run it anytime; run it first. New user flow: setup_xtb.sh -> this -> you're ready.
#
#   Usage:  bash scripts/tests/check_env.sh      (works from any directory)
# ---------------------------------------------------------------------------------------
set -uo pipefail            # NOT -e: run EVERY check and report; don't stop at first failure

# resolve repo root from this file's location, so it runs from anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_DIR"

echo "=== Chameleon pipeline environment check ==="
echo "repo: $REPO_DIR"
source scripts/env.sh || { echo "FAIL  could not source scripts/env.sh"; exit 1; }
echo

pass=0; fail=0; warn=0

# check "description" cmd...    -> PASS if cmd exits 0, else FAIL (required)
check() { local d="$1"; shift
    if "$@" >/dev/null 2>&1; then printf '  PASS  %s\n' "$d"; pass=$((pass+1))
    else printf '  FAIL  %s\n' "$d"; fail=$((fail+1)); fi; }

# check_opt "description" cmd...  -> PASS or WARN (optional; never fails the run)
check_opt() { local d="$1"; shift
    if "$@" >/dev/null 2>&1; then printf '  PASS  %s\n' "$d"; pass=$((pass+1))
    else printf '  WARN  %s (optional)\n' "$d"; warn=$((warn+1)); fi; }

# does this xtb actually have CPCM-X compiled in? (the whole saga — verify, don't assume)
xtb_has_cpcmx() {
    local tmp; tmp=$(mktemp -d) || return 1
    printf '3\n\nO 0 0 0\nH 0.757 0.586 0\nH -0.757 0.586 0\n' > "$tmp/h2o.xyz"
    local rc=1
    (cd "$tmp" && xtb h2o.xyz --gfn 2 --cpcmx water --chrg 0 >out 2>&1) && rc=0
    rm -rf "$tmp"; return $rc
}

echo "--- required ---"
check "python on PATH"                 command -v python
check "rdkit importable (generation)"  python -c "import rdkit"
check "xtb on PATH"                     command -v xtb
check "XTBPATH is set"                  test -n "${XTBPATH:-}"
check "xtb has CPCM-X (--cpcmx water)"  xtb_has_cpcmx

echo "--- optional (only needed for the ORCA SMD/COSMO arm) ---"
check_opt "orca on PATH"                command -v orca

echo
echo "=== summary: $pass passed, $fail failed, $warn warning(s) ==="
if [ "$fail" -gt 0 ]; then
    echo "Not ready. Fix the FAIL lines above:"
    echo "  - xtb / CPCM-X issues  -> run: bash scripts/setup_xtb.sh"
    echo "  - rdkit / python       -> activate the conda env (see scripts/env.sh: CHAMELEON_CONDA_ENV)"
    exit 1
fi
echo "Ready to generate conformers and score energies."
exit 0
