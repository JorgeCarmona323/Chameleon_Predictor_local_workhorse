#!/bin/bash
# run_all.sh — run every check/test in this directory and print ONE combined summary.
# ---------------------------------------------------------------------------------------
# Auto-discovers any `check_*.sh` (environment/preflight) and `test_*.py` (functional) in
# scripts/tests/ — so to add a test, just drop a correctly-named file here; no wiring needed.
# Exits 0 iff every discovered check/test passes. This is the one command a new user runs
# after setup_xtb.sh to confirm they have everything they need to roll.
#
#   Usage:  bash scripts/tests/run_all.sh
# ---------------------------------------------------------------------------------------
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

rc=0
ran=0
shopt -s nullglob
for t in check_*.sh test_*.py; do
    ran=$((ran+1))
    echo "########################################################################"
    echo "# $t"
    echo "########################################################################"
    case "$t" in
        *.sh) bash   "$t" || rc=1 ;;
        *.py) python "$t" || rc=1 ;;
    esac
    echo
done

if [ "$ran" -eq 0 ]; then
    echo "no check_*.sh / test_*.py found in $SCRIPT_DIR"; exit 0
fi

echo "========================================================================"
if [ "$rc" -eq 0 ]; then
    echo "ALL CHECKS PASSED ($ran run) — you're good to roll."
else
    echo "SOME CHECKS FAILED — see the output above."
fi
exit $rc
