#!/bin/bash
# setup_xtb.sh — one-time install of the official grimme-lab xtb release binary, which is
# built WITH the CPCM-X solvation library (the conda/CREST xtb is not). The whole pipeline
# (via env.sh) then uses this binary for BOTH conformer generation and free-energy scoring,
# so everything is one consistent, reproducible xtb. Idempotent: re-running is a no-op if a
# working xtb-dist is already present. Verifies CPCM-X support before declaring success.
#
#   Usage:  bash scripts/setup_xtb.sh
#   Override: XTB_VERSION (default 6.7.1), XTB_DIST (default $HOME/xtb-dist)
set -euo pipefail

XTB_VERSION="${XTB_VERSION:-6.7.1}"
XTB_DIST="${XTB_DIST:-$HOME/xtb-dist}"
TARBALL="xtb-${XTB_VERSION}-linux-x86_64.tar.xz"
URL="https://github.com/grimme-lab/xtb/releases/download/v${XTB_VERSION}/${TARBALL}"

if [ -x "$XTB_DIST/bin/xtb" ]; then
    echo "xtb already installed: $XTB_DIST/bin/xtb"
    "$XTB_DIST/bin/xtb" --version 2>&1 | grep -i version | head -1
else
    cd "$HOME"
    echo "Downloading $URL ..."
    wget -q --show-progress "$URL"
    echo "Extracting $TARBALL ..."
    tar -xf "$TARBALL"
    rm -f "$TARBALL"
    # release extracts to xtb-dist/ by default; if a versioned dir was used, locate it
    if [ ! -x "$XTB_DIST/bin/xtb" ]; then
        found=$(find "$HOME" -maxdepth 2 -type f -path "*/bin/xtb" 2>/dev/null | grep -v miniconda | head -1 || true)
        [ -n "$found" ] && XTB_DIST=$(dirname "$(dirname "$found")")
    fi
    [ -x "$XTB_DIST/bin/xtb" ] || { echo "ERROR: xtb binary not found after extract" >&2; exit 1; }
    echo "Installed: $XTB_DIST/bin/xtb"
    "$XTB_DIST/bin/xtb" --version 2>&1 | grep -i version | head -1
fi

# --- verify CPCM-X is actually compiled in (the whole point of using this binary) ------
echo "Verifying CPCM-X support ..."
export XTBPATH="$XTB_DIST/share/xtb"
tmp=$(mktemp -d)
cat > "$tmp/h2o.xyz" <<'XYZ'
3

O  0.00000  0.00000  0.00000
H  0.75700  0.58600  0.00000
H -0.75700  0.58600  0.00000
XYZ
if (cd "$tmp" && "$XTB_DIST/bin/xtb" h2o.xyz --gfn 2 --cpcmx water --chrg 0 >out 2>&1); then
    echo "  OK: '--cpcmx water' ran and terminated normally — CPCM-X is available."
    echo
    echo "Done. Every wrapper that sources scripts/env.sh will now use this xtb."
    echo "  (env.sh sets XTB_DIST=$XTB_DIST — override with the XTB_DIST env var if it moved.)"
else
    if grep -q "CPCM-X library was not included" "$tmp/out"; then
        echo "  ERROR: this xtb build LACKS CPCM-X — wrong tarball/build." >&2
    else
        echo "  WARNING: '--cpcmx' returned nonzero for another reason; last lines:" >&2
        tail -8 "$tmp/out" >&2
    fi
    rm -rf "$tmp"; exit 1
fi
rm -rf "$tmp"
