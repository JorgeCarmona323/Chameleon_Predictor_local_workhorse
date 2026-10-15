#!/bin/bash
# Stage ONLY the files needed downstream (final CREST ensembles + CPCM-X energy CSVs) for the
# 32 ortho 6-mer hits into ONE folder + a single tarball, for easy one-drag download via
# Termius/SFTP.  Deliberately skips the big CREST/xTB scratch.  Each molecule's LATEST run dir
# is found by matching "name" in metadata.json (latest-timestamp convention), so the charged
# 1-6-4-7 (newest) wins over the old neutral/xylene runs.
#
#   Run on the HPC from the repo root:  bash scripts/stage_hits_transfer.sh
#   Output:  ~/Chameleon_Predictor/hits_transfer.tgz   (download this one file via Termius)
#
# Local unpack (Git Bash / WSL, from the repo root):
#     tar xzf hits_transfer.tgz
#   -> hits_transfer/conformers/<name>/{water,hexane}/ensemble.xyz(+.sdf)+metadata.json
#      hits_transfer/free_energy/fe_hit_<name>.csv
# Folders use the BARE code name (no "DOPC " prefix) on purpose, so unpacking never overwrites
# your existing results/conformers/DOPC 1-6-4-7 (neutral) or DOPC 6-5-8-12.  The 1-6-4-7 here
# is the CHARGED (-1) run.

set -uo pipefail
cd "$HOME/Chameleon_Predictor"
CSV="data/hits_ortho_batch.csv"
[ -f "$CSV" ] || { echo "ERROR: $CSV not found (run from repo root on the HPC)" >&2; exit 1; }

OUT="hits_transfer"
rm -rf "$OUT" "${OUT}.tgz"
mkdir -p "$OUT/conformers" "$OUT/free_energy"

names=$(awk -F',' 'NR>1 { for(i=1;i<=NF;i++) gsub(/\r/,"",$i); if($4==1) print $2 }' "$CSV")
missing=""
n_ok=0

for name in $names; do
    base=""
    for d in $(ls -dt results/runs/*/ 2>/dev/null); do
        if grep -q "\"name\": \"$name\"" "${d}water/metadata.json" 2>/dev/null \
           && [ -f "${d}water/ensemble.xyz" ] && [ -f "${d}hexane/ensemble.xyz" ]; then
            base="${d%/}"; break
        fi
    done
    if [ -z "$base" ]; then
        echo "  MISSING complete water+hexane run: $name"
        missing="$missing $name"; continue
    fi
    dest="$OUT/conformers/$name"
    mkdir -p "$dest/water" "$dest/hexane"
    cp "$base"/water/ensemble.*  "$base"/water/metadata.json  "$dest/water/"  2>/dev/null
    cp "$base"/hexane/ensemble.* "$base"/hexane/metadata.json "$dest/hexane/" 2>/dev/null
    cp results/free_energy/fe_hit_"$name".csv         "$OUT/free_energy/" 2>/dev/null
    cp results/free_energy/fe_hit_"$name".summary.csv "$OUT/free_energy/" 2>/dev/null
    # record which run dir this came from (provenance)
    echo "$name  <-  $base" >> "$OUT/_source_run_dirs.txt"
    echo "  staged $name  <-  $base"
    n_ok=$((n_ok + 1))
done

tar czf "${OUT}.tgz" "$OUT"
echo
echo "===================================================================="
echo "Staged $n_ok / $(echo "$names" | wc -w) hits."
echo "Download this ONE file via Termius:"
echo "    $(pwd)/${OUT}.tgz   ($(du -sh "${OUT}.tgz" | cut -f1))"
[ -n "$missing" ] && echo "MISSING (not staged -- incomplete runs):$missing"
echo "===================================================================="
