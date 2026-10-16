#!/usr/bin/env bash
# Download the POG570 files used in this project into Data/Raw.
# Research use only, do not share or commit these files (see docs/data.md).
set -euo pipefail

URL="https://www.bcgsc.ca/downloads/POG570"
OUT="Data/Raw"
FILES=(
    POG570_README.txt
    POG570_TPM_expression.txt.gz
    POG570_small_mutations.txt.gz
    Table_S1_Demographics.xlsx
    Table_S2_Treatment.xlsx
)

mkdir -p "$OUT"
for f in "${FILES[@]}"; do
    if [[ -s "$OUT/$f" ]]; then
        echo "already have $f, skipping"
    else
        echo "downloading $f"
        curl -fSL --retry 3 -o "$OUT/$f" "$URL/$f"
    fi
done

# check the gzip files are complete (a broken download fails here)
for f in "$OUT"/*.gz; do
    gzip -t "$f" && echo "ok: $f"
done

# record checksums so a later download can be compared with this one
sha256sum "$OUT"/POG570_* "$OUT"/Table_S* > "$OUT/checksums.sha256"
echo "checksums written to $OUT/checksums.sha256"
