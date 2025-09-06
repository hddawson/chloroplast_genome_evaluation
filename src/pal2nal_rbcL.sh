#!/usr/bin/env bash
set -euo pipefail

PAL2NAL=/programs/pal2nal/pal2nal.pl
AA_DIR="data/tmp/alignedGenes"
CDS_DIR="data/tmp/genesToAlign"
OUT_DIR="data/tmp/pal2nal"
GENETIC_CODE=11

mkdir -p "$OUT_DIR"

for aa in "$AA_DIR"/*_AA_aligned.fasta; do
    gene=$(basename "$aa" _AA_aligned.fasta)
    cds="$CDS_DIR/${gene}_CDS.fasta"
    out="$OUT_DIR/${gene}_pal2nal.fa"

    if [[ ! -f "$cds" ]]; then
        echo "Skipping $gene: missing CDS file" >&2
        continue
    fi

    echo "Processing $gene..."
    perl "$PAL2NAL" "$aa" "$cds" -codontable "$GENETIC_CODE" -output fasta > "$out"
done

echo "Done. Results in $OUT_DIR/"
