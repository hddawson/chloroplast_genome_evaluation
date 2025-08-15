#!/bin/bash
set -euo pipefail

mkdir -p data/tmp/genesAligned

export MAFFT_BIN="/programs/mafft/bin/mafft"

find data/tmp/genesToAlign -name "*.fasta" | parallel -j 10 '
    gene=$(basename {} .fasta)
    outfile="data/tmp/genesAligned/${gene}_aligned.fasta"
    echo "Aligning {} -> $outfile"
    $MAFFT_BIN --adjustdirection --reorder --thread 8 --auto {} > "$outfile"
'