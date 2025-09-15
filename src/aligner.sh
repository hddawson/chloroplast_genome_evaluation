#!/bin/bash
set -euo pipefail

mkdir -p data/tmp/alignedGenes

export MAFFT_BIN="/programs/mafft/bin/mafft"

find data/tmp/genesToAlign -name "*CDS.fasta" | parallel -j 20 '
    gene=$(basename {} .fasta)
    outfile="data/tmp/alignedGenes/${gene}_aligned.fasta"
    echo "Aligning {} -> $outfile"
    $MAFFT_BIN --adjustdirection --thread 4 --auto --treeout {} > "$outfile"
'