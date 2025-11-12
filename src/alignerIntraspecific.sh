#!/bin/bash
set -euo pipefail

mkdir -p data/speciesWork/Capsicum/alignedProteins

export MAFFT_BIN="/programs/mafft/bin/mafft"

find data/speciesWork/Capsicum/proteinsByGene/ -name "*.fasta" | parallel -j 20 '
    gene=$(basename {} .fasta)
    outfile="data/speciesWork/Capsicum/alignedProteins/${gene}_aligned.fasta"
    echo "Aligning {} -> $outfile"
    $MAFFT_BIN --adjustdirection --thread 4 --auto --treeout {} > "$outfile"
'