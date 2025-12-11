#!/bin/bash
set -euo pipefail

mkdir -p data/speciesWork/Salix/alignedProteins

export MAFFT_BIN="/programs/mafft/bin/mafft"

find data/speciesWork/Salix/proteinsByGene/ -name "*.fasta" | parallel -j 20 '
    gene=$(basename {} .fasta)
    outfile="data/speciesWork/Salix/alignedProteins/${gene}_aligned.fasta"
    echo "Aligning {} -> $outfile"
    $MAFFT_BIN --adjustdirection --thread 4 --auto --treeout {} > "$outfile"
'