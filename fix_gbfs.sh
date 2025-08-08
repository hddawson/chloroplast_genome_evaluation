#!/bin/bash

# Define working directory
gbf_dir="/workdir/hdd29/chloroplast_genome_evaluation/data/gbfs"

# Loop through each .gbf file
for file in "$gbf_dir"/*.gbf; do
    # Use sed to replace the first line with the new LOCUS line
    sed -i '1s/.*/LOCUS       PLACEHOLDER 00000001 bp   DNA  circular PLN 05-AUG-2025/' "$file"
done
