#!/usr/bin/env bash

# Script to get occurrences of species from GBIF and BIEN using gnu parallel to call pull_species_occurrences.r

#read in the query file 
query_file=$1
output_dir=$2
#echo both
echo "Query file: $query_file"
echo "Output directory: $output_dir"

mkdir -p "$output_dir"


#echo each occurence name to stdout
while read -r line; do
    echo "$line"
done < "$query_file" #| xargs -I {} pixi run Rscript src/occurrenceQuerier.r "{}" "$output_dir"

cat "$query_file" | parallel -j 7 pixi run Rscript src/occurrenceQuerier.r "{}" "$output_dir"
