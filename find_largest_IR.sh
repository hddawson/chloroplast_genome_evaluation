#!/bin/bash
set -euo pipefail

fasta_dir="../theRefseqening/theRefseqening/data/genomes/"
output_dir="data/vmatch_output"
mkdir -p "$output_dir"

export VMATCH_BIN="/programs/vmatch-2.3.0"
export OUTPUT_DIR="$output_dir"

process_fasta() {
  fasta="$1"
  name=$(basename "$fasta" .fa)
  index_dir="$OUTPUT_DIR/${name}_index"
  result_file="$OUTPUT_DIR/${name}_inverted.txt"
  mkdir -p "$index_dir"
  "$VMATCH_BIN"/mkvtree -db "$fasta" -dna -pl -allout -indexname "$index_dir/$name"
  "$VMATCH_BIN"/vmatch -v -l 2000 -p -identity 90 "$index_dir/$name" > "$result_file"
}
export -f process_fasta

find "$fasta_dir" -name "*.fa*" | parallel --eta -j 16 process_fasta
