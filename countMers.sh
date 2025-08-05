#!/bin/bash
set -euo pipefail

MEM=100M
THREADS=35
LIST="data/pathway_luis.txt"
OUTDIR="jellyfish_output"
mkdir -p "$OUTDIR"

export PATH=/programs/jellyfish-2.3.1/bin:$PATH

# Loop through k from 2 to 10
for K in {31..31}; do
  echo "=== Processing k=$K ==="
  KDIR="$OUTDIR/k$K"
  mkdir -p "$KDIR"

  while read -r LUI; do
    FASTA="../theRefseqening/theRefseqening/data/genomes/${LUI}.fa"
    echo "→ Processing $LUI from $FASTA"

    jf="$KDIR/${LUI}.jf"
    counts="$KDIR/${LUI}_counts.tsv"

    if [[ ! -f "$FASTA" ]]; then
      echo "⚠️  File not found: $FASTA" >&2
      continue
    fi

    jellyfish count -m "$K" -s "$MEM" -t "$THREADS" -C "$FASTA" -o "$jf"
    jellyfish dump -c "$jf" > "$counts"
  done < "$LIST"
done
