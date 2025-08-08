#!/bin/bash
set -euo pipefail

MEM=100M
THREADS=10
LIST="${1:-data/n_status_luis.txt}"
OUTDIR="jellyfish_output"
CACHEDIR="/dev/shm/genomes"

mkdir -p "$OUTDIR" "$CACHEDIR"
export PATH=/programs/jellyfish-2.3.1/bin:$PATH

for K in {12..12}; do
  echo "=== Processing k=$K ==="
  KDIR="$OUTDIR/k$K"
  mkdir -p "$KDIR"

  while read -r LUI; do
    SRC="../theRefseqening/theRefseqening/data/genomes/${LUI}.fa"
    CACHED="$CACHEDIR/${LUI}.fa"
    jf="$KDIR/${LUI}.jf"
    counts="$KDIR/${LUI}_counts.tsv"

    if [[ -f "$counts" ]]; then
      echo "✓ Skipping $LUI - output already exists: $counts"
      continue
    fi
    if [[ ! -f "$SRC" ]]; then
      echo "⚠️  File not found: $SRC" >&2
      continue
    fi
    if [[ ! -f "$CACHED" ]]; then
      cp "$SRC" "$CACHED"
    fi

    jellyfish count -m "$K" -s "$MEM" -t "$THREADS" -C "$CACHED" -o "$jf"
    jellyfish dump -c "$jf" > "$counts"
  done < "$LIST"
done
