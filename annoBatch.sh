#!/bin/bash
set -euo pipefail

mkdir -p /workdir/hdd29/tmp
export TMP_BASE_DIR="/workdir/hdd29/tmp"

input_dir="/workdir/hdd29/theRefseqening/data/genomes"
final_dir="/workdir/hdd29/theRefseqening/data/results_batched"
singularity_image="/local/workdir/hdd29/singularity_images/cpgavas2_sandbox"
batch_size=100

mkdir -p "$final_dir" "$TMP_BASE_DIR"

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate singularity

if [ $# -eq 0 ]; then
  echo "Usage: $0 <file_list>"
  exit 1
fi

file_list="$1"
split_dir="$TMP_BASE_DIR/batches"
rm -rf "$split_dir"
mkdir -p "$split_dir"

split -l $batch_size "$file_list" "$split_dir/batch_"

process_batch() {
  batch_file="$1"
  batch_id=$(basename "$batch_file")
  job_tmp="$TMP_BASE_DIR/${batch_id}_tmp"
  mkdir -p "$job_tmp"

  run_log="$final_dir/${batch_id}.log"
  BATCH_FILE_ABS=$(realpath "$batch_file")

  {
    echo "[START] $(date)"
    df -h "$job_tmp"
    SECONDS=0

    SINGULARITYENV_BATCH_FILE="/mnt/$(basename "$BATCH_FILE_ABS")" /usr/bin/singularity exec \
    --bind "$job_tmp:/mnt" \
    --bind "$input_dir:/input" \
    --bind "$final_dir:/output" \
    --bind "$BATCH_FILE_ABS:/mnt/$(basename "$BATCH_FILE_ABS")" \
    "$singularity_image" bash -c '
        while read basename; do
        fasta="/input/${basename}.fa"
        pid=$(echo "$basename" | tr -d "_.")
        outdir="/output/$pid"
        mkdir -p "$outdir"
        cp "$fasta" /mnt/input.fa
        run-cpgavas2 -in /mnt/input.fa -db 1 -pid "$pid" -out "$outdir"
        cp -r /mnt/* "$outdir/"
        rm -rf /mnt/*
        done < "$BATCH_FILE"
    '

    echo "[END] $(date) Duration: ${SECONDS}s"
    df -h "$job_tmp"
    rm -rf "$job_tmp"
  } >> "$run_log" 2>&1
}

export -f process_batch
export singularity_image input_dir final_dir TMP_BASE_DIR

find "$split_dir" -type f | parallel --progress -j 5 process_batch {}
