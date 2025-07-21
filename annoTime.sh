#!/bin/bash
set -e

# Set directories
input_dir="/workdir/hdd29/theRefseqening/data/genomes"
final_dir="/workdir/hdd29/theRefseqening/data/results"
temp_base_dir="/workdir/hdd29/tmp"
singularity_image="cpgavas2_0.03.sif"
bind_dir="/workdir/hdd29/theRefseqening/data/genomes"

mkdir -p "$final_dir"
mkdir -p "$temp_base_dir"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate singularity

export TMPDIR="$temp_base_dir"

# Function to process each file
process_file() {
  fasta_file="$1"
  base_name=$(basename "$fasta_file" .fa)
  pid_name=$(echo "$base_name" | tr -d '_.')
  final_output_dir="$final_dir/${pid_name}"
  
  if [[ -d "$final_output_dir" && -n "$(ls -A "$final_output_dir" 2>/dev/null)" ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [SKIP] $base_name: output already exists."
    return 0
  fi

  echo "$(date '+%Y-%m-%d %H:%M:%S') [START] $base_name"
  mkdir -p "$final_output_dir"

  run_log="$final_output_dir/run.log"
  {
    echo "[START] $(date '+%Y-%m-%d %H:%M:%S')"
    SECONDS=0

    singularity exec \
      --containall \
      --writable-tmpfs \
      --bind "$bind_dir:/mnt" \
      --bind "$final_output_dir:/output" \
      "$singularity_image" run-cpgavas2 \
      -in "/mnt/$(basename "$fasta_file")" \
      -db 1 \
      -pid "$pid_name" \
      -out "/output"

    rc=$?
    duration=$SECONDS
    echo "[END]   $(date '+%Y-%m-%d %H:%M:%S') Duration: ${duration}s"

    if [ $rc -eq 0 ]; then
      echo "[DONE]  $base_name completed successfully."
    else
      echo "[FAIL]  $base_name failed."
    fi
    exit $rc
  } >> "$run_log" 2>&1

  if [ $? -eq 0 ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [DONE] $base_name"
  else
    echo "$(date '+%Y-%m-%d %H:%M:%S') [FAIL] $base_name — check $run_log"
    return 1
  fi
}

export -f process_file
export singularity_image bind_dir final_dir temp_base_dir TMPDIR

# Process in parallel
find "$input_dir" -name "*.fa" | parallel --progress -j 20 process_file {}
