#!/bin/bash
set -euo pipefail

mkdir -p /workdir/hdd29/tmp
export TMP_BASE_DIR="/workdir/hdd29/tmp"

# Set directories
input_dir="/workdir/hdd29/theRefseqening/data/genomes"
final_dir="/workdir/hdd29/theRefseqening/data/results_good"
singularity_image="cpgavas2_0.03.sif"
bind_dir="$input_dir"

mkdir -p "$final_dir"
mkdir -p "$TMP_BASE_DIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate singularity

# Calculate max parallel jobs based on free space in /tmp (assume 10GB per job)
available_kb=$(df --output=avail -k /workdir/hdd29/tmp | tail -1)
max_jobs=$(( available_kb / (10 * 1024 * 1024) ))  # 10GB per job in KB to GB
max_jobs=$(( max_jobs > 20 ? 20 : (max_jobs < 1 ? 1 : max_jobs) ))

echo "Detected free space: $(( available_kb / 1024 / 1024 )) GB"
echo "Setting max parallel jobs to: $max_jobs"

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

  # Create a unique per-job temp directory
  job_tmp="$TMP_BASE_DIR/${pid_name}_tmp"
  mkdir -p "$job_tmp"
  export TMPDIR="$job_tmp"
  export TMP="$job_tmp"

  # Optional: Copy input fasta to tmp for faster access
  cp "$fasta_file" "$job_tmp/input.fa"

  run_log="$final_output_dir/run.log"
  {
    echo "[START] $(date '+%Y-%m-%d %H:%M:%S')"
    echo "TMPDIR: $TMPDIR"
    echo "Disk usage before run:"
    df -h "$TMPDIR"
    SECONDS=0

    singularity exec \
      --bind "$job_tmp:/tmp" \
      --bind "$job_tmp:/mnt" \
      --bind "$final_output_dir:/output" \
      "$singularity_image" run-cpgavas2 \
      -in "/mnt/input.fa" \
      -db 1 \
      -pid "$pid_name" \
      -out "/output"

    rc=$?
    duration=$SECONDS
    echo "[END]   $(date '+%Y-%m-%d %H:%M:%S') Duration: ${duration}s"
    echo "Disk usage after run:"
    df -h "$TMPDIR"

    if [ $rc -eq 0 ]; then
      echo "[DONE]  $base_name completed successfully."
    else
      echo "[FAIL]  $base_name failed."
    fi

    # Clean up tmp dir after run
    rm -rf "$job_tmp"

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
export singularity_image bind_dir final_dir TMP_BASE_DIR

find "$input_dir" -name "*.fa" | parallel --progress -j "$max_jobs" process_file {}
