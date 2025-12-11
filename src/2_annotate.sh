#!/bin/bash
set -euo pipefail

mkdir -p /workdir/hdd29/tmp
export TMP_BASE_DIR="/workdir/hdd29/tmp"

# Set directories
input_dir="/workdir/hdd29/chloroplast_genome_evaluation/data/speciesWork/Salix/genomes"
final_dir="/workdir/hdd29/chloroplast_genome_evaluation/data/speciesWork/Salix/annotationResults/"
singularity_image="cpgavas2_0.03.sif"
bind_dir="$input_dir"

mkdir -p "$final_dir"
mkdir -p "$TMP_BASE_DIR"

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate singularity

# Check if file list is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <file_list>"
    echo "  file_list: Path to a file containing basenames (one per line) to process"
    echo "  Example: $0 data/smallFiles.txt"
    exit 1
fi

file_list="$1"

# Check if file list exists
if [ ! -f "$file_list" ]; then
    echo "Error: File list '$file_list' not found"
    exit 1
fi

# Calculate max parallel jobs based on free space in /tmp (assume 10GB per job)
available_kb=$(df --output=avail -k /workdir/hdd29/tmp | tail -1)
max_jobs=$(( available_kb / (10 * 1024 * 1024) ))  # 10GB per job in KB to GB
max_jobs=$(( max_jobs > 20 ? 20 : (max_jobs < 1 ? 1 : max_jobs) ))

echo "Detected free space: $(( available_kb / 1024 / 1024 )) GB"
echo "Setting max parallel jobs to: $max_jobs"
echo "Processing files from list: $file_list"

process_file() {
  basename="$1"
  fasta_file="$input_dir/${basename}"
  pid_name=$(echo "$basename" | tr -d '_.')
  final_output_dir="$final_dir/${pid_name}"

  # Check if input file exists
  if [ ! -f "$fasta_file" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $basename: Input file '$fasta_file' not found."
    return 1
  fi

  if [[ -d "$final_output_dir" && -n "$(ls -A "$final_output_dir" 2>/dev/null)" ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [SKIP] $basename: output already exists."
    return 0
  fi

  echo "$(date '+%Y-%m-%d %H:%M:%S') [START] $basename"
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

    /usr/bin/singularity exec \
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
      echo "[DONE]  $basename completed successfully."
    else
      echo "[FAIL]  $basename failed."
    fi

    cp -r "$job_tmp"/* "$final_output_dir/"

    # Clean up tmp dir after run
    rm -rf "$job_tmp"

    exit $rc
  } >> "$run_log" 2>&1

  if [ $? -eq 0 ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') [DONE] $basename"
  else
    echo "$(date '+%Y-%m-%d %H:%M:%S') [FAIL] $basename — check $run_log"
    return 1
  fi
}

export -f process_file
export singularity_image bind_dir final_dir TMP_BASE_DIR input_dir

# Read the file list and process each basename
cat "$file_list" | parallel --progress -j 10 process_file {}
