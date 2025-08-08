import os
import subprocess
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

input_dir = Path("/workdir/hdd29/theRefseqening/data/genomes")
final_dir = Path("/workdir/hdd29/theRefseqening/data/results_batched")
singularity_image = "/local/workdir/hdd29/singularity_images/cpgavas2_sandbox"
batch_size = 100
max_workers = 5

final_dir.mkdir(parents=True, exist_ok=True)

def run_batch(batch_files, batch_id):
    with tempfile.TemporaryDirectory(prefix=f"{batch_id}_tmp_") as job_tmp:
        batch_file_path = Path(job_tmp) / f"batch_{batch_id}.txt"
        batch_file_path.write_text("\n".join(batch_files))
        run_log = final_dir / f"{batch_id}.log"

        cmd = [
            "singularity", "exec",
            "--bind", f"{job_tmp}:/mnt",
            "--bind", f"{input_dir}:/input",
            "--bind", f"{final_dir}:/output",
            "--bind", f"{batch_file_path}:/mnt/batch.txt",
            singularity_image,
            "bash", "-c",
            (
                'while read -r basename; do '
                'fasta="/input/${basename}.fa"; '
                '[ ! -f "$fasta" ] && echo "Missing: $fasta" && continue; '
                'pid=$(echo "$basename" | tr -d "_."); '
                'outdir="/output/$pid"; '
                'mkdir -p "$outdir"; '
                'cp "$fasta" /mnt/input.fa; '
                'run-cpgavas2 -in /mnt/input.fa -db 1 -pid "$pid" -out "$outdir"; '
                'find /mnt -mindepth 1 -maxdepth 1 ! -name "input.fa" -exec cp -r {} "$outdir/" \\;; '
                'rm -f /mnt/input.fa; '
                'done < /mnt/batch.txt'
            )
        ]

        with open(run_log, "w") as logf:
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
            proc.check_returncode()

def main(file_list_path):
    with open(file_list_path) as f:
        files = [line.strip() for line in f if line.strip()]
    batches = [files[i:i+batch_size] for i in range(0, len(files), batch_size)]

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for i, batch in enumerate(batches):
            executor.submit(run_batch, batch, f"batch_{i+1}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <file_list>")
        exit(1)
    main(sys.argv[1])
