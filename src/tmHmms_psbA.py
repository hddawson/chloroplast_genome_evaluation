import biolib, glob, os
from Bio import SeqIO
from tqdm import tqdm

in_fasta = "data/AA_seqs/psbA_AA.fasta"
out_dir = "data/deeptmhmm_results_100seqsBatch"
os.makedirs(out_dir, exist_ok=True)

records = list(SeqIO.parse(in_fasta, "fasta"))
#records = records[:100]
batch_size = 100  # adjust to stay <50 MB
for i in tqdm(range(0, len(records), batch_size)):
    batch = records[i:i+batch_size]
    batch_fasta = f"{out_dir}/batch_{i//batch_size+1}.fasta"
    SeqIO.write(batch, batch_fasta, "fasta")
    print("Running", batch_fasta)
    deeptmhmm = biolib.load("DTU/DeepTMHMM")
    job = deeptmhmm.cli(args=f"--fasta {batch_fasta}")
    job.save_files(f"{out_dir}/batch_{i//batch_size+1}_results")