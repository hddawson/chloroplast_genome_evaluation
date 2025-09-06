#!/usr/bin/env python3
import os
import glob
import subprocess
from tqdm import tqdm
from Bio import SeqIO

# --- config ---
PAL2NAL = "/programs/pal2nal/pal2nal.pl"
AA_DIR = "data/tmp/alignedGenes"
CDS_DIR = "data/tmp/genesToAlign"
OUT_DIR = "data/tmp/pal2nal"
FILTER_DIR = os.path.join(OUT_DIR, "filtered_inputs")
GENETIC_CODE = 11

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FILTER_DIR, exist_ok=True)


def run_pal2nal(aa_file, cds_file, gene, out_file):
    """Run pal2nal, capture stdout and stderr."""
    cmd = [
        "perl", PAL2NAL,
        aa_file, cds_file,
        "-codontable", str(GENETIC_CODE),
        "-output", "fasta"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    print("PRINTINF STDOURT!")
    print(stdout)

    # pal2nal prints errors to stdout
    if "#---  inconsistency between the following pep and nuc seqs" in stdout:
        bad_ids = []
        for line in stdout.splitlines():
            if line.startswith(">"):
                bad_ids.append(line[1:].strip().split()[0])  # take ID after ">"
        return False, bad_ids, stdout + stderr
    else:
        # success, write fasta
        with open(out_file, "w") as fout:
            fout.write(stdout)
        return True, [], ""


def filter_fastas(aa_file, cds_file, gene, bad_ids):
    """Remove bad IDs from fasta files and write filtered versions."""
    gene_dir = os.path.join(FILTER_DIR, gene)
    os.makedirs(gene_dir, exist_ok=True)

    def keep(rec): return rec.id not in bad_ids

    aa_records = [r for r in SeqIO.parse(aa_file, "fasta") if keep(r)]
    cds_records = [r for r in SeqIO.parse(cds_file, "fasta") if keep(r)]

    aa_out = os.path.join(gene_dir, f"{gene}_AA_aligned.filtered.fasta")
    cds_out = os.path.join(gene_dir, f"{gene}_CDS.filtered.fasta")
    SeqIO.write(aa_records, aa_out, "fasta")
    SeqIO.write(cds_records, cds_out, "fasta")

    return aa_out, cds_out


def main():
    aa_files = glob.glob(os.path.join(AA_DIR, "*_AA_aligned.fasta"))

    for aa_file in tqdm(aa_files, desc="Running pal2nal"):
        gene = os.path.basename(aa_file).replace("_AA_aligned.fasta", "")
        cds_file = os.path.join(CDS_DIR, f"{gene}_CDS.fasta")
        out_file = os.path.join(OUT_DIR, f"{gene}_pal2nal.fa")

        if not os.path.exists(cds_file):
            print(f"Skipping {gene}: missing CDS file", flush=True)
            continue

        print(f"Processing {gene}...")

        cur_aa, cur_cds = aa_file, cds_file
        attempt = 0
        while True:
            success, bad_ids, log = run_pal2nal(cur_aa, cur_cds, gene, out_file)
            if success:
                if attempt == 0:
                    print(f"{gene}: succeeded")
                else:
                    print(f"{gene}: succeeded after filtering")
                break
            elif bad_ids:
                attempt += 1
                print(f"{gene}: removing {len(bad_ids)} bad seq(s): {', '.join(bad_ids)}")
                cur_aa, cur_cds = filter_fastas(cur_aa, cur_cds, gene, bad_ids)
                # if everything is filtered, stop
                if not list(SeqIO.parse(cur_aa, "fasta")):
                    print(f"{gene}: all sequences filtered out, giving up")
                    break
            else:
                print(f"{gene}: unknown error, log follows:\n{log}")
                break

    print(f"Done. Results in {OUT_DIR}/")


if __name__ == "__main__":
    main()
