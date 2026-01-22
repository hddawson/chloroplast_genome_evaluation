import subprocess
from Bio import SeqIO
from collections import defaultdict
import glob
import re
import os
import random

aa = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")

def calc_pairwise_identity(aln_file, n_sample=1000):
    """Sample pairs from alignment, return list of identity scores"""
    recs = [r for r in SeqIO.parse(aln_file, "fasta") if aa.match(str(r.seq).replace("-", ""))]
    if len(recs) < 2:
        return [1.0]
    
    n_pairs = min(n_sample, len(recs) * (len(recs) - 1) // 2)
    identities = []
    
    for _ in range(n_pairs):
        r1, r2 = random.sample(recs, 2)
        s1, s2 = str(r1.seq), str(r2.seq)
        assert len(s1) == len(s2), f"Alignment length mismatch in {aln_file}"
        
        matches = sum(a == b and a != "-" for a, b in zip(s1, s2))
        aligned_pos = sum(a != "-" and b != "-" for a, b in zip(s1, s2))
        
        if aligned_pos > 0:
            identities.append(matches / aligned_pos)
    
    return identities

def get_threshold(identities, percentile=90):
    """Return identity threshold at given percentile"""
    sorted_ids = sorted(identities)
    idx = int(len(sorted_ids) * percentile / 100)
    return sorted_ids[min(idx, len(sorted_ids) - 1)]

def cluster_fasta(input_fasta, identity=0.9):
    """Run CD-HIT, return dict of cluster_id -> [seq_ids]"""
    out_prefix = input_fasta.replace(".fasta", "_cdhit")
    cmd = [
        "/programs/cd-hit-4.8.1/cd-hit", "-i", input_fasta, "-o", out_prefix,
        "-c", str(identity), "-n", "4" if identity >= 0.7 else "2",
        "-M", "0", "-T", "0", "-d", "0"
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    
    clusters = defaultdict(list)
    current_cluster = None
    with open(f"{out_prefix}.clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                current_cluster = int(line.strip().split()[1])
            else:
                match = re.search(r">(.+?)\.\.\.", line)
                assert match, f"Could not parse cluster line: {line}"
                seq_id = match.group(1)
                clusters[current_cluster].append(seq_id)
    
    assert len(clusters) > 0, f"No clusters found for {input_fasta}"
    return clusters

def sample_one_per_cluster(clusters, seqs):
    """Sample 1 sequence per cluster"""
    sampled = []
    for cl in clusters.values():
        for seq_id in cl:
            if seq_id in seqs:
                sampled.append(seqs[seq_id])
                break
    return sampled

# Setup
os.makedirs("data/netsurfp_input", exist_ok=True)
gene_files = glob.glob("data/AA_seqs/*fasta")
assert len(gene_files) > 0, "No fasta files found"

all_sampled = []

for f in gene_files:
    gene = os.path.basename(f).replace(".fasta", "")
    aln_file = f"data/tmp/alignedGenes/{gene}_aligned.fasta"
    
    if not os.path.exists(aln_file):
        print(f"{gene}: no alignment found, using default threshold 0.9")
        threshold = 0.9
    else:
        identities = calc_pairwise_identity(aln_file)
        threshold = get_threshold(identities, percentile=90)
        threshold = max(0.5, min(0.99, threshold))
    
    recs = [r for r in SeqIO.parse(f, "fasta") if aa.match(str(r.seq))]
    assert len(recs) > 0, f"No valid sequences in {f}"
    
    seqs = {r.id: r for r in recs}
    clusters = cluster_fasta(f, identity=threshold)
    
    sampled = sample_one_per_cluster(clusters, seqs)
    all_sampled.extend(sampled)
    print(f"{gene}: threshold={threshold:.2f}, {len(clusters)} clusters, sampled {len(sampled)}")

# Batch into ~2500 seq files
print(f"\nTotal sampled: {len(all_sampled)}")
batch_size = 2500
for i in range(0, len(all_sampled), batch_size):
    batch = all_sampled[i:i+batch_size]
    batch_num = i // batch_size + 1
    out_path = f"data/netsurfp_input/batch_{batch_num}.fasta"
    SeqIO.write(batch, out_path, "fasta")
    print(f"Batch {batch_num}: {len(batch)} sequences -> {out_path}")

    
"""old version, random sampling 
import random, glob, re
from Bio import SeqIO

aa = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
all_recs = []

for f in glob.glob("data/AA_seqs/*fasta"):
    recs = [r for r in SeqIO.parse(f,"fasta") if aa.match(str(r.seq))]
    all_recs.extend(random.sample(recs,100))

#how many total recordd?
print(len(all_recs))
#divide into 3 equal batches
batch_size = len(all_recs) // 3
batch1 = all_recs[:batch_size]
batch2 = all_recs[batch_size:2*batch_size]
batch3 = all_recs[2*batch_size:]


SeqIO.write(batch1, "data/netsurfp_batch_1.fasta", "fasta")
SeqIO.write(batch2, "data/netsurfp_batch_2.fasta", "fasta")
SeqIO.write(batch3, "data/netsurfp_batch_3.fasta", "fasta")
"""