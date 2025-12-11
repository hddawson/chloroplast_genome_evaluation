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
