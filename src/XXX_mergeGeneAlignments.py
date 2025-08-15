import os, glob, argparse, pandas as pd, numpy as np
from Bio import SeqIO
from tqdm.auto import tqdm

ambig = str.maketrans({c: '-' for c in 'WYKSMNR'})

def aln_df(path):
    recs = list(SeqIO.parse(path, 'fasta'))
    L = len(recs[0].seq)
    assert all(len(r.seq) == L for r in recs)
    ids = [r.id.split('|')[0] for r in recs]
    seqs = [str(r.seq).upper().translate(ambig) for r in recs]
    A = np.array([list(s) for s in seqs], dtype='<U1')

    # Filter columns within this gene alignment
    gap_counts = (A == '-').sum(axis=0)
    max_gaps = int(np.floor(0.05 * A.shape[0]))  # 5% threshold
    keep_mask = gap_counts <= max_gaps
    A = A[:, keep_mask]

    gene = os.path.splitext(os.path.basename(path))[0].split('_')[0]
    cols = [f'{gene}_alnSite_{i}' for i in range(keep_mask.sum())]
    return pd.DataFrame(A, index=ids, columns=cols)

def build_supermatrix(aln_dir, pattern='*.fasta'):
    files = sorted(glob.glob(os.path.join(aln_dir, pattern)))
    dfs = [aln_df(f) for f in tqdm(files)]
    M = pd.concat(dfs, axis=1).sort_index()
    return M.fillna('-')

def gap_counts(df):
    return (df == '-').sum(axis=0)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--aln_dir', required=True)
    p.add_argument('--pattern', default='*.fasta')
    p.add_argument('--out_parquet', required=True)
    p.add_argument('--out_gapcounts', required=True)
    a = p.parse_args()

    supermatrix = build_supermatrix(a.aln_dir, a.pattern)
    print("alignments merged, saving")
    supermatrix.to_parquet(a.out_parquet)
    gap_counts(supermatrix).to_frame('gaps').to_csv(a.out_gapcounts)
