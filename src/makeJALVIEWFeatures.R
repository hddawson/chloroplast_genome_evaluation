library(data.table)

# Load mapping
map <- fread("data/tmp/aa_supermatrix.fasta_mapping.csv")

# Jalview uses 1-based coords
map[, site := site + 1]

# Collapse consecutive runs of the same gene
map[, run_id := rleid(gene)]
collapsed <- map[, .(start = min(site), end = max(site)), by = .(gene, run_id)]

# Build feature table
features <- data.table(
  SeqID = "*",      # use "alignment" if feature applies to all seqs
  start = collapsed$start,
  end   = collapsed$end,
  type  = "gene",
  label = collapsed$gene
)

# Save as TSV
fwrite(features, "data/tmp/aa_supermatrix_features.tsv", sep="\t", col.names=FALSE)
