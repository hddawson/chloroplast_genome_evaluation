#!/usr/bin/env Rscript

library(Biostrings)

aln_dir <- "data/tmp/alignedGenes"
out_dir <- "gene_trees"
dir.create(out_dir, showWarnings = FALSE)

aa_files <- list.files(aln_dir, pattern = "_AA_aligned\\.fasta$", full.names = TRUE)

for (f in aa_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(f))
  aln <- readAAStringSet(f)
  names(aln) <- sub("\\|.*$", "", names(aln))
  if (length(aln) < 4) next
  writeXStringSet(aln, file.path(out_dir, paste0(gene, ".faa")))
}

