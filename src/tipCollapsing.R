#!/usr/bin/env Rscript
# RAxML-NG phylogeny input preparation pipeline
# Handles QC, tip collapsing, and constraint tree generation

library(data.table)
library(Biostrings)
library(ape)
library(arrow)

# PARAMETERS
MIN_NON_GAP <- 8500
AA_DISTANCE_THRESHOLD <- 30  # amino acids

# ---- LOAD DATA ----
cat("Loading data...\n")
data <- as.data.table(read_parquet("data/processed_data.parquet"))
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds_with_mds)

stopifnot("ManualOutlier" %in% colnames(embeds_with_mds))
stopifnot("Gene" %in% colnames(embeds_with_mds))
stopifnot("ID" %in% colnames(embeds_with_mds))
stopifnot("Order" %in% colnames(data))

# Get clean IDs (exclude manual outliers)
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)

cat("Total clean species:", length(unique(clean_ids_by_gene$ID)), "\n")

# ---- PART 1: BUILD CDS SUPERALIGNMENT ----
cat("\nBuilding CDS superalignment...\n")
aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_CDS_aligned\\.fasta$", 
                        full.names = TRUE)
stopifnot(length(aln_files) > 0)

get_gene <- function(path) sub("_CDS_aligned\\.fasta", "", basename(path))

aligned_gene_list <- list()
for (file in aln_files) {
  gene <- get_gene(file)
  aln <- readDNAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  
  gene_clean_ids <- clean_ids_by_gene[Gene == gene, ID]
  aln <- aln[names(aln) %in% gene_clean_ids]
  
  if (length(aln) > 0) {
    aligned_gene_list[[gene]] <- aln
  }
}

stopifnot(length(aligned_gene_list) > 0)
all_species <- Reduce(union, lapply(aligned_gene_list, names))
cat("Total species:", length(all_species), "\n")

# Build supermatrix
super_list <- list()
for (gene in names(aligned_gene_list)) {
  aln <- aligned_gene_list[[gene]]
  aln_mat <- as.matrix(aln)
  gene_len <- ncol(aln_mat)
  
  mat <- matrix("-", nrow = length(all_species), ncol = gene_len,
                dimnames = list(all_species, NULL))
  present <- match(names(aln), all_species)
  mat[present[!is.na(present)], ] <- aln_mat
  
  # QC columns
  keep <- sapply(seq_len(gene_len), function(j) {
    col <- mat[, j]
    non_gap <- col[col != "-"]
    length(non_gap) >= MIN_NON_GAP && length(unique(non_gap)) >= 2
  })
  
  kept_mat <- mat[, keep, drop = FALSE]
  if (ncol(kept_mat) > 0) {
    super_list[[gene]] <- kept_mat
  }
}

stopifnot(length(super_list) > 0)
supermat <- do.call(cbind, super_list)
cat("Superalignment length:", ncol(supermat), "bp\n")

# ---- PART 2: COLLAPSE TIPS BY AA DISTANCE ----
cat("\nCollapsing tips by AA distance...\n")
aa_files <- list.files("data/tmp/alignedGenes/", 
                       pattern = "_AA_aligned\\.fasta$", 
                       full.names = TRUE)
stopifnot(length(aa_files) > 0)

# Build combined AA superalignment
aa_super_list <- list()
for (file in aa_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(file))
  if (!gene %in% names(super_list)) next
  
  aln <- readAAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  aln <- aln[names(aln) %in% all_species]
  
  if (length(aln) > 0) {
    mat <- as.matrix(aln)
    full_mat <- matrix("-", nrow = length(all_species), ncol = ncol(mat),
                       dimnames = list(all_species, NULL))
    present <- match(names(aln), all_species)
    full_mat[present[!is.na(present)], ] <- mat
    aa_super_list[[gene]] <- full_mat
  }
}

for (gene in names(aa_super_list)) {
  gene_clean_ids <- clean_ids_by_gene[Gene == gene, ID]
  mat <- aa_super_list[[gene]]
  keep_rows <- rownames(mat) %in% gene_clean_ids
  mat[!keep_rows, ] <- "-"  # Set non-clean rows to gaps
  aa_super_list[[gene]] <- mat
}

stopifnot(length(aa_super_list) > 0)
aa_supermat <- do.call(cbind, aa_super_list)
cat("Combined AA alignment length:", ncol(aa_supermat), "sites\n")

# Extract genus from Organism
setDT(data)
genus_map <- data[ID %in% all_species, .(ID, Genus = sub(" .*$", "", Organism))]
stopifnot(nrow(genus_map) == length(all_species))

set.seed(123)
# Collapse within genera
collapse_groups <- list()
for (gen in unique(genus_map$Genus)) {
  genus_ids <- genus_map[Genus == gen, ID]
  if (length(genus_ids) == 1) {
    collapse_groups[[genus_ids]] <- genus_ids
    next
  }
  
  remaining <- genus_ids
  while (length(remaining) > 0) {
    rep <- sample(remaining, 1)
    group <- rep
    
    to_check <- setdiff(remaining, rep)
    for (sp2 in to_check) {
      # Count AA differences across full alignment
      seq1 <- aa_supermat[rep, ]
      seq2 <- aa_supermat[sp2, ]
      diff <- sum(seq1 != seq2 & seq1 != "-" & seq2 != "-")
      
      if (diff <= AA_DISTANCE_THRESHOLD) {
        group <- c(group, sp2)
      }
    }
    
    collapse_groups[[rep]] <- group
    remaining <- setdiff(remaining, group)
  }
}

cat("Collapsed", length(all_species), "->", length(collapse_groups), "representatives\n")

# Keep one representative per group
reps <- names(collapse_groups)
supermat_collapsed <- supermat[reps, , drop = FALSE]

# ---- WRITE AA SUPERALIGNMENT (COLLAPSED) ----
aa_supermat_collapsed <- aa_supermat[reps, , drop = FALSE]

aa_strings <- apply(aa_supermat_collapsed, 1, paste0, collapse = "")
aa_out <- "raxml_input/superaa_collapsed.fasta"

dir.create("raxml_input", showWarnings = FALSE, recursive = TRUE)
writeXStringSet(AAStringSet(aa_strings), aa_out)

cat("Wrote:", aa_out, "\n")


# Write superalignment
dir.create("raxml_input", showWarnings = FALSE, recursive = TRUE)
out_fasta <- "raxml_input/supercds_collapsed.fasta"
dna_strings <- apply(supermat_collapsed, 1, paste0, collapse = "")
writeXStringSet(DNAStringSet(dna_strings), out_fasta)
cat("Wrote:", out_fasta, "\n")

# ---- PART 3: BUILD CONSTRAINT TREE ----
cat("\nBuilding constraint tree...\n")
backbone <- read.tree("data/2_global_order_level.tre")
backbone$edge.length[is.na(backbone$edge.length)] <- 1

# Map representatives to orders
rep_orders <- data[ID %in% reps, .(ID, Order)]
stopifnot(nrow(rep_orders) == length(reps))

groups <- split(rep_orders$ID, rep_orders$Order)
backbone <- drop.tip(backbone, setdiff(backbone$tip.label, names(groups)))
plot(backbone)
for (ord in names(groups)) {
  sub <- stree(length(groups[[ord]]), tip.label = groups[[ord]])
  backbone <- bind.tree(backbone, sub, where = which(backbone$tip.label == ord))
}

# Validate constraint tree
stopifnot(length(backbone$tip.label) == length(unique(backbone$tip.label)))
stopifnot(all(backbone$tip.label %in% reps))
stopifnot(all(sapply(groups, function(x) is.monophyletic(backbone, x))))

write.tree(backbone, "raxml_input/constraint_tree.tre")
cat("Wrote: raxml_input/constraint_tree.tre\n")

# ---- PART 4: BUILD PARTITION FILE FOR RAXML ----
cat("\nBuilding GWAS-based partition file...\n")

# Load GWAS results
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_with_pcs = m$P_aa_with_pcs
    )
  }))
}))

# Define GWAS hits
gwas_thresh <- quantile(gwas_results$P_aa_with_pcs, 0.05)
gwas_results[, gwas_hit := P_aa_with_pcs < gwas_thresh]

# Map gene positions to global superalignment positions
partition_map <- data.table()
global_pos <- 1

for (gene in names(super_list)) {
  gene_mat <- super_list[[gene]]
  gene_len <- ncol(gene_mat)
  
  partition_map <- rbind(partition_map, data.table(
    Gene = gene,
    GenePos = seq_len(gene_len),
    GlobalPos = seq(global_pos, global_pos + gene_len - 1)
  ))
  
  global_pos <- global_pos + gene_len
}

# Merge with GWAS status
partition_map <- merge(partition_map, gwas_results, 
                       by.x = c("Gene", "GenePos"), 
                       by.y = c("Gene", "Position"),
                       all.x = TRUE)
partition_map[is.na(gwas_hit), gwas_hit := FALSE]

# Build RAxML partition file
partition_file <- "raxml_input/gwas_partition.txt"

gwas_sites <- partition_map[gwas_hit == TRUE, GlobalPos]
background_sites <- partition_map[gwas_hit == FALSE, GlobalPos]

cat("AA, GWAS =", paste(gwas_sites, collapse = ","), "\n", 
    file = partition_file)
cat("AA, Background =", paste(background_sites, collapse = ","), "\n", 
    file = partition_file, append = TRUE)

cat("Wrote:", partition_file, "\n")
cat("GWAS sites:", length(gwas_sites), "\n")
cat("Background sites:", length(background_sites), "\n")

# RAxML command with partition
raxml_partition_cmd <- sprintf(
  "/programs/raxml-ng_v1.2.0/raxml-ng --all --msa %s --model GTR+G --tree-constraint %s --bs-trees 100 --threads 20 --subs-rates %s --prefix gwas_partition",
  "raxml_input/supercds_collapsed.fasta",
  "raxml_input/constraint_tree.tre",
  partition_file
)

cat("\nRAxML command with partition:\n")
cat(raxml_partition_cmd, "\n")
writeLines(raxml_partition_cmd, "raxml_input/run_raxml_partition.sh")

cat("\nPipeline complete. Run RAxML-NG with:\n")
cat("/programs/raxml-ng_v1.2.0/raxml-ng --all --msa raxml_input/supercds_collapsed.fasta --model GTR+G",
    "--tree-constraint raxml_input/constraint_tree.tre --bs-trees 100 --threads 40\n")

cmd <- "/programs/raxml-ng_v1.2.0/raxml-ng --all --msa supercds_collapsed.fasta --model GTR+G
--tree-constraint constraint_tree.tre --bs-trees 100 --threads 20"

"/programs/raxml-ng_v1.2.0/raxml-ng --check --msa supercds_collapsed.fasta --model GTR+G --prefix T1"
"/programs/raxml-ng_v1.2.0/raxml-ng --parse --msa supercds_collapsed.fasta --model GTR+G --prefix T2"

"/programs/raxml-ng_v1.2.0/raxml-ng --search1 --msa supercds_collapsed.fasta --model GTR+G --prefix search1tree --threads 11 --seed 1103"

"/programs/raxml-ng_v1.2.0/raxml-ng \
  --all \
  --msa supercds_collapsed.fasta \
  --model GTR+G \
  --tree-constraint constraint_tree.tre \
  --bs-trees 100 \
  --threads 11 \
  --prefix supercds_run1 \
  --checkpoint \
  --checkpoint-interval 60"