#!/usr/bin/env Rscript
# Map GWAS variants to phylogenetic tree and visualize

library(ape)
library(phytools)
library(ggplot2)
library(data.table)
library(Biostrings)
library(arrow)

data <- read_parquet("data/processed_data.parquet")
orders <- unique(data$Order)

# ---- LOAD TREE ----
tree_file <- "raxml_input/search1tree_seed2.raxml.bestTree"
stopifnot(file.exists(tree_file))
tree <- read.tree(tree_file)

is.rooted(tree)
is.monophyletic(tree,data$ID[grep("Pinales",data$Order)])

pinales <- intersect(tree$tip.label, data$ID[data$Order == "Pinales"])

tree <- root(tree,outgroup=pinales)

cat("Tree loaded:", length(tree$tip.label), "tips\n")

# ---- LOAD GWAS RESULTS ----
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", 
                          full.names = TRUE)
stopifnot(length(model_files) > 0)

all_results <- rbindlist(lapply(model_files, function(f) {
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P = m$P_aa_with_pcs,
      IDs = list(m$IDs)
    )
  }))
}))

threshold <- quantile(all_results$P, 0.05)
all_results[, is_hit := P < threshold]

cat("GWAS sites:", nrow(all_results), "\n")
cat("Significant sites:", sum(all_results$is_hit), "\n")

# ---- MAP VARIANTS TO TIPS ----
aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", 
                        full.names = TRUE)

tip_variant_counts <- data.table(tip = tree$tip.label, n_variants = 0, n_hits = 0)
setkey(tip_variant_counts, tip)

for (file in aln_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(file))
  gene_hits <- all_results[Gene == gene & is_hit == TRUE]
  
  if (nrow(gene_hits) == 0) next
  
  aln <- readAAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  aln_mat <- as.matrix(aln)
  
  for (hit in seq_len(nrow(gene_hits))) {
    pos <- gene_hits$Position[hit]
    
    # Check position is valid
    if (pos > ncol(aln_mat)) next
    
    tips_with_data <- names(aln)[aln_mat[, pos] != "-"]
    tips_in_tree <- intersect(tips_with_data, tree$tip.label)
    
    tip_variant_counts[tip %in% tips_in_tree, n_hits := n_hits + 1]
  }
}

cat("Tips with variant data:", sum(tip_variant_counts$n_hits > 0), "\n")

# ---- MAP VARIANTS TO NODES ----
# For each variant site, find the node where that variant arose
is.unro

node_variants <- data.table(node = (length(tree$tip.label) + 1):max(tree$edge), 
                            gene = NA_character_, 
                            position = NA_integer_)

for (file in aln_files) {
  gene <- sub("_AA_aligned\\.fasta", "", basename(file))
  gene_hits <- all_results[Gene == gene & is_hit == TRUE]
  
  if (nrow(gene_hits) == 0) next
  
  aln <- readAAStringSet(file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  aln_mat <- as.matrix(aln)
  
  for (hit in seq_len(nrow(gene_hits))) {
    pos <- gene_hits$Position[hit]
    if (pos > ncol(aln_mat)) next
    
    common_tips <- intersect(tree$tip.label, rownames(aln_mat))
    alleles <- aln_mat[common_tips, pos]
    names(alleles) <- common_tips
    
    # Find the MRCA of tips with each allele
    allele_states <- unique(alleles[alleles != "-"])
    for (allele in allele_states) {
      tips_with_allele <- tree$tip.label[alleles == allele & !is.na(alleles)]
      if (length(tips_with_allele) < 10) next
      print("tips with allele:")
      print(length(tips_with_allele))
      if (length(tips_with_allele) > 1) {
        mrca_node <- getMRCA(tree, tips_with_allele)
        node_variants <- rbind(node_variants, 
                               data.table(node = mrca_node, gene = gene, position = pos))
      }
    }
  }
}


plot(tree, show.tip.label = FALSE)
edgelabels(
  round(tree$edge.length, 3),
  cex = 1,
  frame = "none"
)

hist(tree$edge.length)
# Plot tree with variant nodes highlighted
pdf("results/tree_with_variant_nodes.pdf", width = 12, height = 16)

plot(tree, type = "f", show.tip.label = FALSE, 
     edge.width = 0.5, edge.color = "grey50")

# Highlight nodes with variants
variant_nodes <- unique(node_variants[!is.na(gene), node])
nodelabels(node = variant_nodes, pch = 16, cex = 0.5, col = "red")

title("Nodes where GWAS variants map (red dots)")

dev.off()
cat("Saved: results/tree_with_variant_nodes.pdf\n")

library(arrow)
data <- read_parquet("data/processed_data.parquet")
orders <- unique(data$Order)

for (ord in unique(data$Order)) {
  tips <- taxon_map[Order == ord, Tip]
  if (length(tips) < 2) next
  
  node <- getMRCA(tree, tips)
  
  p <- ggtree(tree) +
    geom_hilight(node=node, fill="red", alpha=0.3) +
    geom_tippoint(aes(subset=label %in% tips), color="red", size=1)
  
  ggsave(
    filename=paste0("variantMap_", ord, ".pdf"),
    plot=p,
    width=8,
    height=10
  )
}


library(phytools)

for (ord in unique(data$Order)) {
  ids <- data$ID[which(data$Order==ord)]
  tips <- intersect(tree$tip.label, ids)
  
  if (length(tips) < 10) next
  
  node <- getMRCA(tree, tips)
  if (is.null(node)) next
  
  print("ord")
  
  pdf(paste0("variantMap_", ord, ".pdf"), width = 8, height = 10)
  
  edge_cols <- rep("grey60", nrow(tree$edge))
  
  desc <- extract.clade(tree, node)$tip.label
  clade_tips <- match(desc, tree$tip.label)
  clade_edges <- which(tree$edge[,2] %in% c(clade_tips, getDescendants(tree, node)))
  
  edge_cols[clade_edges] <- "red"
  
  plot(
    tree,
    show.tip.label = FALSE,
    edge.color = edge_cols,
    edge.width = 0.6
  )
  
  tiplabels(
    pch = 16,
    col = "red",
    cex = 0.4,
    tip = match(tips, tree$tip.label)
  )
  
  
  title(ord)
  
  dev.off()
}

library(MASS)


hist(tree$edge.length)
x <- tree$edge.length
x <- x[x > 0]
fit <- fitdistr(x, "gamma")
fit
hist(x, breaks = 50, freq = FALSE, col = "grey80")
curve(dgamma(x, shape = fit$estimate["shape"], rate = fit$estimate["rate"]),
      add = TRUE, lwd = 2)


x <- tree$edge.length
x <- x[x > 0]

fit_g <- fitdistr(x, "gamma")
fit_ln <- fitdistr(x, "lognormal")
lambda <- mean(x)

hist(x, breaks = 50, freq = FALSE, col = "grey85", border = "black")

curve(dgamma(x, fit_g$estimate["shape"], fit_g$estimate["rate"]),
      add = TRUE, lwd = 2)

curve(dlnorm(x, fit_ln$estimate["meanlog"], fit_ln$estimate["sdlog"]),
      add = TRUE, lwd = 2, lty = 2)

legend(
  "topright",
  legend = c(
    paste0("Gamma (k=", round(fit_g$estimate["shape"],2), ")"),
    paste0("Log-normal (sd=", round(fit_ln$estimate["sdlog"],2), ")"),
    ""
  ),
  lwd = 2,
  lty = c(1,2,3),
  bty = "n"
)


x <- tree$edge.length
x <- x[x > 0]

fit_g <- fitdistr(x, "gamma")
fit_ln <- fitdistr(x, "lognormal")

d <- density(x)
y_g <- dgamma(d$x, fit_g$estimate["shape"], fit_g$estimate["rate"])
y_ln <- dlnorm(d$x, fit_ln$estimate["meanlog"], fit_ln$estimate["sdlog"])

cor_gamma <- cor(d$y, y_g)
cor_lognorm <- cor(d$y, y_ln)

c(gamma = cor_gamma, lognormal = cor_lognorm)


# ---- VISUALIZATION 1: Full tree with variant density ----
pdf("results/tree_variant_overview.pdf", width = 12, height = 16)
hist(tip_variant_counts$n_hits)
tip_colors <- ifelse(tip_variant_counts$n_hits == 0, "grey80",
                     ifelse(tip_variant_counts$n_hits < 400, "lightblue",
                            ifelse(tip_variant_counts$n_hits < 500, "steelblue", "darkblue")))

plot(tree, type = "fan", show.tip.label = FALSE, 
     edge.width = 0.5, edge.color = "grey50")
tiplabels(pch = 16, cex = 0.3, col = tip_colors)

legend("topleft", 
       legend = c("No variants", "1-4 variants", "5-9 variants", "10+ variants"),
       col = c("grey80", "lightblue", "steelblue", "darkblue"),
       pch = 16, bty = "n")
title("GWAS variant density across phylogeny")

dev.off()
cat("Saved: results/tree_variant_overview.pdf\n")

# ---- VISUALIZATION 2: Clade-level summaries ----
clades <- cutree(as.hclust(chronos(tree)), k = 20)
tip_variant_counts[, clade := clades[match(tip, names(clades))]]

clade_summary <- tip_variant_counts[, .(
  n_tips = .N,
  mean_variants = mean(n_hits),
  max_variants = max(n_hits),
  total_variants = sum(n_hits)
), by = clade]

setorder(clade_summary, -total_variants)

pdf("results/clade_variant_barplot.pdf", width = 10, height = 6)
par(mar = c(5, 4, 4, 2))
barplot(clade_summary$total_variants, 
        names.arg = paste0("C", clade_summary$clade),
        col = "steelblue",
        ylab = "Total variant sites",
        xlab = "Clade",
        main = "Variant distribution across major clades",
        las = 2)
dev.off()
cat("Saved: results/clade_variant_barplot.pdf\n")

# ---- VISUALIZATION 3: Subtree zoom on high-variant clade ----
top_clade <- clade_summary$clade[1]
top_clade_tips <- tip_variant_counts[clade == top_clade, tip]

if (length(top_clade_tips) > 2) {
  subtree <- keep.tip(tree, top_clade_tips)
  
  pdf("results/subtree_high_variant.pdf", width = 10, height = 12)
  
  sub_colors <- tip_colors[match(subtree$tip.label, tree$tip.label)]
  
  plot(subtree, cex = 0.5, edge.width = 1)
  tiplabels(pch = 16, cex = 0.8, col = sub_colors)
  title(paste0("High-variant clade (", length(top_clade_tips), " taxa)"))
  
  dev.off()
  cat("Saved: results/subtree_high_variant.pdf\n")
}

# ---- VISUALIZATION 4: Gene-specific variant mapping ----
gene_counts <- all_results[is_hit == TRUE, .N, by = Gene]
setorder(gene_counts, -N)
top_genes <- head(gene_counts$Gene, 5)

pdf("results/gene_specific_variants.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

for (gene in top_genes) {
  gene_tips <- data.table(tip = tree$tip.label, has_gene_variant = 0)
  
  gene_hits <- all_results[Gene == gene & is_hit == TRUE]
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  if (!file.exists(aln_file)) next
  
  aln <- readAAStringSet(aln_file)
  names(aln) <- sub("\\|.*$", "", names(aln))
  aln_mat <- as.matrix(aln)
  
  for (pos in gene_hits$Position) {
    if (pos > ncol(aln_mat)) next
    
    tips_with_data <- names(aln)[aln_mat[, pos] != "-"]
    tips_in_tree <- intersect(tips_with_data, tree$tip.label)
    gene_tips[tip %in% tips_in_tree, has_gene_variant := has_gene_variant + 1]
  }
  
  gene_colors <- ifelse(gene_tips$has_gene_variant == 0, "grey80", "red")
  
  plot(tree, type = "fan", show.tip.label = FALSE, 
       edge.width = 0.3, edge.color = "grey50")
  tiplabels(pch = 16, cex = 0.2, col = gene_colors)
  title(paste(gene, "-", sum(gene_tips$has_gene_variant > 0), "taxa"))
}

dev.off()
cat("Saved: results/gene_specific_variants.pdf\n")

# ---- SUMMARY STATS ----
cat("\n=== SUMMARY ===\n")
cat("Total tips:", length(tree$tip.label), "\n")
cat("Tips with variant data:", sum(tip_variant_counts$n_hits > 0), "\n")
cat("Mean variants per tip:", round(mean(tip_variant_counts$n_hits), 2), "\n")
cat("Max variants per tip:", max(tip_variant_counts$n_hits), "\n")
cat("\nTop clades by variant count:\n")
print(head(clade_summary, 10))

cat("\nDone! Check results/ for PDFs\n")