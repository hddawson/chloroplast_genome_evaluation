#!/usr/bin/env Rscript
# Amino acid phylogeny input preparation: QC, tip collapsing, and partitioning

library(data.table)
library(Biostrings)
library(arrow)

# ---- PARAMETERS ----
MIN_NON_GAP <- 8500
AA_DISTANCE_THRESHOLD <- 30

# ---- LOAD DATA ----
cat("Loading data...\n")
data <- as.data.table(read_parquet("data/processed_data.parquet"))
embeds_with_mds <- readRDS("data/tmp/embeds_with_mds.rds")
setDT(embeds_with_mds)

stopifnot(all(c("ManualOutlier","Gene","ID") %in% colnames(embeds_with_mds)))
stopifnot("Order" %in% colnames(data))

# ---- CLEAN ID SET ----
clean_ids_by_gene <- embeds_with_mds[ManualOutlier == FALSE, .(ID, Gene)]
stopifnot(nrow(clean_ids_by_gene) > 0)
cat("Total clean species:", length(unique(clean_ids_by_gene$ID)), "\n")

# ---- PHENO IDS (to match GWAS) ----
pheno_col <- "pheno_wc2.1_2.5m_bio_8_p50"
pheno_ids <- data$ID[!is.na(data[[pheno_col]])]

# ---- BUILD AA SUPERALIGNMENT ----
cat("\nBuilding AA superalignment with GWAS-matching QC...\n")
aa_files <- list.files("data/tmp/alignedGenes/", pattern="_AA_aligned\\.fasta$", full.names=TRUE)
stopifnot(length(aa_files) > 0)

aligned_gene_list <- list()
for(file in aa_files) {
  gene <- sub("_AA_aligned\\.fasta","",basename(file))
  aln <- readAAStringSet(file)
  names(aln) <- sub("\\|.*$","",names(aln))
  
  # Intersection of clean IDs and phenotype IDs
  gene_clean_ids <- intersect(clean_ids_by_gene[Gene==gene, ID], pheno_ids)
  aln <- aln[names(aln) %in% gene_clean_ids]
  
  if(length(aln) > 0) aligned_gene_list[[gene]] <- aln
}
stopifnot(length(aligned_gene_list) > 0)

all_species <- Reduce(union, lapply(aligned_gene_list, names))
cat("Total species:", length(all_species), "\n")

# ---- BUILD SUPERMATRIX (GWAS-style QC) ----
super_list <- list()
kept_positions <- list()  # track original positions

for(gene in names(aligned_gene_list)) {
  aln <- aligned_gene_list[[gene]]
  aln_mat <- as.matrix(aln)
  
  gene_clean_ids <- intersect(names(aln), intersect(clean_ids_by_gene[Gene==gene, ID], pheno_ids))
  if(length(gene_clean_ids) < MIN_NON_GAP) next
  
  aln_mat <- aln_mat[match(gene_clean_ids, names(aln)), , drop=FALSE]
  
  keep_cols <- sapply(seq_len(ncol(aln_mat)), function(pos){
    residues <- aln_mat[,pos]
    non_gap <- residues != "-"
    if(sum(non_gap) < MIN_NON_GAP) return(FALSE)
    residue_table <- table(residues[non_gap])
    if(length(residue_table) < 2L) return(FALSE)
    TRUE
  })
  
  kept_mat <- aln_mat[,keep_cols, drop=FALSE]
  if(ncol(kept_mat) > 0){
    full_mat <- matrix("-", nrow=length(all_species), ncol=ncol(kept_mat),
                       dimnames=list(all_species,NULL))
    present <- match(rownames(kept_mat), all_species)
    full_mat[present[!is.na(present)], ] <- kept_mat
    super_list[[gene]] <- full_mat
    kept_positions[[gene]] <- which(keep_cols)  # original alignment positions
  }
}

stopifnot(length(super_list) > 0)
aa_supermat <- do.call(cbind, super_list)
aa_supermat_pre <- aa_supermat
cat("AA superalignment length:", ncol(aa_supermat), "sites\n")

# ---- BUILD GWAS-BASED PARTITION (PRE-COLLAPSE) ----
cat("\nBuilding GWAS-based partition...\n")
model_files <- list.files("results/residue_models_triple/", pattern="_effects\\.rds$", full.names=TRUE)
stopifnot(length(model_files) > 0)

gwas_results <- rbindlist(lapply(model_files, function(f){
  models <- readRDS(f)
  rbindlist(lapply(models, function(m){
    data.table(Gene=m$Gene, Position=m$Aligned_Position, P_aa_with_pcs=m$P_aa_with_pcs)
  }))
}))

gwas_thresh <- quantile(gwas_results$P_aa_with_pcs, 0.05, na.rm=TRUE)
gwas_results[, gwas_hit := P_aa_with_pcs < gwas_thresh]

# Build partition map
partition_map <- data.table()
global_pos <- 1
for(gene in names(super_list)){
  gene_len <- ncol(super_list[[gene]])
  orig_positions <- kept_positions[[gene]]
  stopifnot(length(orig_positions) == gene_len)
  
  partition_map <- rbind(partition_map,
                         data.table(Gene=gene, 
                                    GenePos=orig_positions,  # original alignment positions
                                    GlobalPos=global_pos:(global_pos+gene_len-1)))
  global_pos <- global_pos + gene_len
}
partition_map <- merge(partition_map, gwas_results,
                       by.x=c("Gene","GenePos"),
                       by.y=c("Gene","Position"), all.x=TRUE)
sum(is.na(partition_map$gwas_hit))
table(partition_map$gwas_hit)
partition_map[is.na(gwas_hit), gwas_hit := FALSE]
saveRDS(partition_map, file="raxml_input/partitionMap.rds")
quit()

# ---- COLLAPSE TIPS BY AA DISTANCE ----
cat("\nCollapsing tips by AA distance...\n")
genus_map <- data[ID %in% all_species, .(ID, Genus=sub(" .*$","",Organism))]
stopifnot(nrow(genus_map) == length(all_species))

set.seed(123)
collapse_groups <- list()
for(gen in unique(genus_map$Genus)){
  genus_ids <- genus_map[Genus==gen, ID]
  if(length(genus_ids)==1){ collapse_groups[[genus_ids]] <- genus_ids; next }
  remaining <- genus_ids
  while(length(remaining) > 0){
    rep <- sample(remaining,1)
    group <- rep
    to_check <- setdiff(remaining,rep)
    for(sp2 in to_check){
      seq1 <- aa_supermat_pre[rep,]; seq2 <- aa_supermat_pre[sp2,]
      diff <- sum(seq1 != seq2 & seq1 != "-" & seq2 != "-")
      if(diff <= AA_DISTANCE_THRESHOLD) group <- c(group, sp2)
    }
    collapse_groups[[rep]] <- group
    remaining <- setdiff(remaining, group)
  }
}
cat("Collapsed", length(all_species), "->", length(collapse_groups), "representatives\n")
reps <- names(collapse_groups)
aa_supermat_collapsed <- aa_supermat_pre[reps, , drop=FALSE]

# ---- WRITE AA SUPERALIGNMENT ----
dir.create("raxml_input", showWarnings=FALSE, recursive=TRUE)
aa_strings <- apply(aa_supermat_collapsed,1,paste0,collapse="")
writeXStringSet(AAStringSet(aa_strings), "raxml_input/superaa_collapsed.fasta")
cat("Wrote: raxml_input/superaa_collapsed.fasta\n")

# ---- WRITE PARTITION FILE ----
gwas_sites <- partition_map[gwas_hit==TRUE, GlobalPos]
background_sites <- partition_map[gwas_hit==FALSE, GlobalPos]
stopifnot(length(gwas_sites) + length(background_sites) == ncol(aa_supermat_pre))

cat("AA, GWAS =", paste(gwas_sites, collapse=","), "\n",
    file="raxml_input/gwas_partition.txt")
cat("AA, Background =", paste(background_sites, collapse=","), "\n",
    file="raxml_input/gwas_partition.txt", append=TRUE)
cat("Wrote: raxml_input/gwas_partition.txt\n")
cat("GWAS sites:", length(gwas_sites), "\nBackground sites:", length(background_sites), "\n")
cat("\nPipeline complete.\n")
nrow(aa_supermat_collapsed)

# ---- cds ---- 
# ---- BUILD CDS SUPERALIGNMENT (matching collapsed AA species) ----
cat("\nBuilding CDS superalignment for collapsed species...\n")

cds_files <- list.files("data/tmp/alignedGenes/", pattern="_CDS_aligned\\.fasta$", full.names=TRUE)
stopifnot(length(cds_files) > 0)

# Only use genes that passed AA QC
genes_with_aa <- names(super_list)

cds_super_list <- list()
for(file in cds_files) {
  gene <- sub("_CDS_aligned\\.fasta","",basename(file))
  if(!gene %in% genes_with_aa) next
  
  cds_aln <- readDNAStringSet(file)
  names(cds_aln) <- sub("\\|.*$","",names(cds_aln))
  cds_aln <- cds_aln[names(cds_aln) %in% reps]
  if(length(cds_aln) == 0) next
  
  cds_mat <- as.matrix(cds_aln)
  
  # Get AA positions that passed QC for this gene
  aa_kept <- kept_positions[[gene]]
  # Convert to codon positions (each AA = 3 nucleotides)
  cds_kept <- as.vector(sapply(aa_kept, function(p) (p-1)*3 + 1:3))
  cds_kept <- cds_kept[cds_kept <= ncol(cds_mat)]
  
  if(length(cds_kept) == 0) next
  cds_mat <- cds_mat[, cds_kept, drop=FALSE]
  
  full_cds_mat <- matrix("-", nrow=length(reps), ncol=ncol(cds_mat),
                         dimnames=list(reps, NULL))
  present <- match(rownames(cds_mat), reps)
  full_cds_mat[present[!is.na(present)], ] <- cds_mat[!is.na(present), ]
  
  cds_super_list[[gene]] <- full_cds_mat
}

cds_supermat <- do.call(cbind, cds_super_list)
cat("CDS superalignment:", nrow(cds_supermat), "taxa,", ncol(cds_supermat), "sites\n")

stopifnot(length(cds_super_list) > 0)
cat("CDS genes included:", length(cds_super_list), "\n")

cds_supermat <- do.call(cbind, cds_super_list)
cat("CDS superalignment:", nrow(cds_supermat), "taxa,", ncol(cds_supermat), "sites\n")

# ---- WRITE CDS SUPERALIGNMENT ----
cds_strings <- apply(cds_supermat, 1, paste0, collapse="")
writeXStringSet(DNAStringSet(cds_strings), "raxml_input/supercds_collapsed.fasta")
cat("Wrote: raxml_input/supercds_collapsed.fasta\n")

# ---- CHECK ----
# CDS should be ~3x AA length (allowing for some variation due to filtering)
cat("AA sites:", ncol(aa_supermat_collapsed), "\n")
cat("CDS sites:", ncol(cds_supermat), "\n")
cat("Ratio:", ncol(cds_supermat) / ncol(aa_supermat_collapsed), "\n")


# ---- commands ---- 

#cd raxml_input

#/programs/raxml-ng_v1.2.0/raxml-ng --check --msa superaa_collapsed.fasta --model LG+G4 --prefix check


# /programs/raxml-ng_v1.2.0/raxml-ng --search --msa superaa_collapsed.fasta --model LG+G4 --prefix aa_tree_fast --threads auto

# If check passes, run the tree
#/programs/raxml-ng_v1.2.0/raxml-ng --all --msa superaa_collapsed.fasta --model LG+G4 --prefix aa_tree --threads auto --bs-trees 100

#cds

# Check alignment
#/programs/raxml-ng_v1.2.0/raxml-ng --check --msa supercds_collapsed.fasta --model GTR+G4 --prefix check_cds

# Fast tree search
#/programs/raxml-ng_v1.2.0/raxml-ng --search --msa supercds_collapsed.fasta --model GTR+G4 --prefix cds_tree_fast --threads auto


get_majority_aa <- function(residues) {
  residues <- residues[residues != "-"]
  if(length(residues) == 0) return(NA_character_)
  tt <- table(residues)
  names(tt)[which.max(tt)]
}

gwas_majority <- sapply(gwas_sites, function(i) get_majority_aa(aa_supermat_collapsed[,i]))
bg_majority <- sapply(background_sites, function(i) get_majority_aa(aa_supermat_collapsed[,i]))

gwas_majority <- gwas_majority[!is.na(gwas_majority)]
bg_majority <- bg_majority[!is.na(bg_majority)]

gwas_freq <- table(gwas_majority)
bg_freq <- table(bg_majority)

all_aa <- union(names(gwas_freq), names(bg_freq))
gwas_freq <- gwas_freq[all_aa]; gwas_freq[is.na(gwas_freq)] <- 0
bg_freq <- bg_freq[all_aa]; bg_freq[is.na(bg_freq)] <- 0

# Chi-square test
chisq.test(rbind(gwas_freq, bg_freq))

# Look at enrichment/depletion per AA
gwas_prop <- gwas_freq / sum(gwas_freq)
bg_prop <- bg_freq / sum(bg_freq)
log2_enrichment <- log2((gwas_prop + 0.001) / (bg_prop + 0.001))
sort(log2_enrichment)

table(partition_map[gwas_hit == TRUE, Gene])

# ---- effect ananlysis wrt residue count per sit e
effect_dt <- rbindlist(lapply(model_files, function(f){
  models <- readRDS(f)
  rbindlist(lapply(models, function(m){
    if(is.null(m$effects)) return(NULL)
    counts <- unlist(m$residue_counts)
    dt <- m$effects[, .(Gene = m$Gene, Position = m$Aligned_Position, 
                        Residue, Effect, P_value)]
    # Match residue to count
    dt[, Residue_clean := gsub("^X_aa", "", Residue)]
    dt[, count := counts[Residue_clean], by = .I]
    dt[count >= 10]
  }))
}))
stopifnot(nrow(effect_dt) > 0)

effect_dt[, Residue := Residue_clean]
effect_dt[, Residue_clean := NULL]

# Merge gwas_hit
effect_dt <- merge(effect_dt, 
                   gwas_results[, .(Gene, Position, gwas_hit)],
                   by = c("Gene", "Position"))

# Remove ambiguous codes
effect_dt_clean <- effect_dt[!Residue %in% c("-", "X", "B", "J")]

cat("Retained", nrow(effect_dt_clean), "effects with count >= 10\n")

# Compare effect magnitudes by AA at GWAS vs background sites
effect_summary <- effect_dt_clean[, .(
  mean_effect = mean(Effect),
  mean_abs_effect = mean(abs(Effect)),
  n = .N
), by = .(Residue, gwas_hit)]

# Wide format for comparison
effect_wide <- dcast(effect_summary, Residue ~ gwas_hit, 
                     value.var = c("mean_effect", "mean_abs_effect", "n"))
setnames(effect_wide, c("Residue", "mean_eff_bg", "mean_eff_gwas", 
                        "abs_eff_bg", "abs_eff_gwas", "n_bg", "n_gwas"))
effect_wide[order(-abs_eff_gwas)]

wilcox.test(abs(Effect) ~ gwas_hit, data = effect_dt_clean)

effect_dt_clean[, .(mean_abs = mean(abs(Effect))), by = gwas_hit]

effect_wide[, effect_ratio := abs_eff_gwas / abs_eff_bg]
effect_wide[!Residue %in% c("-", "X", "B", "J")][order(-effect_ratio)]

effect_dt_clean[gwas_hit == TRUE, .(
  mean_eff = mean(Effect),
  se = sd(Effect)/sqrt(.N),
  t = mean(Effect) / (sd(Effect)/sqrt(.N)),
  n = .N
), by = Residue][order(mean_eff)]

# Add biochemical properties
aa_properties <- data.table(
  Residue = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
  Property = c("Hydrophobic","Cys","Charged_neg","Charged_neg","Aromatic","Small",
               "Charged_pos","Hydrophobic","Charged_pos","Hydrophobic","Hydrophobic",
               "Polar","Proline","Polar","Charged_pos","Polar","Polar","Hydrophobic",
               "Aromatic","Aromatic"),
  Hydropathy = c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, 
                 -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3)  # Kyte-Doolittle
)

effect_dt_clean <- merge(effect_dt_clean, aa_properties, by = "Residue", all.x = TRUE)

# Summary tables
effect_summary <- effect_dt_clean[, .(
  mean_eff = mean(Effect),
  mean_abs_eff = mean(abs(Effect)),
  se = sd(Effect)/sqrt(.N),
  n = .N
), by = .(Residue, Property, Hydropathy, gwas_hit)]

# Save
saveRDS(list(
  effect_dt = effect_dt_clean,
  effect_summary = effect_summary,
  effect_wide = effect_wide
), "results/aa_enrichment_analysis.rds")

# Property-level summary
prop_summary <- effect_dt_clean[gwas_hit == TRUE, .(
  mean_eff = mean(Effect),
  se = sd(Effect)/sqrt(.N),
  n = .N
), by = Property][order(mean_eff)]
print(prop_summary)

# Correlation of hydropathy with effect at GWAS sites
gwas_by_aa <- effect_dt_clean[gwas_hit == TRUE, .(mean_eff = mean(Effect)), by = .(Residue, Hydropathy)]
cor.test(gwas_by_aa$Hydropathy, gwas_by_aa$mean_eff)


# Save analysis results
analysis_results <- list(
  effect_dt = effect_dt_clean,
  effect_summary = effect_wide[!Residue %in% c("-", "X", "B", "J")],
  gwas_vs_bg_test = wilcox.test(abs(Effect) ~ gwas_hit, data = effect_dt_clean),
  mean_abs_by_class = effect_dt_clean[, .(mean_abs = mean(abs(Effect)), n = .N), by = gwas_hit],
  directional_effects = effect_dt_clean[gwas_hit == TRUE, .(
    mean_eff = mean(Effect),
    se = sd(Effect)/sqrt(.N),
    n = .N
  ), by = Residue][order(mean_eff)]
)

# Set up output
pdf("results/aa_effect_plots.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

# Plot 1: Effect magnitude at GWAS vs background
boxplot(abs(Effect) ~ gwas_hit, data = effect_dt_clean,
        names = c("Background", "GWAS"),
        ylab = "|Effect|",
        main = "Effect magnitude by site class",
        col = c("gray80", "steelblue"))

# Plot 2: Effect ratio by AA
effect_wide_clean <- effect_wide[!Residue %in% c("-", "X", "B", "J")]
effect_wide_clean <- effect_wide_clean[order(effect_ratio)]
barplot(effect_wide_clean$effect_ratio, 
        names.arg = effect_wide_clean$Residue,
        las = 2,
        ylab = "Effect ratio (GWAS / Background)",
        main = "Effect amplification at GWAS sites",
        col = "steelblue")
abline(h = 1, lty = 2)

# Plot 3: Hydropathy vs mean effect at GWAS sites
gwas_by_aa <- effect_dt_clean[gwas_hit == TRUE, .(
  mean_eff = mean(Effect),
  Hydropathy = Hydropathy[1]
), by = Residue]
plot(gwas_by_aa$Hydropathy, gwas_by_aa$mean_eff,
     xlab = "Hydropathy (Kyte-Doolittle)",
     ylab = "Mean effect at GWAS sites",
     main = "Hydropathy vs climate effect",
     pch = 19, col = "steelblue")
text(gwas_by_aa$Hydropathy, gwas_by_aa$mean_eff, 
     labels = gwas_by_aa$Residue, pos = 3, cex = 0.8)
abline(h = 0, lty = 2)
abline(lm(mean_eff ~ Hydropathy, data = gwas_by_aa), col = "red", lty = 2)

dev.off()
cat("Wrote: results/aa_effect_plots.pdf\n")


# Mean effect by AA at GWAS sites with error bars
gwas_by_aa <- effect_dt_clean[gwas_hit == TRUE, .(
  mean_eff = mean(Effect),
  se = sd(Effect)/sqrt(.N),
  n = .N
), by = Residue][order(mean_eff)]

plot(effect_dt_clean$Effect, -log10(effect_dt_clean$P_value))
pdf("results/aa_mean_effect_gwas.pdf", width = 8, height = 5)
par(mar = c(5, 4, 3, 1))

x <- barplot(gwas_by_aa$mean_eff,
             names.arg = gwas_by_aa$Residue,
             las = 2,
             ylab = "Mean effect at GWAS sites",
             main = "Amino acid effects on bio8",
             col = ifelse(gwas_by_aa$mean_eff < 0, "steelblue", "coral"),
             ylim = range(c(gwas_by_aa$mean_eff - gwas_by_aa$se, 
                            gwas_by_aa$mean_eff + gwas_by_aa$se)) * 1.1)

arrows(x, gwas_by_aa$mean_eff - gwas_by_aa$se,
       x, gwas_by_aa$mean_eff + gwas_by_aa$se,
       angle = 90, code = 3, length = 0.05)

abline(h = 0, lty = 2)

dev.off()
cat("Wrote: results/aa_mean_effect_gwas.pdf\n")

saveRDS(analysis_results, "results/aa_enrichment_analysis.rds")

# Classify AAs by biochemical property
aa_properties <- data.table(
  Residue = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
  Property = c("Hydrophobic","Cys","Neg_charged","Neg_charged","Aromatic","Small",
               "Pos_charged","Hydrophobic","Pos_charged","Hydrophobic","Hydrophobic",
               "Polar","Proline","Polar","Pos_charged","Polar","Polar","Hydrophobic",
               "Aromatic","Aromatic")
)
effect_dt_clean <- merge(effect_dt_clean, aa_properties, by = "Residue", all.x = TRUE)
saveRDS(effect_dt_clean, "results/effect_dt_with_properties.rds")

chisq_result <- chisq.test(rbind(gwas_freq, bg_freq))
print(chisq_result)

# Per-AA enrichment with fisher test for each
aa_enrichment <- data.table(
  AA = all_aa,
  gwas = as.numeric(gwas_freq[all_aa]),
  bg = as.numeric(bg_freq[all_aa])
)
aa_enrichment[, gwas_prop := gwas / sum(gwas)]
aa_enrichment[, bg_prop := bg / sum(bg)]
aa_enrichment[, log2_enrich := log2((gwas_prop + 0.001) / (bg_prop + 0.001))]

# Fisher test per AA
aa_enrichment[, fisher_p := sapply(seq_len(.N), function(i) {
  mat <- matrix(c(gwas[i], sum(gwas) - gwas[i],
                  bg[i], sum(bg) - bg[i]), nrow = 2)
  fisher.test(mat)$p.value
})]
aa_enrichment[, fisher_padj := p.adjust(fisher_p, method = "BH")]
aa_enrichment <- aa_enrichment[order(-log2_enrich)]
print(aa_enrichment)

# Plot
pdf("results/aa_enrichment_gwas_sites.pdf", width = 8, height = 5)
par(mar = c(5, 4, 3, 1))

aa_enrichment <- aa_enrichment[order(log2_enrich)]
cols <- ifelse(aa_enrichment$fisher_padj < 0.05, 
               ifelse(aa_enrichment$log2_enrich > 0, "coral", "steelblue"),
               "gray70")

x <- barplot(aa_enrichment$log2_enrich,
             names.arg = aa_enrichment$AA,
             las = 2,
             ylab = expression(log[2]~"enrichment (GWAS / Background)"),
             main = "Amino acid enrichment at GWAS sites",
             col = cols)
abline(h = 0, lty = 2)
legend("topleft", legend = c("Enriched (FDR < 0.05)", "Depleted (FDR < 0.05)", "NS"),
       fill = c("coral", "steelblue", "gray70"), bty = "n", cex = 0.8)

dev.off()
cat("Wrote: results/aa_enrichment_gwas_sites.pdf\n")

# Print significant
cat("\nSignificant enrichment/depletion (FDR < 0.05):\n")
print(aa_enrichment[fisher_padj < 0.05])
