library(data.table)
library(Biostrings)
library(Peptides)

# Load PCA results
pca <- readRDS("data/tmp/pca_results_cleaned.rds")
clean_embeds <- readRDS("data/tmp/embeds_with_mds.rds")
clean_embeds <- clean_embeds[abs(clean_embeds$MDS_zscore) < 3]  # Remove outliers

# Get sequence files
seq_files <- list.files("data/AA_seqs", pattern = "_AA\\.fasta$", full.names = TRUE)
stopifnot(length(seq_files) > 0)

# Helper functions for sequence properties
calc_hydrophobicity <- function(seq_string) {
  # Kyte-Doolittle scale
  hydro_scale <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, 
                   G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, 
                   P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
  aa_freq <- table(strsplit(seq_string, "")[[1]])
  mean(hydro_scale[names(aa_freq)] * as.numeric(aa_freq), na.rm=TRUE)
}

calc_charge <- function(seq_string) {
  charge_scale <- c(R=1, K=1, H=0.5, D=-1, E=-1)  # Simplified at pH 7
  aa_freq <- table(strsplit(seq_string, "")[[1]])
  sum(charge_scale[names(aa_freq)] * as.numeric(aa_freq), na.rm=TRUE)
}

calc_gc_rich_aa <- function(seq_string) {
  gc_rich_aa <- c("G", "C", "A", "P", "R")  # AA from GC-rich codons
  aa_chars <- strsplit(seq_string, "")[[1]]
  sum(aa_chars %in% gc_rich_aa) / length(aa_chars)
}

# Individual amino acid proportions
calc_aa_proportions <- function(seq_string) {
  aa_chars <- strsplit(seq_string, "")[[1]]
  total_len <- length(aa_chars)
  
  # All 20 standard amino acids
  aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  proportions <- sapply(aa_list, function(aa) sum(aa_chars == aa) / total_len)
  names(proportions) <- paste0("prop_", aa_list)
  
  return(as.list(proportions))
}

calc_helix_propensity <- function(seq_string) {
  # Chou-Fasman helix propensities (simplified)
  helix_scale <- c(A=1.42, R=0.98, N=0.67, D=1.01, C=0.7, Q=1.11, E=1.51, 
                   G=0.57, H=1.0, I=1.08, L=1.21, K=1.16, M=1.45, F=1.13, 
                   P=0.57, S=0.77, T=0.83, W=1.08, Y=0.69, V=1.06)
  aa_freq <- table(strsplit(seq_string, "")[[1]])
  weighted.mean(helix_scale[names(aa_freq)], as.numeric(aa_freq), na.rm=TRUE)
}

# Transmembrane domain prediction (simplified sliding window)
predict_tm_domains <- function(seq_string, window_size = 20, hydro_threshold = 1.6) {
  aa_chars <- strsplit(seq_string, "")[[1]]
  seq_len <- length(aa_chars)
  
  if (seq_len < window_size) return(0)
  
  # Kyte-Doolittle scale
  hydro_scale <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, 
                   G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, 
                   P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
  
  tm_count <- 0
  in_tm <- FALSE
  
  for (i in 1:(seq_len - window_size + 1)) {
    window <- aa_chars[i:(i + window_size - 1)]
    window_hydro <- mean(hydro_scale[window], na.rm = TRUE)
    
    if (window_hydro > hydro_threshold && !in_tm) {
      tm_count <- tm_count + 1
      in_tm <- TRUE
    } else if (window_hydro <= hydro_threshold) {
      in_tm <- FALSE
    }
  }
  
  return(tm_count)
}

# Amino acid composition functions
calc_aa_composition <- function(seq_string) {
  aa_chars <- strsplit(seq_string, "")[[1]]
  total_len <- length(aa_chars)
  
  # Functional groups
  hydrophobic <- sum(aa_chars %in% c("A", "V", "I", "L", "M", "F", "W", "P")) / total_len
  polar <- sum(aa_chars %in% c("S", "T", "Y", "N", "Q", "C")) / total_len
  charged <- sum(aa_chars %in% c("R", "K", "D", "E", "H")) / total_len
  tiny <- sum(aa_chars %in% c("A", "G", "C", "S")) / total_len
  small <- sum(aa_chars %in% c("A", "G", "C", "S", "T", "P", "N", "D", "V")) / total_len
  aromatic <- sum(aa_chars %in% c("F", "W", "Y")) / total_len
  aliphatic <- sum(aa_chars %in% c("I", "L", "V")) / total_len
  sulfur <- sum(aa_chars %in% c("C", "M")) / total_len
  branched <- sum(aa_chars %in% c("I", "L", "V", "T")) / total_len
  
  return(list(
    hydrophobic = hydrophobic,
    polar = polar, 
    charged = charged,
    tiny = tiny,
    small = small,
    aromatic = aromatic,
    aliphatic = aliphatic,
    sulfur = sulfur,
    branched = branched
  ))
}

# Read sequences and calculate properties
seq_data <- data.table()
for(file in seq_files) {
  seqs <- readAAStringSet(file)
  
  # Extract sample IDs and gene names from headers
  headers <- names(seqs)
  sample_ids <- gsub("\\|.*", "", headers)
  gene_names <- gsub(".*Gene_([^|]+)\\|.*", "\\1", headers)
  seq_strings <- as.character(seqs)
  
  # Calculate all properties including AA composition and individual AA proportions
  aa_comp_list <- lapply(seq_strings, calc_aa_composition)
  aa_comp_df <- do.call(rbind, lapply(aa_comp_list, as.data.frame))
  
  aa_prop_list <- lapply(seq_strings, calc_aa_proportions)
  aa_prop_df <- do.call(rbind, lapply(aa_prop_list, as.data.frame))
  
  file_data <- data.table(
    sample = sample_ids,
    gene = gene_names,
    seq_length = width(seqs),
    hydrophobicity = sapply(seq_strings, function(x) hydrophobicity(x, scale = "KyteDoolittle")),
    net_charge = sapply(seq_strings, function(x) charge(x, pH = 7, pKscale = "EMBOSS")),
    gc_rich_aa = sapply(seq_strings, calc_gc_rich_aa),
    helix_propensity = sapply(seq_strings, function(x) calc_helix_propensity(x)),
    aromatic_content = sapply(seq_strings, function(x) sum(strsplit(x,"")[[1]] %in% c("F","W","Y"))/nchar(x)),
    tm_domains = sapply(seq_strings, predict_tm_domains),
    molecular_weight = sapply(seq_strings, mw),
    isoelectric_point = sapply(seq_strings, pI),
    aliphatic_index = sapply(seq_strings, aIndex),
    instability_index = sapply(seq_strings, instaIndex)
  )
  
  # Add AA composition and proportion columns
  file_data <- cbind(file_data, aa_comp_df, aa_prop_df)
  
  seq_data <- rbind(seq_data, file_data)
}
saveRDS(seq_data, "data/tmp/bigSeqData.rds")

stopifnot(nrow(seq_data) > 0)

# Merge with PCA data
# Assuming clean_embeds has sample identifiers that match seq_data$sample
pc1_data <- data.table(
  sample = clean_embeds$Sample,  # Adjust column name as needed
  gene = clean_embeds$Gene,
  PC1 = pca$x[,1]
)

merged_data <- merge(seq_data, pc1_data, by = c("sample", "gene"), all.x = FALSE, all.y = FALSE)
stopifnot(nrow(merged_data) > 0)

# Test all properties including new ones
properties <- c("seq_length", "hydrophobicity", "net_charge", "gc_rich_aa", "helix_propensity", 
                "aromatic_content", "tm_domains", "molecular_weight", "isoelectric_point", 
                "aliphatic_index", "instability_index", "hydrophobic", "polar", "charged", 
                "tiny", "small", "aromatic", "aliphatic", "sulfur", "branched",
                paste0("prop_", c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                                  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")))

cat("Correlations with PC1:\n")
cat("=====================\n")
all_cors <- data.table()
for(prop in properties) {
  cor_val <- cor(merged_data$PC1, merged_data[[prop]], use = "complete.obs")
  cor_test <- cor.test(merged_data$PC1, merged_data[[prop]])
  cat(sprintf("%-20s r = %6.3f, p = %.2e\n", prop, cor_val, cor_test$p.value))
  all_cors <- rbind(all_cors, data.table(property = prop, correlation = cor_val, abs_cor = abs(cor_val)))
}

# Find top correlations
top_cors <- all_cors[order(-abs_cor)][1:5]
cat("\nTop 5 correlations:\n")
print(top_cors)

# Gene analysis for transmembrane domains
cat("\nTransmembrane domain analysis by gene:\n")
tm_analysis <- merged_data[, .(
  n = .N,
  mean_pc1 = mean(PC1),
  mean_tm = mean(tm_domains),
  max_tm = max(tm_domains),
  prop_with_tm = mean(tm_domains > 0)
), by = gene][n >= 100][order(mean_pc1)]

print(tm_analysis)

# Create comprehensive plots
png("plots/pc1_comprehensive.png", width = 1600, height = 1200)
par(mfrow = c(3, 4))

# Plot top correlates
for(prop in head(all_cors[order(-abs_cor)]$property, 12)) {
  cor_val <- all_cors[property == prop]$correlation
  plot(merged_data[[prop]], merged_data$PC1,
       pch = 16, cex = 0.4, col = rainbow(length(unique(merged_data$gene)))[as.numeric(as.factor(merged_data$gene))],
       xlab = prop, ylab = "PC1",
       main = paste0(prop, "\nr = ", round(cor_val, 3)))
  if(abs(cor_val) > 0.1) {
    abline(lm(merged_data$PC1 ~ merged_data[[prop]]), col = "red", lwd = 2)
  }
}
dev.off()

# Special TM analysis
cat("\nPSB genes with TM domains:\n")
psb_tm <- merged_data[grepl("^psb", gene), .(
  gene = gene[1],
  mean_tm = mean(tm_domains),
  mean_pc1 = mean(PC1)
), by = gene][order(mean_pc1)]
print(psb_tm)

cat("\nNon-PSB genes with high TM content:\n")
non_psb_tm <- merged_data[!grepl("^psb", gene) & tm_domains > 2, .(
  gene = gene[1],
  mean_tm = mean(tm_domains),
  mean_pc1 = mean(PC1)
), by = gene][order(-mean_tm)]
print(head(non_psb_tm, 10))