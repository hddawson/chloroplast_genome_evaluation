# ---- plastome mapping /circos plot ---- 
library(data.table)

# Load mapped results
gwas <- fread("results/gwas_arabidopsis_mapped.csv")
stopifnot(nrow(gwas) > 0)

# Basic info
message("Positions: ", nrow(gwas))
message("Genes: ", uniqueN(gwas$Gene))
message("GWAS hits: ", sum(gwas$gwas_hit))
message("Genomic range: ", min(gwas$Genomic_Position), " - ", max(gwas$Genomic_Position))

# ---- CIRCULAR PLOT PREP ----
# Arabidopsis chloroplast is ~154 kb
GENOME_SIZE <- 153153  # adjust to actual size from GenBank

# Convert to radians for circular plotting
gwas[, theta := 2 * pi * Genomic_Position / GENOME_SIZE]

# ---- LINEAR MANHATTAN ----
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))

# Color by gene (alternating)
genes <- unique(gwas$Gene)
gene_cols <- setNames(rep(c("steelblue", "grey50"), length.out = length(genes)), genes)

plot(gwas$Genomic_Position / 1000, gwas$neg_log10_p,
     pch = 16, cex = 0.6,
     col = ifelse(gwas$gwas_hit, "red", gene_cols[gwas$Gene]),
     xlab = "Arabidopsis chloroplast position (kb)",
     ylab = "-log10(P)",
     main = "GWAS mapped to At chloroplast - linear")

# Gene labels
gene_mids <- gwas[, .(mid = median(Genomic_Position) / 1000), by = Gene]

axis(1, at = gene_mids$mid, labels = gene_mids$Gene, las = 2, cex.axis = 0.7)


# ---- CIRCULAR PLOT ----
par(mar = c(1, 1, 2, 1))
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), 
     asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "AA GWAS hits mapped to At chloroplast")

# Draw chromosome circle
theta_seq <- seq(0, 2 * pi, length.out = 200)
lines(cos(theta_seq), sin(theta_seq), col = "grey70", lwd = 2)

# Plot points - radius scaled by -log10(p)
max_p <- max(gwas$neg_log10_p, na.rm = TRUE)
gwas[, radius := 1 + 0.4 * neg_log10_p / max_p]

points(gwas$radius * cos(gwas$theta), 
       gwas$radius * sin(gwas$theta),
       pch = 16, cex = 0.4,
       col = ifelse(gwas$gwas_hit, "red", "grey50"))

# Gene labels on inner circle
gene_theta <- gwas[, .(theta = median(theta)), by = Gene]
text(0.75 * cos(gene_theta$theta), 
     0.75 * sin(gene_theta$theta),
     gene_theta$Gene, cex = 0.5, srt = 0)

# Add scale
text(0, 0, paste0(round(GENOME_SIZE/1000), " kb"), cex = 0.8)

plot(gwas$Genomic_Position, -log10(gwas$P_aa_with_pcs))
