library(Biostrings)
library(data.table)
library(Peptides)

# SINGLE MSA: RESIDUE-RESIDUE COVARIANCE

calc_residue_covariance <- function(aln_file, min_samples = 100) {
  aln <- readAAStringSet(aln_file)
  stopifnot(length(aln) > 0)
  
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  
  # Encode residues numerically (simple: convert to integer codes)
  # Treat gaps as NA
  aln_num <- apply(aln_mat, 2, function(col) {
    col[col == "-"] <- NA
    as.numeric(factor(col))
  })
  
  stopifnot(nrow(aln_num) >= min_samples)
  
  # Covariance matrix (pairwise complete obs)
  cov_mat <- cov(aln_num, use = "pairwise.complete.obs")
  
  list(cov = cov_mat, ids = rownames(aln_mat), n = nrow(aln_mat))
}

plot_residue_covariance <- function(cov_mat, title = "Residue Covariance") {
  image(1:nrow(cov_mat), 1:ncol(cov_mat), cov_mat,
        xlab = "Position", ylab = "Position", main = title,
        col = hcl.colors(50, "RdBu", rev = TRUE))
}

# ---------------------------------------------------------------------
# TWO MSAs: MERGED CROSS-GENE COVARIANCE
# ---------------------------------------------------------------------

calc_cross_gene_covariance <- function(aln_file1, aln_file2, min_samples = 100) {
  
  read_and_encode <- function(f) {
    aln <- readAAStringSet(f)
    stopifnot(length(aln) > 0)
    names(aln) <- sub("\\|.*", "", names(aln))
    aln_mat <- as.matrix(aln)
    aln_num <- apply(aln_mat, 2, function(col) {
      col[col == "-"] <- NA
      as.numeric(factor(col))
    })
    rownames(aln_num) <- names(aln)
    aln_num
  }
  
  mat1 <- read_and_encode(aln_file1)
  mat2 <- read_and_encode(aln_file2)
  
  common_ids <- intersect(rownames(mat1), rownames(mat2))
  stopifnot(length(common_ids) >= min_samples)
  
  mat1 <- mat1[common_ids, , drop = FALSE]
  mat2 <- mat2[common_ids, , drop = FALSE]
  
  # Concatenate
  merged <- cbind(mat1, mat2)
  
  cov_mat <- cov(merged, use = "pairwise.complete.obs")
  
  list(
    cov = cov_mat,
    n_pos_gene1 = ncol(mat1),
    n_pos_gene2 = ncol(mat2),
    ids = common_ids,
    n = length(common_ids)
  )
}

plot_cross_gene_covariance <- function(res, gene1_name = "Gene1", gene2_name = "Gene2") {
  cov_mat <- res$cov
  n1 <- res$n_pos_gene1
  n2 <- res$n_pos_gene2
  
  image(1:nrow(cov_mat), 1:ncol(cov_mat), cov_mat,
        xlab = "Position", ylab = "Position",
        main = paste("Cross-gene covariance:", gene1_name, "&", gene2_name),
        col = hcl.colors(50, "RdBu", rev = TRUE))
  
  # Mark gene boundaries
  abline(v = n1 + 0.5, col = "black", lwd = 2)
  abline(h = n1 + 0.5, col = "black", lwd = 2)
}

# ---------------------------------------------------------------------
# EXAMPLE USAGE
# ---------------------------------------------------------------------

# Single gene
res <- calc_residue_covariance("data/tmp/alignedGenes/rbcL_AA_aligned.fasta")
plot_residue_covariance(res$cov, title = "GENE1")
hist(res$cov)
summary(as.vector(res$cov))
440*440
# Two genes
# res2 <- calc_cross_gene_covariance(
#   "data/tmp/alignedGenes/GENE1_AA_aligned.fasta",
#   "data/tmp/alignedGenes/GENE2_AA_aligned.fasta"
# )
# plot_cross_gene_covariance(res2, "GENE1", "GENE2")
# ---- kidera factor ----
library(Biostrings)
library(Peptides)  # for AAdata$kideraFactors

# KIDERA FACTOR ENCODING

get_kidera_matrix <- function() {
  kf <- AAdata$kideraFactors
  # Convert to matrix: 20 AA x 10 factors
  aa <- names(kf[[1]])
  mat <- sapply(kf, function(x) x[aa])
  rownames(mat) <- aa
  mat
}

encode_kidera <- function(aln_mat) {
  kidera <- get_kidera_matrix()  # 20 AA x 10 factors
  valid_aa <- rownames(kidera)
  
  n_samples <- nrow(aln_mat)
  n_pos <- ncol(aln_mat)
  
  out <- matrix(NA_real_, nrow = n_samples, ncol = n_pos * 10)
  
  for (p in seq_len(n_pos)) {
    residues <- as.character(aln_mat[, p])
    valid_idx <- residues %in% valid_aa
    
    for (k in 1:10) {
      col_idx <- (p - 1) * 10 + k
      out[valid_idx, col_idx] <- kidera[residues[valid_idx], k]
    }
  }
  
  colnames(out) <- paste0("P", rep(seq_len(n_pos), each = 10), "_KF", 1:10)
  rownames(out) <- rownames(aln_mat)
  out
}

# ---------------------------------------------------------------------
# SINGLE MSA: KIDERA COVARIANCE
# ---------------------------------------------------------------------

calc_kidera_covariance <- function(aln_file, min_samples = 100,
                                   max_gap_frac = 0.5) {
  aln <- readAAStringSet(aln_file)
  stopifnot(length(aln) > 0)
  
  names(aln) <- sub("\\|.*", "", names(aln))
  aln_mat <- as.matrix(aln)
  
  # Filter high-gap positions
  gap_frac <- apply(aln_mat, 2, function(x) mean(x == "-"))
  keep <- gap_frac <= max_gap_frac
  stopifnot(sum(keep) > 1)
  aln_mat <- aln_mat[, keep, drop = FALSE]
  
  aln_mat[aln_mat == "-"] <- NA
  stopifnot(nrow(aln_mat) >= min_samples)
  
  kidera_mat <- encode_kidera(aln_mat)
  cov_mat <- cov(kidera_mat, use = "pairwise.complete.obs")
  
  list(cov = cov_mat, ids = rownames(aln_mat), n = nrow(aln_mat),
       n_pos = sum(keep), kept_positions = which(keep))
}

# ---------------------------------------------------------------------
# POSITION-LEVEL COVARIANCE (aggregate across Kidera factors)
# ---------------------------------------------------------------------

aggregate_to_position_cov <- function(cov_mat, n_pos) {
  # Collapse 10 Kidera factors per position to single position-position covariance
  # Using Frobenius norm of 10x10 block
  
  pos_cov <- matrix(NA_real_, n_pos, n_pos)
  
  for (i in seq_len(n_pos)) {
    for (j in seq_len(n_pos)) {
      idx_i <- ((i - 1) * 10 + 1):(i * 10)
      idx_j <- ((j - 1) * 10 + 1):(j * 10)
      block <- cov_mat[idx_i, idx_j]
      pos_cov[i, j] <- sqrt(sum(block^2, na.rm = TRUE))  # Frobenius norm
    }
  }
  pos_cov
}

plot_kidera_covariance <- function(res, title = "Residue Covariance (Kidera)") {
  pos_cov <- aggregate_to_position_cov(res$cov, res$n_pos)
  image(1:nrow(pos_cov), 1:ncol(pos_cov), pos_cov,
        xlab = "Position", ylab = "Position", main = title,
        col = hcl.colors(50, "YlOrRd"))
}

# ---------------------------------------------------------------------
# TWO MSAs: CROSS-GENE KIDERA COVARIANCE
# ---------------------------------------------------------------------

calc_cross_gene_kidera_covariance <- function(aln_file1, aln_file2, 
                                              min_samples = 100, max_gap_frac = 0.5) {
  
  read_and_filter <- function(f) {
    aln <- readAAStringSet(f)
    stopifnot(length(aln) > 0)
    names(aln) <- sub("\\|.*", "", names(aln))
    aln_mat <- as.matrix(aln)
    gap_frac <- apply(aln_mat, 2, function(x) mean(x == "-"))
    keep <- gap_frac <= max_gap_frac
    list(mat = aln_mat[, keep, drop = FALSE], n_pos = sum(keep))
  }
  
  r1 <- read_and_filter(aln_file1)
  r2 <- read_and_filter(aln_file2)
  
  common_ids <- intersect(rownames(r1$mat), rownames(r2$mat))
  stopifnot(length(common_ids) >= min_samples)
  
  mat1 <- r1$mat[common_ids, , drop = FALSE]
  mat2 <- r2$mat[common_ids, , drop = FALSE]
  
  mat1[mat1 == "-"] <- NA
  mat2[mat2 == "-"] <- NA
  
  kidera1 <- encode_kidera(mat1)
  kidera2 <- encode_kidera(mat2)
  merged <- cbind(kidera1, kidera2)
  
  cov_mat <- cov(merged, use = "pairwise.complete.obs")
  
  list(cov = cov_mat, n_pos_gene1 = r1$n_pos, n_pos_gene2 = r2$n_pos,
       ids = common_ids, n = length(common_ids))
}

plot_cross_gene_kidera <- function(res, gene1_name = "Gene1", gene2_name = "Gene2") {
  n1 <- res$n_pos_gene1
  n2 <- res$n_pos_gene2
  pos_cov <- aggregate_to_position_cov(res$cov, n1 + n2)
  
  image(1:nrow(pos_cov), 1:ncol(pos_cov), pos_cov,
        xlab = "Position", ylab = "Position",
        main = paste("Cross-gene covariance:", gene1_name, "&", gene2_name),
        col = hcl.colors(50, "YlOrRd"))
  
  abline(v = n1 + 0.5, col = "black", lwd = 2)
  abline(h = n1 + 0.5, col = "black", lwd = 2)
}


# ---------------------------------------------------------------------
# SINGLE GENE
# ---------------------------------------------------------------------

aln_file <- "data/tmp/alignedGenes/rbcL_AA_aligned.fasta"

res <- calc_kidera_covariance(aln_file)
res$n       # number of samples
res$n_pos   # number of positions retained

plot_kidera_covariance(res, title = "rbcL Kidera Covariance")

# ---------------------------------------------------------------------
# TWO GENES
# ---------------------------------------------------------------------

res2 <- calc_cross_gene_kidera_covariance(
  "data/tmp/alignedGenes/GENE1_AA_aligned.fasta",
  "data/tmp/alignedGenes/GENE2_AA_aligned.fasta"
)

res2$n            # shared samples
res2$n_pos_gene1  # positions in gene 1
res2$n_pos_gene2  # positions in gene 2

plot_cross_gene_kidera(res2, "GENE1", "GENE2")

# ---------------------------------------------------------------------
# LOOP OVER ALL GENE PAIRS
# ---------------------------------------------------------------------

aln_files <- list.files("data/tmp/alignedGenes/", 
                        pattern = "_AA_aligned\\.fasta$", full.names = TRUE)
genes <- sub("_AA_aligned\\.fasta", "", basename(aln_files))

# All pairwise combinations
pairs <- combn(genes, 2, simplify = FALSE)

pdf("results/cross_gene_kidera_covariance.pdf", width = 8, height = 8)
for (p in pairs) {
  f1 <- file.path("data/tmp/alignedGenes", paste0(p[1], "_AA_aligned.fasta"))
  f2 <- file.path("data/tmp/alignedGenes", paste0(p[2], "_AA_aligned.fasta"))
  
  res <- tryCatch(
    calc_cross_gene_kidera_covariance(f1, f2),
    error = function(e) NULL
  )
  
  if (!is.null(res)) {
    plot_cross_gene_kidera(res, p[1], p[2])
    message("Plotted: ", p[1], " x ", p[2], " (n=", res$n, ")")
  }
}
dev.off()
