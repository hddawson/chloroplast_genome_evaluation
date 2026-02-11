# ---- RUBISCO 3D STRUCTURE ANALYSIS ----
# Step 1: Get PDB structure and map alignment positions

library(bio3d)
library(data.table)

download.file(
  "https://files.rcsb.org/download/1RCX.cif",
  "1RCX.cif",
  method = "libcurl"
)
pdb_1rcx <- read.pdb("1RCX.cif")

pdb_1rcx <- read.pdb("1RCX", type = "cif")

options(download.file.method = "libcurl")
library(bio3d)
pdb_1rcx <- read.pdb("1RCX")


# Download Arabidopsis Rubisco (1.5 Å resolution)
pdb <- read.pdb("5IU0")

# Examine structure
print(pdb)

# Chain A = RbcL (large subunit), Chain B = RbcS (small subunit)
# The asymmetric unit has A2B2, full complex is L8S8
table(pdb$atom$chain)

# Extract just RbcL (chain A)
rbcL_pdb <- trim.pdb(pdb, chain = "A")

# Get CA atoms for spatial analysis
ca_atoms <- atom.select(rbcL_pdb, elety = "CA")
rbcL_ca <- trim.pdb(rbcL_pdb, ca_atoms)

# Extract coordinates and sequence
coords <- rbcL_ca$atom[, c("resno", "resid", "x", "y", "z")]
coords <- as.data.table(coords)
setnames(coords, c("pdb_resno", "aa", "x", "y", "z"))

summary(coords)
pairs(coords[, c("x","y","z")])

# Check
stopifnot(nrow(coords) > 400)  # rbcL is ~475 residues
message("Extracted ", nrow(coords), " CA atoms from rbcL")
print(head(coords, 10))
print(tail(coords, 10))

# Get the sequence
pdb_seq <- paste(coords$aa, collapse = "")
message("PDB sequence length: ", nchar(pdb_seq))


# ---- Step 2: Get rbcL GWAS sites and alignment ----
model_files <- list.files("results/residue_models_triple/", 
                          pattern = "_effects\\.rds$", full.names = TRUE)
stopifnot(length(model_files) > 0)

gwas_list <- lapply(model_files, function(f) {
  cat(f)
  models <- readRDS(f)
  rbindlist(lapply(models, function(m) {
    data.table(
      Gene = m$Gene,
      Position = m$Aligned_Position,
      P_aa_only = m$P_aa_only,
      P_aa_with_pcs = m$P_aa_with_pcs,
      N = m$N
    )
  }))
})
sites_df <- rbindlist(gwas_list)
# Filter to rbcL
rbcL_sites <- sites_df[Gene == "rbcL"]
stopifnot(nrow(rbcL_sites) > 0)
message("rbcL GWAS sites: ", nrow(rbcL_sites), " positions")
message("Position range: ", min(rbcL_sites$Position), " - ", max(rbcL_sites$Position))

# Check p-value distribution
message("Significant sites (P_aa_with_pcs < 0.05): ", 
        sum(rbcL_sites$P_aa_with_pcs < 0.05, na.rm = TRUE))

arabidopsis_th_id <- data$ID[grep("thaliana", data$Organism)]

# ---- Step 3: Map alignment positions to PDB positions ----

library(Biostrings)

# Read alignment
rbcL_aln <- readAAStringSet("data/tmp/alignedGenes/rbcL_AA_aligned.fasta")

# Get Arabidopsis sequence from alignment
at_idx <- grep(arabidopsis_th_id, names(rbcL_aln))
stopifnot(length(at_idx) == 1)
at_aln_seq <- as.character(rbcL_aln[[at_idx]])

# Convert aligned Arabidopsis to ungapped sequence
at_ungapped <- gsub("-", "", at_aln_seq)
message("Arabidopsis ungapped length: ", nchar(at_ungapped))

# PDB sequence (from coords)
pdb_seq <- paste(coords$aa, collapse = "")
# Convert 3-letter to 1-letter
aa_map <- c(ALA="A", CYS="C", ASP="D", GLU="E", PHE="F", GLY="G", HIS="H", 
            ILE="I", LYS="K", LEU="L", MET="M", ASN="N", PRO="P", GLN="Q", 
            ARG="R", SER="S", THR="T", VAL="V", TRP="W", TYR="Y")
pdb_seq_1letter <- paste(aa_map[coords$aa], collapse = "")
message("PDB sequence length: ", nchar(pdb_seq_1letter))

# PDB starts at residue 13, so the full protein sequence would have 12 missing N-terminal residues
# Align Arabidopsis ungapped to PDB sequence
# The PDB sequence should be a near-perfect substring of the Arabidopsis sequence

# Find where PDB sequence matches in Arabidopsis
pdb_in_at <- regexpr(pdb_seq_1letter, at_ungapped, fixed = TRUE)
if (pdb_in_at == -1) {
  # Try allowing mismatches with pairwise alignment
  message("Exact match not found, doing pairwise alignment...")
  pw_aln <- pairwiseAlignment(pdb_seq_1letter, at_ungapped, type = "local")
  print(pw_aln)
} else {
  message("PDB sequence found at position ", pdb_in_at, " in Arabidopsis sequence")
}

# Create mapping: alignment position -> PDB residue number
# Step through Arabidopsis aligned sequence
aln_chars <- strsplit(at_aln_seq, "")[[1]]
ungapped_pos <- 0L
aln_to_pdb <- data.table(
  aln_pos = seq_along(aln_chars),
  aln_char = aln_chars,
  at_ungapped_pos = NA_integer_,
  pdb_resno = NA_integer_
)

for (i in seq_along(aln_chars)) {
  if (aln_chars[i] != "-") {
    ungapped_pos <- ungapped_pos + 1L
    aln_to_pdb$at_ungapped_pos[i] <- ungapped_pos
    
    # PDB residues 13-475 correspond to ungapped positions 13-475 of full protein
    # (assuming Arabidopsis protein starts at Met and PDB numbering matches protein numbering)
    if (ungapped_pos >= min(coords$pdb_resno) & ungapped_pos <= max(coords$pdb_resno)) {
      aln_to_pdb$pdb_resno[i] <- ungapped_pos
    }
  }
}

# Check mapping
message("Alignment positions with PDB mapping: ", sum(!is.na(aln_to_pdb$pdb_resno)))
print(aln_to_pdb[!is.na(pdb_resno)][1:10])
print(aln_to_pdb[!is.na(pdb_resno)][(nrow(aln_to_pdb[!is.na(pdb_resno)])-9):nrow(aln_to_pdb[!is.na(pdb_resno)])])

# ---- Step 4: Merge GWAS sites with PDB coordinates ----

# Add PDB coordinates to the mapping
aln_to_pdb <- merge(aln_to_pdb, coords[, .(pdb_resno, x, y, z)], 
                    by = "pdb_resno", all.x = TRUE)

# Merge with rbcL GWAS sites
rbcL_sites <- sites_df[Gene == "rbcL"]
rbcL_sites <- merge(rbcL_sites, aln_to_pdb[, .(aln_pos, pdb_resno, x, y, z)], 
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

# Check coverage
message("rbcL sites with PDB coordinates: ", sum(!is.na(rbcL_sites$pdb_resno)), 
        " of ", nrow(rbcL_sites))

# Get genome-wide thresholds from full sites_df
thresh_aa_only <- quantile(sites_df$P_aa_only, 0.25, na.rm = TRUE)
thresh_aa_pcs <- quantile(sites_df$P_aa_with_pcs, 0.05, na.rm = TRUE)
message("Genome-wide thresholds: P_aa_only < ", signif(thresh_aa_only, 3), 
        ", P_aa_with_pcs < ", signif(thresh_aa_pcs, 3))

# Classify sites using genome-wide thresholds
rbcL_sites[, site_class := "not_sig"]
rbcL_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
rbcL_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
rbcL_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(rbcL_sites$site_class))

# Filter to sites with 3D coordinates
rbcL_3d <- rbcL_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(rbcL_3d$site_class))

# ---- Step 5: Spatial autocorrelation of significant sites ----

# Calculate pairwise distance matrix for all sites with coords
dist_mat <- as.matrix(dist(rbcL_3d[, .(x, y, z)]))

hist(dist_mat)

# Binary indicator: is site significant?
rbcL_3d[, is_sig := site_class == "sig_both"]

message("Significant sites: ", sum(rbcL_3d$is_sig), " of ", nrow(rbcL_3d))

# Moran's I for spatial autocorrelation
# Using inverse distance weights (with cutoff to avoid tiny weights)
calc_morans_i <- function(values, dist_mat, max_dist = 20) {
  n <- length(values)
  stopifnot(n == nrow(dist_mat))
  
  # Inverse distance weights, zero beyond max_dist
  W <- 1 / dist_mat
  W[dist_mat > max_dist | dist_mat == 0] <- 0
  
  # Row-standardize
  W <- W / rowSums(W)
  W[is.nan(W)] <- 0
  
  # Moran's I
  x <- values - mean(values)
  I <- (n / sum(W)) * (sum(W * outer(x, x)) / sum(x^2))
  
  # Expected value under null
  E_I <- -1 / (n - 1)
  
  return(list(I = I, E_I = E_I))
}

# Observed Moran's I
obs_moran <- calc_morans_i(as.numeric(rbcL_3d$is_sig), dist_mat)
message("Observed Moran's I: ", round(obs_moran$I, 4), 
        " (expected under null: ", round(obs_moran$E_I, 4), ")")

# Permutation test
set.seed(42)
n_perm <- 999
perm_I <- replicate(n_perm, {
  perm_sig <- sample(rbcL_3d$is_sig)
  calc_morans_i(as.numeric(perm_sig), dist_mat)$I
})

p_value <- (sum(perm_I >= obs_moran$I) + 1) / (n_perm + 1)
message("Permutation p-value (one-sided, clustering): ", p_value)

# Plot permutation distribution
hist(perm_I, breaks = 30, main = "Moran's I Permutation Distribution",
     xlab = "Moran's I", col = "gray80")
abline(v = obs_moran$I, col = "red", lwd = 2)
abline(v = obs_moran$E_I, col = "blue", lwd = 2, lty = 2)
legend("topright", c("Observed", "Expected (null)"), col = c("red", "blue"), lwd = 2, lty = c(1, 2))
summary(perm_I)

# ---- Step 6: Distance-based analysis and functional site proximity ----

# Test multiple distance cutoffs
rbcL_3d[, is_sig := site_class == "sig_both"]  # strictest class

cutoffs <- c(8, 10, 15, 20, 30, 40)
moran_by_cutoff <- sapply(cutoffs, function(d) {
  calc_morans_i(as.numeric(rbcL_3d$is_sig), dist_mat, max_dist = d)$I
})
names(moran_by_cutoff) <- cutoffs
print(round(moran_by_cutoff, 4))

# Alternative: Ripley's K / nearest neighbor approach
# Are sig sites closer to each other than expected?

sig_idx <- which(rbcL_3d$is_sig)
nonsig_idx <- which(!rbcL_3d$is_sig)

# Mean nearest-neighbor distance among sig sites
if (length(sig_idx) > 1) {
  sig_dist <- dist_mat[sig_idx, sig_idx]
  diag(sig_dist) <- NA
  mean_nn_sig <- mean(apply(sig_dist, 1, min, na.rm = TRUE))
  
  # Permutation test: random sets of same size
  set.seed(123)
  n_perm <- 999
  perm_nn <- replicate(n_perm, {
    rand_idx <- sample(nrow(rbcL_3d), length(sig_idx))
    rand_dist <- dist_mat[rand_idx, rand_idx]
    diag(rand_dist) <- NA
    mean(apply(rand_dist, 1, min, na.rm = TRUE))
  })
  
  p_nn <- (sum(perm_nn <= mean_nn_sig) + 1) / (n_perm + 1)
  
  message("\nNearest-neighbor analysis (sig_both sites):")
  message("Mean NN distance among sig sites: ", round(mean_nn_sig, 2), " Å")
  message("Expected (permutation mean): ", round(mean(perm_nn), 2), " Å")
  message("P-value (clustering): ", p_nn)
  
  hist(perm_nn, breaks = 30, main = "Nearest Neighbor Distance (sig_both)",
       xlab = "Mean NN distance (Å)", col = "gray80")
  abline(v = mean_nn_sig, col = "red", lwd = 2)
}

# ---- Step 7: Build L8 multimer from biological assembly ----
# Re-read PDB with biological assembly info
pdb <- read.pdb("5IU0")

# Check what chains we have in asymmetric unit
table(pdb$atom$chain)

# 5IU0 asymmetric unit is A2B2 (2 RbcL + 2 RbcS)
# Full L8S8 is generated by crystallographic symmetry
# We need to fetch the biological assembly

# Download biological assembly (symmetry-expanded)
# ---- Step 7: Build L8 multimer ----

# Download biological assembly PDB file directly
bio_url <- "https://files.rcsb.org/download/5IU0-assembly1.cif"
download.file(bio_url, "5IU0_bio.cif", mode = "wb")

# Try reading as mmCIF
pdb_bio <- read.cif("5IU0_bio.cif")

message("Chains in biological assembly:")
print(table(pdb_bio$atom$chain))
message("Total atoms: ", nrow(pdb_bio$atom))
# Check chains in biological assembly
message("Chains in biological assembly:")
print(table(pdb_bio$atom$chain))

# If biounit doesn't work, we can build it manually using BIOMT records
# Let's check what we got
message("Total atoms: ", nrow(pdb_bio$atom))

cat(pdb$remark$biomat, sep = "\n")
pdb$remark$biomat

# ---- Try 1RCX - spinach Rubisco L8S8 ----

pdb_1rcx <- read.pdb("1RCX")
print(pdb_1rcx)

# Check chains
message("Chains in 1RCX:")
print(table(pdb_1rcx$atom$chain))

# Check biomat
pdb_1rcx$remark$biomat


# 
# ---- Extract L8 from 1RCX ----

# Identify large subunit chains (RbcL) by atom count
chain_counts <- table(pdb_1rcx$atom$chain)
rbcL_chains <- names(chain_counts[chain_counts > 3000])
rbcS_chains <- names(chain_counts[chain_counts < 2000])

message("RbcL chains: ", paste(rbcL_chains, collapse = ", "))
message("RbcS chains: ", paste(rbcS_chains, collapse = ", "))

# Extract CA atoms from first RbcL chain (L) to get reference sequence
chain_L <- trim.pdb(pdb_1rcx, chain = "L")
ca_L <- atom.select(chain_L, elety = "CA")
chain_L_ca <- trim.pdb(chain_L, ca_L)

coords_1rcx <- as.data.table(chain_L_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_1rcx, c("pdb_resno", "aa", "x", "y", "z"))

message("Chain L residues: ", nrow(coords_1rcx))
message("Residue range: ", min(coords_1rcx$pdb_resno), "-", max(coords_1rcx$pdb_resno))

# Get 1-letter sequence
aa_map <- c(ALA="A", CYS="C", ASP="D", GLU="E", PHE="F", GLY="G", HIS="H", 
            ILE="I", LYS="K", LEU="L", MET="M", ASN="N", PRO="P", GLN="Q", 
            ARG="R", SER="S", THR="T", VAL="V", TRP="W", TYR="Y")
seq_1rcx <- paste(aa_map[coords_1rcx$aa], collapse = "")
message("1RCX sequence length: ", nchar(seq_1rcx))

# Align to Arabidopsis sequence
pw_aln <- pairwiseAlignment(seq_1rcx, at_ungapped, type = "global")
print(pw_aln)

# ---- EDA: Visualize Rubisco L8S8 structure ----

# 1. 3D scatter of all CA atoms, colored by chain type
library(ggplot2)

all_ca <- atom.select(pdb_1rcx, elety = "CA")
all_ca_pdb <- trim.pdb(pdb_1rcx, all_ca)
ca_df <- as.data.table(all_ca_pdb$atom[, c("chain", "resno", "x", "y", "z")])
ca_df[, subunit := ifelse(chain %in% rbcL_chains, "RbcL", "RbcS")]

# 2D projections
p1 <- ggplot(ca_df, aes(x = x, y = y, color = subunit)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("RbcL" = "steelblue", "RbcS" = "coral")) +
  coord_fixed() +
  labs(title = "Rubisco L8S8 - XY projection", x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p2 <- ggplot(ca_df, aes(x = x, y = z, color = chain)) +
  geom_point(alpha = 0.4, size = 0.5) +
  coord_fixed() +
  labs(title = "Rubisco L8S8 - XZ projection (colored by chain)", x = "X (Å)", y = "Z (Å)") +
  theme_classic() +
  theme(legend.position = "none")
p1
p2
# 2. Distance from center of mass by subunit type
ca_df[, dist_from_center := sqrt((x - mean(x))^2 + (y - mean(y))^2 + (z - mean(z))^2)]

p3 <- ggplot(ca_df, aes(x = dist_from_center, fill = subunit)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("RbcL" = "steelblue", "RbcS" = "coral")) +
  labs(title = "Distance from complex center", x = "Distance (Å)", y = "Count") +
  theme_classic()
p3
# 3. Per-residue distance from center (shows which parts are core vs peripheral)
chain_L_df <- ca_df[chain == "L"]
chain_L_df[, dist_from_center := sqrt((x - mean(ca_df$x))^2 + 
                                        (y - mean(ca_df$y))^2 + 
                                        (z - mean(ca_df$z))^2)]

p4 <- ggplot(chain_L_df, aes(x = resno, y = dist_from_center)) +
  geom_line(color = "steelblue") +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, color = "red") +
  labs(title = "RbcL: residue position vs distance from complex center",
       x = "Residue number", y = "Distance from center (Å)") +
  theme_classic()
p4
# 4. Contact map within one RbcL chain
dist_within_L <- as.matrix(dist(chain_L_df[, .(x, y, z)]))
rownames(dist_within_L) <- chain_L_df$resno
colnames(dist_within_L) <- chain_L_df$resno

# Plot as heatmap (contacts < 10Å)
contact_df <- as.data.table(which(dist_within_L < 10 & dist_within_L > 0, arr.ind = TRUE))
contact_df[, res1 := chain_L_df$resno[row]]
contact_df[, res2 := chain_L_df$resno[col]]

p5 <- ggplot(contact_df, aes(x = res1, y = res2)) +
  geom_point(size = 0.1, alpha = 0.5) +
  coord_fixed() +
  labs(title = "RbcL contact map (<10Å)", x = "Residue", y = "Residue") +
  theme_classic()
p5
# Display
library(patchwork)
(p1 + p2) / (p3 + p4)
p5

# ---- Biochemist-relevant plots ----

# 1. Active site location - RuBP binding site
# 1RCX has RUB (ribulose-1,5-bisphosphate analog) bound - find it
ligand_atoms <- pdb_1rcx$atom[pdb_1rcx$atom$resid == "RUB", ]
message("RUB ligands found: ", length(unique(ligand_atoms$chain)))

# Get ligand centers for each chain
ligand_centers <- as.data.table(ligand_atoms)[, .(
  lig_x = mean(x), lig_y = mean(y), lig_z = mean(z)
), by = chain]
print(ligand_centers)

# 2. Distance of each residue to nearest active site
chain_L_df <- ca_df[chain == "L"]

# Find closest RUB to chain L
chain_L_df[, dist_to_active := sqrt((x - ligand_centers[1, lig_x])^2 + 
                                      (y - ligand_centers[1, lig_y])^2 + 
                                      (z - ligand_centers[1, lig_z])^2)]

p_active <- ggplot(chain_L_df, aes(x = resno, y = dist_to_active)) +
  geom_line() +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  labs(title = "Distance to active site (RuBP)", 
       subtitle = "Red line = 10Å (typical interaction cutoff)",
       x = "Residue", y = "Distance (Å)") +
  theme_classic()
p_active
# 3. Inter-subunit contacts - where do L-L and L-S interfaces occur?
# Get CA coords for neighboring chain
chain_B <- ca_df[chain == "B"]  # adjacent RbcL
chain_S <- ca_df[chain == "S"]  # adjacent RbcS

# Distance from each L residue to closest B residue (L-L interface)
chain_L_df[, dist_to_L2 := sapply(1:.N, function(i) {
  min(sqrt((x[i] - chain_B$x)^2 + (y[i] - chain_B$y)^2 + (z[i] - chain_B$z)^2))
})]

# Distance to closest RbcS (L-S interface)
chain_L_df[, dist_to_S := sapply(1:.N, function(i) {
  min(sqrt((x[i] - chain_S$x)^2 + (y[i] - chain_S$y)^2 + (z[i] - chain_S$z)^2))
})]

p_interface <- ggplot(chain_L_df, aes(x = resno)) +
  geom_line(aes(y = dist_to_L2, color = "L-L interface")) +
  geom_line(aes(y = dist_to_S, color = "L-S interface")) +
  geom_hline(yintercept = 8, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("L-L interface" = "darkblue", "L-S interface" = "darkred")) +
  labs(title = "Subunit interface distances",
       subtitle = "Dashed = 8Å contact threshold",
       x = "Residue", y = "Distance to nearest neighbor subunit (Å)") +
  theme_classic() +
  theme(legend.position = "top")

p_interface
# 4. Classify residues by structural context
chain_L_df[, context := "bulk"]
chain_L_df[dist_to_active < 15, context := "near_active_site"]
chain_L_df[dist_to_L2 < 8, context := "L-L_interface"]
chain_L_df[dist_to_S < 8, context := "L-S_interface"]
chain_L_df[dist_to_L2 < 8 & dist_to_active < 15, context := "active_site_interface"]  # L-L interface forms active site!

message("\nResidue context breakdown:")
print(table(chain_L_df$context))

p_context <- ggplot(chain_L_df, aes(x = resno, y = 1, fill = context)) +
  geom_tile() +
  scale_fill_manual(values = c("bulk" = "gray80", 
                               "near_active_site" = "red",
                               "L-L_interface" = "blue",
                               "L-S_interface" = "purple",
                               "active_site_interface" = "darkred")) +
  labs(title = "RbcL structural context by residue", x = "Residue", y = "") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_context
# 5. Key: Active site is at L-L dimer interface!
# Show this with a 2D projection
p_complex <- ggplot() +
  geom_point(data = ca_df[subunit == "RbcL"], aes(x = x, y = y), 
             color = "gray70", size = 0.3, alpha = 0.5) +
  geom_point(data = ca_df[subunit == "RbcS"], aes(x = x, y = y), 
             color = "gray40", size = 0.3, alpha = 0.5) +
  geom_point(data = ligand_centers, aes(x = lig_x, y = lig_y), 
             color = "red", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "Rubisco L8S8 with active sites (red diamonds)",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()
p_complex
# Display
p_complex + p_active
p_interface / p_context

# ---- Validate structural context classification ----

# 1. 3D-ish view: show chain L colored by context, with ligand
p_validate_xy <- ggplot() +
  geom_point(data = chain_L_df, aes(x = x, y = y, color = context), size = 1.5) +
  geom_point(data = ligand_centers[chain == "L"], aes(x = lig_x, y = lig_y), 
             color = "black", size = 5, shape = 18) +
  geom_point(data = chain_B[, .(x, y)], aes(x = x, y = y), 
             color = "gray70", size = 0.5, alpha = 0.3) +
  geom_point(data = chain_S[, .(x, y)], aes(x = x, y = y), 
             color = "gray40", size = 0.5, alpha = 0.3) +
  scale_color_manual(values = c("bulk" = "gray80", 
                                "near_active_site" = "orange",
                                "L-L_interface" = "blue",
                                "L-S_interface" = "purple",
                                "active_site_interface" = "red")) +
  coord_fixed() +
  labs(title = "Chain L (colored) with neighbors (gray)", 
       subtitle = "Black diamond = RuBP ligand",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_validate_xz <- ggplot() +
  geom_point(data = chain_L_df, aes(x = x, y = z, color = context), size = 1.5) +
  geom_point(data = ligand_centers[chain == "L"], aes(x = lig_x, y = lig_z), 
             color = "black", size = 5, shape = 18) +
  scale_color_manual(values = c("bulk" = "gray80", 
                                "near_active_site" = "orange",
                                "L-L_interface" = "blue",
                                "L-S_interface" = "purple",
                                "active_site_interface" = "red")) +
  coord_fixed() +
  labs(title = "XZ projection", x = "X (Å)", y = "Z (Å)") +
  theme_classic()

# 2. Distance distributions by context - sanity check
p_dist_active <- ggplot(chain_L_df, aes(x = context, y = dist_to_active, fill = context)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("bulk" = "gray80", 
                               "near_active_site" = "orange",
                               "L-L_interface" = "blue",
                               "L-S_interface" = "purple",
                               "active_site_interface" = "red")) +
  labs(title = "Distance to active site by context", 
       subtitle = "Red dashed = 15Å threshold",
       y = "Distance to RuBP (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_dist_L2 <- ggplot(chain_L_df, aes(x = context, y = dist_to_L2, fill = context)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "blue") +
  scale_fill_manual(values = c("bulk" = "gray80", 
                               "near_active_site" = "orange",
                               "L-L_interface" = "blue",
                               "L-S_interface" = "purple",
                               "active_site_interface" = "red")) +
  labs(title = "Distance to adjacent RbcL by context",
       subtitle = "Blue dashed = 8Å threshold",
       y = "Distance to chain B (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_dist_S <- ggplot(chain_L_df, aes(x = context, y = dist_to_S, fill = context)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "purple") +
  scale_fill_manual(values = c("bulk" = "gray80", 
                               "near_active_site" = "orange",
                               "L-L_interface" = "blue",
                               "L-S_interface" = "purple",
                               "active_site_interface" = "red")) +
  labs(title = "Distance to RbcS by context",
       subtitle = "Purple dashed = 8Å threshold",
       y = "Distance to chain S (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# 3. Sequence position view - where are these regions in primary structure?
p_seq <- ggplot(chain_L_df, aes(x = resno, y = 0)) +
  geom_point(aes(color = context), size = 2, shape = 15) +
  scale_color_manual(values = c("bulk" = "gray80", 
                                "near_active_site" = "orange",
                                "L-L_interface" = "blue",
                                "L-S_interface" = "purple",
                                "active_site_interface" = "red")) +
  labs(title = "Structural context along RbcL sequence", x = "Residue number") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), legend.position = "top")

# Display
(p_validate_xy + p_validate_xz) / p_seq
(p_dist_active + p_dist_L2 + p_dist_S)

# ---- XY projection with context categories ----

p_xy_context <- ggplot() +
  # Background: other chains in light gray
  geom_point(data = ca_df[chain != "L"], aes(x = x, y = y), 
             color = "gray85", size = 0.3, alpha = 0.4) +
  # Chain L colored by context
  geom_point(data = chain_L_df, aes(x = x, y = y, color = context), size = 2) +
  # Active site ligands
  geom_point(data = ligand_centers, aes(x = lig_x, y = lig_y), 
             color = "black", size = 4, shape = 18) +
  scale_color_manual(values = c("bulk" = "gray50", 
                                "near_active_site" = "orange",
                                "L-L_interface" = "blue",
                                "L-S_interface" = "purple",
                                "active_site_interface" = "red"),
                     name = "Context") +
  coord_fixed() +
  labs(title = "RbcL structural context in L8S8 complex",
       subtitle = "Chain L highlighted; black diamonds = active sites (RuBP)",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic() +
  theme(legend.position = "right")

p_xy_context

# ---- Map GWAS sites to 1RCX and test context enrichment ----

# First, align 1RCX (spinach) to your alignment via Arabidopsis
# 1RCX residue numbering starts at 1, chain L

# Check alignment between spinach (1RCX) and Arabidopsis
pw_aln <- pairwiseAlignment(seq_1rcx, at_ungapped, type = "global")
print(pw_aln)
message("Alignment score: ", score(pw_aln))
message("Percent identity: ", round(pid(pw_aln), 1), "%")

# Create mapping: 1RCX resno -> Arabidopsis ungapped position
# Extract aligned sequences
aln_pattern <- as.character(pattern(pw_aln))  # 1RCX (with gaps)
aln_subject <- as.character(subject(pw_aln))  # Arabidopsis (with gaps)

stopifnot(nchar(aln_pattern) == nchar(aln_subject))

# Build position correspondence
rcx_pos <- 0L
at_pos <- 0L
rcx_to_at <- data.table(
  rcx_resno = integer(),
  at_ungapped = integer()
)

for (i in 1:nchar(aln_pattern)) {
  p_char <- substr(aln_pattern, i, i)
  s_char <- substr(aln_subject, i, i)
  
  if (p_char != "-") rcx_pos <- rcx_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    rcx_to_at <- rbind(rcx_to_at, data.table(rcx_resno = rcx_pos, at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(rcx_to_at), " positions between 1RCX and Arabidopsis")

# Now link to alignment positions via aln_to_pdb (which has at_ungapped_pos)
# aln_to_pdb has: aln_pos, at_ungapped_pos
at_to_aln <- aln_to_pdb[!is.na(at_ungapped_pos), .(aln_pos, at_ungapped_pos)]

# Merge: 1RCX -> Arabidopsis ungapped -> alignment position
rcx_to_aln <- merge(rcx_to_at, at_to_aln, 
                    by.x = "at_ungapped", by.y = "at_ungapped_pos")

message("Mapped ", nrow(rcx_to_aln), " 1RCX positions to alignment positions")

# Add structural context from chain_L_df
chain_L_df[, rcx_resno := resno]
rcx_to_aln <- merge(rcx_to_aln, chain_L_df[, .(rcx_resno, context, dist_to_active, dist_to_L2, dist_to_S)],
                    by = "rcx_resno", all.x = TRUE)

# Merge with GWAS sites
rbcL_context <- merge(rbcL_sites, rcx_to_aln[, .(aln_pos, context, dist_to_active, dist_to_L2, dist_to_S)],
                      by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\nGWAS sites with structural context: ", sum(!is.na(rbcL_context$context)), " of ", nrow(rbcL_context))
print(table(rbcL_context$context, useNA = "ifany"))

# ---- Enrichment test: are sig sites enriched in specific contexts? ----

# Filter to sites with context
rbcL_ctx <- rbcL_context[!is.na(context)]

message("\nSite class by structural context:")
print(table(rbcL_ctx$site_class, rbcL_ctx$context))

# Fisher's exact test for each context vs each site class
contexts <- unique(rbcL_ctx$context)
site_classes <- c("sig_no_control", "sig_with_control", "sig_both")

enrichment_results <- rbindlist(lapply(site_classes, function(sc) {
  rbindlist(lapply(contexts, function(ctx) {
    # 2x2 table: sig/not_sig vs in_context/not_in_context
    sig_in_ctx <- sum(rbcL_ctx$site_class == sc & rbcL_ctx$context == ctx)
    sig_not_ctx <- sum(rbcL_ctx$site_class == sc & rbcL_ctx$context != ctx)
    notsig_in_ctx <- sum(rbcL_ctx$site_class == "not_sig" & rbcL_ctx$context == ctx)
    notsig_not_ctx <- sum(rbcL_ctx$site_class == "not_sig" & rbcL_ctx$context != ctx)
    
    mat <- matrix(c(sig_in_ctx, sig_not_ctx, notsig_in_ctx, notsig_not_ctx), nrow = 2)
    ft <- fisher.test(mat)
    
    # Proportions
    prop_sig <- sig_in_ctx / (sig_in_ctx + sig_not_ctx)
    prop_notsig <- notsig_in_ctx / (notsig_in_ctx + notsig_not_ctx)
    
    data.table(
      site_class = sc,
      context = ctx,
      n_sig_in_ctx = sig_in_ctx,
      n_sig_total = sig_in_ctx + sig_not_ctx,
      prop_sig = prop_sig,
      prop_background = prop_notsig,
      fold_enrichment = prop_sig / prop_notsig,
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  }))
}))

enrichment_results[, p_adj := p.adjust(p_value, method = "BH")]
enrichment_results <- enrichment_results[order(p_adj)]

message("\n=== Context Enrichment Results ===")
print(enrichment_results[, .(site_class, context, n_sig_in_ctx, fold_enrichment, p_adj)])

# ---- Plot enrichment ----

enrichment_results[, sig_label := ifelse(p_adj < 0.001, "***",
                                         ifelse(p_adj < 0.01, "**",
                                                ifelse(p_adj < 0.05, "*", "")))]

p_enrich <- ggplot(enrichment_results, aes(x = context, y = fold_enrichment, fill = site_class)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = sig_label, y = fold_enrichment + 0.1), 
            position = position_dodge(0.8), vjust = 0, size = 4) +
  scale_fill_manual(values = c("sig_no_control" = "gold", 
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "GWAS hit enrichment by structural context",
       subtitle = "Fold enrichment vs not_sig background",
       x = "Structural context", y = "Fold enrichment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

p_enrich

# ---- Visualize sig sites on structure ----

# Add site class to spatial data
chain_L_df <- merge(chain_L_df, rcx_to_aln[, .(rcx_resno, aln_pos)], 
                    by.x = "resno", by.y = "rcx_resno", all.x = TRUE)
chain_L_df <- merge(chain_L_df, rbcL_context[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_L_df[is.na(site_class), site_class := "no_data"]

p_gwas_on_structure <- ggplot() +
  geom_point(data = ca_df[chain != "L"], aes(x = x, y = y), 
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_L_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_L_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_L_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_L_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_L_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers, aes(x = lig_x, y = lig_y),
             color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "GWAS significant sites on RbcL structure",
       subtitle = "Gold = sig_no_control, Blue = sig_with_control, Red = sig_both",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_gwas_on_structure


plot(rbcL_context$dist_to_active, rbcL_context$dist_to_L2)
plot(rbcL_context$dist_to_active, rbcL_context$dist_to_S)
plot(rbcL_context$dist_to_L2, rbcL_context$dist_to_S)

plot(rbcL_context$dist_to_active, -log10(rbcL_context$P_aa_with_pcs))
plot(rbcL_context$dist_to_active, -log10(rbcL_context$P_aa_only))

plot(rbcL_context$dist_to_L2, -log10(rbcL_context$P_aa_with_pcs))
plot(rbcL_context$dist_to_L2, -log10(rbcL_context$P_aa_only))

plot(rbcL_context$dist_to_S, -log10(rbcL_context$P_aa_with_pcs))
plot(rbcL_context$dist_to_S, -log10(rbcL_context$P_aa_only))

# ---- Continuous distance vs significance relationships ----

# Filter to sites with structural data
rbcL_ctx <- rbcL_context[!is.na(context)]

# Add -log10 p-values
rbcL_ctx[, neglog_p_only := -log10(P_aa_only)]
rbcL_ctx[, neglog_p_pcs := -log10(P_aa_with_pcs)]

p1 <- ggplot(rbcL_ctx, aes(x = dist_to_active, y = neglog_p_only)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(title = "Distance to active site vs -log10(P_aa_only)",
       x = "Distance to RuBP (Å)", y = "-log10(P)") +
  theme_classic()

p2 <- ggplot(rbcL_ctx, aes(x = dist_to_active, y = neglog_p_pcs)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "Distance to active site vs -log10(P_aa_with_pcs)",
       x = "Distance to RuBP (Å)", y = "-log10(P)") +
  theme_classic()

p3 <- ggplot(rbcL_ctx, aes(x = dist_to_L2, y = neglog_p_only)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(title = "Distance to L-L interface vs -log10(P_aa_only)",
       x = "Distance to adjacent RbcL (Å)", y = "-log10(P)") +
  theme_classic()

p4 <- ggplot(rbcL_ctx, aes(x = dist_to_L2, y = neglog_p_pcs)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "Distance to L-L interface vs -log10(P_aa_with_pcs)",
       x = "Distance to adjacent RbcL (Å)", y = "-log10(P)") +
  theme_classic()

p5 <- ggplot(rbcL_ctx, aes(x = dist_to_S, y = neglog_p_only)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  labs(title = "Distance to L-S interface vs -log10(P_aa_only)",
       x = "Distance to RbcS (Å)", y = "-log10(P)") +
  theme_classic()

p6 <- ggplot(rbcL_ctx, aes(x = dist_to_S, y = neglog_p_pcs)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = "Distance to L-S interface vs -log10(P_aa_with_pcs)",
       x = "Distance to RbcS (Å)", y = "-log10(P)") +
  theme_classic()

(p1 + p2) / (p3 + p4) / (p5 + p6)

# 2. Correlations
message("\n=== Correlations (Spearman) ===")
cor_results <- data.table(
  comparison = c("active_site vs P_only", "active_site vs P_pcs",
                 "L-L_interface vs P_only", "L-L_interface vs P_pcs",
                 "L-S_interface vs P_only", "L-S_interface vs P_pcs"),
  rho = c(
    cor(rbcL_ctx$dist_to_active, rbcL_ctx$neglog_p_only, method = "spearman", use = "complete"),
    cor(rbcL_ctx$dist_to_active, rbcL_ctx$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(rbcL_ctx$dist_to_L2, rbcL_ctx$neglog_p_only, method = "spearman", use = "complete"),
    cor(rbcL_ctx$dist_to_L2, rbcL_ctx$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(rbcL_ctx$dist_to_S, rbcL_ctx$neglog_p_only, method = "spearman", use = "complete"),
    cor(rbcL_ctx$dist_to_S, rbcL_ctx$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(rbcL_ctx$dist_to_active, rbcL_ctx$neglog_p_only, method = "spearman")$p.value,
    cor.test(rbcL_ctx$dist_to_active, rbcL_ctx$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(rbcL_ctx$dist_to_L2, rbcL_ctx$neglog_p_only, method = "spearman")$p.value,
    cor.test(rbcL_ctx$dist_to_L2, rbcL_ctx$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(rbcL_ctx$dist_to_S, rbcL_ctx$neglog_p_only, method = "spearman")$p.value,
    cor.test(rbcL_ctx$dist_to_S, rbcL_ctx$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_results)

# 3. Stratified by site class - colored scatterplots
p_strat1 <- ggplot(rbcL_ctx, aes(x = dist_to_active, y = neglog_p_only, color = site_class)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("not_sig" = "gray60", 
                                "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue",
                                "sig_both" = "darkred")) +
  labs(title = "P_aa_only by distance to active site",
       x = "Distance to RuBP (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "top")

p_strat2 <- ggplot(rbcL_ctx, aes(x = dist_to_L2, y = neglog_p_only, color = site_class)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("not_sig" = "gray60", 
                                "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue",
                                "sig_both" = "darkred")) +
  labs(title = "P_aa_only by distance to L-L interface",
       x = "Distance to adjacent RbcL (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_strat3 <- ggplot(rbcL_ctx, aes(x = dist_to_S, y = neglog_p_only, color = site_class)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("not_sig" = "gray60", 
                                "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue",
                                "sig_both" = "darkred")) +
  labs(title = "P_aa_only by distance to L-S interface",
       x = "Distance to RbcS (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_strat1 / (p_strat2 + p_strat3)

# 4. Boxplots: distance distributions by site class
p_box1 <- ggplot(rbcL_ctx, aes(x = site_class, y = dist_to_active, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", 
                               "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "Distance to active site", y = "Distance (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box2 <- ggplot(rbcL_ctx, aes(x = site_class, y = dist_to_L2, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", 
                               "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "Distance to L-L interface", y = "Distance (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box3 <- ggplot(rbcL_ctx, aes(x = site_class, y = dist_to_S, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", 
                               "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "Distance to L-S interface", y = "Distance (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box1 + p_box2 + p_box3

# 5. Wilcoxon tests: are sig sites at different distances?
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- rbcL_ctx[site_class == "not_sig"]

for (dist_var in c("dist_to_active", "dist_to_L2", "dist_to_S")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- rbcL_ctx[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑", "↓")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.3f",
                    sc, median(test_vals, na.rm = TRUE), direction, 
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}


# ---- Merge secondary structure with 3D context ----

# Get rbcL structural predictions from your netsurf analysis
# Recall: struct_var has consensus secondary structure by alignment position

# Check if we have rbcL in struct_var
message("Genes in struct_var:")
print(unique(struct_var$Gene))

# Filter to rbcL
rbcL_struct <- struct_var[Gene == "rbcL"]
message("rbcL positions with secondary structure: ", nrow(rbcL_struct))

# Merge with the 3D context data
rbcL_full <- merge(rbcL_ctx, rbcL_struct[, .(Position = Position, consensus, consensus_freq, 
                                             prop_H, prop_E, prop_C, entropy)],
                   by = "Position", all.x = TRUE)

message("Sites with both 3D and secondary structure: ", sum(!is.na(rbcL_full$consensus)))

# Cross-tabulate: secondary structure vs 3D context
message("\n=== Secondary structure by 3D context ===")
print(table(rbcL_full$consensus, rbcL_full$context, useNA = "ifany"))

# Cross-tabulate: secondary structure vs site class
message("\n=== Secondary structure by site class ===")
print(table(rbcL_full$consensus, rbcL_full$site_class, useNA = "ifany"))

# Proportions
message("\n=== Proportion in each secondary structure by site class ===")
ss_by_class <- rbcL_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class[, total := sum(N), by = site_class]
ss_by_class[, prop := N / total]
print(dcast(ss_by_class, site_class ~ consensus, value.var = "prop"))

# Fisher's test: C (coil) enrichment in sig sites?
message("\n=== Coil (C) enrichment tests ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(rbcL_full$site_class == sc & rbcL_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(rbcL_full$site_class == sc & rbcL_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(rbcL_full$site_class == "not_sig" & rbcL_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(rbcL_full$site_class == "not_sig" & rbcL_full$consensus != "C", na.rm = TRUE)
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / (sig_C + sig_notC)
  prop_bg <- bg_C / (bg_C + bg_notC)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualize: secondary structure vs 3D context ----

# Are coils preferentially at interfaces or in bulk?
p_ss_context <- ggplot(rbcL_full[!is.na(consensus) & !is.na(context)], 
                       aes(x = context, fill = consensus)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A"),
                    name = "Secondary\nStructure") +
  labs(title = "Secondary structure composition by 3D context",
       x = "Structural context", y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Are coils more/less significant?
p_ss_pval <- ggplot(rbcL_full[!is.na(consensus)], 
                    aes(x = consensus, y = -log10(P_aa_only), fill = consensus)) +
  geom_boxplot() +
  scale_fill_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A")) +
  labs(title = "P_aa_only by secondary structure", 
       x = "Secondary structure", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_ss_pval2 <- ggplot(rbcL_full[!is.na(consensus)], 
                     aes(x = consensus, y = -log10(P_aa_with_pcs), fill = consensus)) +
  geom_boxplot() +
  scale_fill_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A")) +
  labs(title = "P_aa_with_pcs by secondary structure",
       x = "Secondary structure", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_ss_context / (p_ss_pval + p_ss_pval2)

# ---- Combined view: site class colored by both SS and 3D ----

p_combined <- ggplot(rbcL_full[!is.na(consensus) & !is.na(context)],
                     aes(x = context, fill = site_class)) +
  geom_bar(position = "fill") +
  facet_wrap(~consensus, ncol = 3) +
  scale_fill_manual(values = c("not_sig" = "gray60",
                               "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "Site class distribution by secondary structure and 3D context",
       x = "3D context", y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_combined

# ---- rbcL summary ----

# Plot 1: GWAS sites on structure colored by significance
p_summary1 <- ggplot() +
  geom_point(data = ca_df[chain != "L"], aes(x = x, y = y), 
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_L_df[site_class == "no_data" | is.na(site_class)], 
             aes(x = x, y = y), color = "gray80", size = 1) +
  geom_point(data = chain_L_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5, alpha = 0.6) +
  geom_point(data = chain_L_df[site_class %in% c("sig_no_control", "sig_with_control", "sig_both")], 
             aes(x = x, y = y, color = site_class), size = 3) +
  geom_point(data = ligand_centers, aes(x = lig_x, y = lig_y),
             color = "black", size = 5, shape = 18) +
  scale_color_manual(values = c("sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", 
                                "sig_both" = "darkred"),
                     name = "Site class") +
  coord_fixed() +
  labs(title = "RbcL: GWAS significant sites in L8S8 complex",
       subtitle = "Black diamonds = active sites (RuBP); gray = L8S8 complex",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic() +
  theme(legend.position = "right")

# Plot 2: Distance distributions by site class
rbcL_ctx_long <- melt(rbcL_ctx[, .(site_class, dist_to_active, dist_to_L2, dist_to_S)],
                      id.vars = "site_class",
                      variable.name = "distance_type",
                      value.name = "distance")
rbcL_ctx_long[, distance_type := factor(distance_type, 
                                        levels = c("dist_to_active", "dist_to_L2", "dist_to_S"),
                                        labels = c("Active site", "L-L interface", "L-S interface"))]
rbcL_ctx_long[, site_class := factor(site_class,
                                     levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

p_summary2 <- ggplot(rbcL_ctx_long, aes(x = site_class, y = distance, fill = site_class)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~distance_type, scales = "free_y") +
  scale_fill_manual(values = c("not_sig" = "gray60",
                               "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue",
                               "sig_both" = "darkred")) +
  labs(title = "RbcL: Distance to functional regions by site class",
       subtitle = "No significant enrichment near active sites or interfaces",
       x = "", y = "Distance (Å)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p_summary1
p_summary2

# Summary message
message("\n", paste(rep("=", 70), collapse = ""), "\n")
message("SUMMARY: RbcL 3D Structure Analysis")
message(paste(rep("=", 70), collapse = ""), "\n")

message("SAMPLE SIZE:")
message(sprintf("  - Total GWAS sites mapped to structure: %d", nrow(rbcL_ctx)))
message(sprintf("  - sig_both: %d | sig_with_control: %d | sig_no_control: %d | not_sig: %d",
                sum(rbcL_ctx$site_class == "sig_both"),
                sum(rbcL_ctx$site_class == "sig_with_control"),
                sum(rbcL_ctx$site_class == "sig_no_control"),
                sum(rbcL_ctx$site_class == "not_sig")))

message("\nSPATIAL AUTOCORRELATION:")
message(sprintf("  - Moran's I = %.4f (expected = %.4f), p = %.3f", 
                obs_moran$I, obs_moran$E_I, p_value))
message("  - No significant spatial clustering of significant sites")

message("\nFUNCTIONAL REGION PROXIMITY:")
message("  - Active site: no enrichment (Wilcoxon p > 0.3 for all classes)")
message("  - L-L interface: marginal trend for sig_both being FARTHER (p = 0.056)")
message("  - L-S interface: no enrichment (p > 0.17)")

message("\nSECONDARY STRUCTURE:")
message("  - No significant coil (C) enrichment in sig sites")
message("  - Weak trend toward helix (H) enrichment in sig_with_control/sig_both (64-65% vs 44% bg)")

message("\nINTERPRETATION:")
message("  RbcL shows no clear spatial or functional enrichment pattern for GWAS hits.")
message("  This may reflect:")
message("    1. LOW POWER: Only 17 sig_both sites limits detection of subtle patterns")
message("    2. STRONG CONSTRAINT: RbcL is highly conserved; variable sites may be")
message("       uniformly distributed across 'permitted' positions")
message("    3. POPULATION STRUCTURE: sig_no_control sites (n=28) may reflect")
message("       phylogenetic signal rather than functional selection")
message("\n", paste(rep("=", 70), collapse = ""))
# ---- Check ATP synthase subunits in your data ----

# ATP synthase subunits are typically: atpA, atpB, atpE, atpF, atpH, atpI
atp_genes <- grep("^atp", unique(sites_df$Gene), value = TRUE)
message("ATP synthase genes in GWAS:")
print(atp_genes)

# Count sites per gene
atp_counts <- sites_df[Gene %in% atp_genes, .N, by = Gene]
# Fix: merge instead of direct assignment
sig_both_counts <- sites_df[Gene %in% atp_genes & 
                              P_aa_only < thresh_aa_only & 
                              P_aa_with_pcs < thresh_aa_pcs, .N, by = Gene]
setnames(sig_both_counts, "N", "n_sig_both")

atp_counts <- merge(atp_counts, sig_both_counts, by = "Gene", all.x = TRUE)
atp_counts[is.na(n_sig_both), n_sig_both := 0]
print(atp_counts)
# Also check: which genes had strongest C enrichment in your original analysis?
# Let's look at gene-level enrichment from high_cons
gene_ss_enrichment <- high_cons[, .(
  n_total = .N,
  n_sig_both = sum(site_class == "sig_both"),
  prop_C_sig = mean(consensus == "C" & site_class == "sig_both", na.rm = TRUE) / mean(site_class == "sig_both"),
  prop_C_bg = mean(consensus == "C" & site_class == "not_sig", na.rm = TRUE) / mean(site_class == "not_sig")
), by = Gene]

gene_ss_enrichment[, C_fold := prop_C_sig / prop_C_bg]
gene_ss_enrichment <- gene_ss_enrichment[n_sig_both >= 5]  # filter to genes with enough sig sites
gene_ss_enrichment <- gene_ss_enrichment[order(-C_fold)]

message("\nGenes with highest coil enrichment in sig_both (n >= 5):")
print(head(gene_ss_enrichment, 15))


# Add other sig classes
sig_no_ctrl <- sites_df[Gene %in% atp_genes & P_aa_only < thresh_aa_only, .N, by = Gene]
setnames(sig_no_ctrl, "N", "n_sig_no_ctrl")

sig_with_ctrl <- sites_df[Gene %in% atp_genes & P_aa_with_pcs < thresh_aa_pcs, .N, by = Gene]
setnames(sig_with_ctrl, "N", "n_sig_with_ctrl")

atp_counts <- merge(atp_counts, sig_no_ctrl, by = "Gene", all.x = TRUE)
atp_counts <- merge(atp_counts, sig_with_ctrl, by = "Gene", all.x = TRUE)
atp_counts[is.na(n_sig_no_ctrl), n_sig_no_ctrl := 0]
atp_counts[is.na(n_sig_with_ctrl), n_sig_with_ctrl := 0]

print(atp_counts[order(-n_sig_both)])

# Check secondary structure enrichment for ATP genes from your original analysis
message("\n=== Secondary structure in ATP genes (from high_cons) ===")
atp_ss <- high_cons[Gene %in% atp_genes, .(
  n_total = .N,
  n_sig_both = sum(site_class == "sig_both"),
  prop_C_all = mean(consensus == "C"),
  prop_C_sig = mean(consensus[site_class == "sig_both"] == "C"),
  prop_C_bg = mean(consensus[site_class == "not_sig"] == "C")
), by = Gene]
atp_ss[, C_fold := prop_C_sig / prop_C_bg]
print(atp_ss[order(-n_sig_both)])

# ---- Get ATP synthase structure ----

# Search for a good plant/chloroplast ATP synthase structure
# Recent cryo-EM structures are excellent

# Let's try spinach chloroplast ATP synthase (there are good recent structures)
# PDB 6FKF, 6FKH, or 6VOX are spinach CF1-CF0

pdb_atp <- read.pdb("6FKF")
print(pdb_atp)

message("\nChains in ATP synthase structure:")
print(table(pdb_atp$atom$chain))

# Check what's in there
message("\nUnique residue types (to find ligands):")
print(unique(pdb_atp$atom$resid[!pdb_atp$atom$resid %in% c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")]))

# ---- Identify ATP synthase subunits ----

# Get chain info
chain_info <- pdb_atp$atom[pdb_atp$atom$elety == "CA", ]
chain_summary <- as.data.table(chain_info)[, .(
  n_residues = .N,
  first_resno = min(resno),
  last_resno = max(resno)
), by = chain]
chain_summary <- chain_summary[order(-n_residues)]
print(chain_summary)

# The large chains should be alpha (atpA) and beta (atpB)
# Beta is the catalytic subunit - let's check sequence to identify
# Extract sequence from largest chains

large_chains <- chain_summary[n_residues > 400]$chain
message("\nLarge chains (likely α/β): ", paste(large_chains, collapse = ", "))

# Get sequence of chain A
chain_A <- trim.pdb(pdb_atp, chain = "A")
ca_A <- atom.select(chain_A, elety = "CA")
chain_A_ca <- trim.pdb(chain_A, ca_A)
seq_A <- paste(aa_map[chain_A_ca$atom$resid], collapse = "")
message("\nChain A length: ", nchar(seq_A))
message("Chain A first 50 aa: ", substr(seq_A, 1, 50))

# Get sequence of chain B
chain_B_atp <- trim.pdb(pdb_atp, chain = "B")
ca_B <- atom.select(chain_B_atp, elety = "CA")
chain_B_ca <- trim.pdb(chain_B_atp, ca_B)
seq_B <- paste(aa_map[chain_B_ca$atom$resid], collapse = "")
message("\nChain B length: ", nchar(seq_B))
message("Chain B first 50 aa: ", substr(seq_B, 1, 50))

# Align to your atpB alignment to identify which is beta
atpB_aln <- readAAStringSet("data/tmp/alignedGenes/atpB_AA_aligned.fasta")
message("\natpB alignment: ", length(atpB_aln), " sequences, length ", width(atpB_aln)[1])

# Get Arabidopsis atpB sequence
at_idx_atpB <- grep(arabidopsis_th_id, names(atpB_aln))
stopifnot(length(at_idx_atpB) == 1)
at_atpB_seq <- gsub("-", "", as.character(atpB_aln[[at_idx_atpB]]))
message("Arabidopsis atpB length: ", nchar(at_atpB_seq))

# Align PDB chains to Arabidopsis atpB
pw_A <- pairwiseAlignment(seq_A, at_atpB_seq, type = "local")
pw_B <- pairwiseAlignment(seq_B, at_atpB_seq, type = "local")

message("\nChain A vs atpB: score = ", round(score(pw_A), 1), ", identity = ", round(pid(pw_A), 1), "%")
message("Chain B vs atpB: score = ", round(score(pw_B), 1), ", identity = ", round(pid(pw_B), 1), "%")

# ---- Map atpB GWAS sites to structure ----

# Chain B = atpB (β subunit)
# Extract coordinates
chain_B_atp <- trim.pdb(pdb_atp, chain = "B")
ca_B <- atom.select(chain_B_atp, elety = "CA")
chain_B_ca <- trim.pdb(chain_B_atp, ca_B)

coords_atpB <- as.data.table(chain_B_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpB, c("pdb_resno", "aa", "x", "y", "z"))
message("atpB (chain B) residues: ", nrow(coords_atpB))

# Get sequence
seq_B <- paste(aa_map[coords_atpB$aa], collapse = "")

# Create position mapping: PDB -> Arabidopsis -> alignment
pw_B <- pairwiseAlignment(seq_B, at_atpB_seq, type = "global")
print(pw_B)

# Build mapping
aln_pattern_B <- as.character(pattern(pw_B))  # PDB seq (with gaps)
aln_subject_B <- as.character(subject(pw_B))  # Arabidopsis (with gaps)

pdb_pos <- 0L
at_pos <- 0L
pdb_to_at_atpB <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_B)) {
  p_char <- substr(aln_pattern_B, i, i)
  s_char <- substr(aln_subject_B, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpB <- rbind(pdb_to_at_atpB, 
                            data.table(pdb_resno = coords_atpB$pdb_resno[pdb_pos], 
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpB), " positions between PDB and Arabidopsis")

# Now map Arabidopsis ungapped -> alignment position
at_atpB_aln_seq <- as.character(atpB_aln[[at_idx_atpB]])
aln_chars_atpB <- strsplit(at_atpB_aln_seq, "")[[1]]

at_to_aln_atpB <- data.table(
  aln_pos = seq_along(aln_chars_atpB),
  aln_char = aln_chars_atpB,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpB)) {
  if (aln_chars_atpB[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpB$at_ungapped[i] <- ungapped
  }
}

# Merge: PDB -> Arabidopsis ungapped -> alignment position
pdb_to_aln_atpB <- merge(pdb_to_at_atpB, at_to_aln_atpB[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
message("Mapped ", nrow(pdb_to_aln_atpB), " PDB positions to alignment positions")

# Add 3D coordinates
pdb_to_aln_atpB <- merge(pdb_to_aln_atpB, coords_atpB[, .(pdb_resno, x, y, z)], by = "pdb_resno")

# ---- Get atpB GWAS sites ----
atpB_sites <- sites_df[Gene == "atpB"]
atpB_sites <- merge(atpB_sites, pdb_to_aln_atpB[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpB GWAS sites: ", nrow(atpB_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpB_sites$x)))

# Classify sites
atpB_sites[, site_class := "not_sig"]
atpB_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpB_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpB_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpB_sites$site_class))

# Filter to sites with 3D coords
atpB_3d <- atpB_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpB_3d$site_class))

# ---- EDA: ATP synthase structure ----

# Get all CA atoms for the full complex
all_ca_atp <- atom.select(pdb_atp, elety = "CA")
all_ca_atp_pdb <- trim.pdb(pdb_atp, all_ca_atp)
ca_df_atp <- as.data.table(all_ca_atp_pdb$atom[, c("chain", "resno", "x", "y", "z")])

# Classify chains by subunit type
ca_df_atp[, subunit := fcase(
  chain %in% c("A", "C", "E"), "alpha",
  chain %in% c("B", "D", "F"), "beta",
  chain %in% c("G","H","I","J","K","L","M","N","O","P","Q","R","S","T"), "c-ring",
  chain == "g", "gamma",
  chain == "a", "a",
  chain == "b", "b",
  chain == "d", "delta",
  chain == "e", "epsilon",
  chain == "p", "OSCP",
  default = "other"
)]

message("Subunit composition:")
print(table(ca_df_atp$subunit))

# Plot 1: XY projection - top view of F1 head
p_atp_xy <- ggplot(ca_df_atp, aes(x = x, y = y, color = subunit)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("alpha" = "steelblue", "beta" = "coral",
                                "c-ring" = "gray70", "gamma" = "darkgreen",
                                "a" = "purple", "b" = "orange", 
                                "delta" = "brown", "epsilon" = "pink",
                                "OSCP" = "cyan", "other" = "gray50")) +
  coord_fixed() +
  labs(title = "ATP synthase - XY projection (top view)",
       subtitle = "α (blue), β (coral), γ (green), c-ring (gray)",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

# Plot 2: XZ projection - side view
p_atp_xz <- ggplot(ca_df_atp, aes(x = x, y = z, color = subunit)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("alpha" = "steelblue", "beta" = "coral",
                                "c-ring" = "gray70", "gamma" = "darkgreen",
                                "a" = "purple", "b" = "orange",
                                "delta" = "brown", "epsilon" = "pink",
                                "OSCP" = "cyan", "other" = "gray50")) +
  coord_fixed() +
  labs(title = "ATP synthase - XZ projection (side view)",
       subtitle = "F1 head on top, c-ring below",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atp_xy + p_atp_xz

# Plot 3: Just the F1 head (α3β3) - XY view with chain B highlighted
f1_chains <- c("A", "B", "C", "D", "E", "F")
ca_f1 <- ca_df_atp[chain %in% f1_chains]

p_f1_xy <- ggplot() +
  geom_point(data = ca_f1[chain != "B"], aes(x = x, y = y, color = subunit), 
             size = 1, alpha = 0.4) +
  geom_point(data = ca_f1[chain == "B"], aes(x = x, y = y), 
             color = "darkred", size = 1.5) +
  scale_color_manual(values = c("alpha" = "steelblue", "beta" = "coral")) +
  coord_fixed() +
  labs(title = "F1 head (α3β3) - chain B (atpB) highlighted",
       subtitle = "Dark red = chain B; coral = other β; blue = α",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

# Plot 4: Find and show ATP/ADP binding sites
ligand_atoms_atp <- pdb_atp$atom[pdb_atp$atom$resid %in% c("ATP", "ADP", "MG"), ]
ligand_centers_atp <- as.data.table(ligand_atoms_atp)[, .(
  lig_x = mean(x), lig_y = mean(y), lig_z = mean(z),
  ligand = unique(resid)
), by = chain]

message("\nLigands by chain:")
print(ligand_centers_atp)

p_f1_ligands <- ggplot() +
  geom_point(data = ca_f1, aes(x = x, y = y, color = subunit), size = 0.8, alpha = 0.4) +
  geom_point(data = ligand_centers_atp[ligand == "ATP"], aes(x = lig_x, y = lig_y),
             color = "red", size = 5, shape = 18) +
  geom_point(data = ligand_centers_atp[ligand == "ADP"], aes(x = lig_x, y = lig_y),
             color = "orange", size = 5, shape = 17) +
  geom_point(data = ligand_centers_atp[ligand == "MG"], aes(x = lig_x, y = lig_y),
             color = "darkgreen", size = 3, shape = 16) +
  scale_color_manual(values = c("alpha" = "steelblue", "beta" = "coral")) +
  coord_fixed() +
  labs(title = "F1 head with nucleotide binding sites",
       subtitle = "Red diamond = ATP, Orange triangle = ADP, Green = Mg",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

(p_atp_xy + p_atp_xz) / (p_f1_xy + p_f1_ligands)

# ---- Calculate structural context for atpB (chain B) ----

chain_B_df <- ca_df_atp[chain == "B"]

# Distance to nucleotide binding site (ATP/ADP in chain B or closest)
# Find which ligand is in/near chain B
message("\nFinding closest nucleotide to each residue in chain B...")

# Get all nucleotide positions
nuc_atoms <- pdb_atp$atom[pdb_atp$atom$resid %in% c("ATP", "ADP"), ]
nuc_coords <- as.data.table(nuc_atoms)[, .(x = mean(x), y = mean(y), z = mean(z)), by = .(chain, resid)]

# Distance to closest nucleotide for each chain B residue
chain_B_df[, dist_to_nucleotide := sapply(1:.N, function(i) {
  min(sqrt((x[i] - nuc_coords$x)^2 + (y[i] - nuc_coords$y)^2 + (z[i] - nuc_coords$z)^2))
})]

# Distance to alpha subunits (β-α interface)
alpha_coords <- ca_df_atp[subunit == "alpha"]
chain_B_df[, dist_to_alpha := sapply(1:.N, function(i) {
  min(sqrt((x[i] - alpha_coords$x)^2 + (y[i] - alpha_coords$y)^2 + (z[i] - alpha_coords$z)^2))
})]

# Distance to gamma (central stalk - β-γ interface)
gamma_coords <- ca_df_atp[subunit == "gamma"]
chain_B_df[, dist_to_gamma := sapply(1:.N, function(i) {
  min(sqrt((x[i] - gamma_coords$x)^2 + (y[i] - gamma_coords$y)^2 + (z[i] - gamma_coords$z)^2))
})]

# Distance to other beta subunits (β-β interface, though less common)
other_beta <- ca_df_atp[subunit == "beta" & chain != "B"]
chain_B_df[, dist_to_other_beta := sapply(1:.N, function(i) {
  min(sqrt((x[i] - other_beta$x)^2 + (y[i] - other_beta$y)^2 + (z[i] - other_beta$z)^2))
})]

# Plot distance profiles
p_dist_nuc <- ggplot(chain_B_df, aes(x = resno, y = dist_to_nucleotide)) +
  geom_line() +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  labs(title = "Distance to nucleotide (ATP/ADP)", x = "Residue", y = "Distance (Å)") +
  theme_classic()

p_dist_alpha <- ggplot(chain_B_df, aes(x = resno, y = dist_to_alpha)) +
  geom_line() +
  geom_hline(yintercept = 8, color = "blue", linetype = "dashed") +
  labs(title = "Distance to α subunits", x = "Residue", y = "Distance (Å)") +
  theme_classic()

p_dist_gamma <- ggplot(chain_B_df, aes(x = resno, y = dist_to_gamma)) +
  geom_line() +
  geom_hline(yintercept = 8, color = "green", linetype = "dashed") +
  labs(title = "Distance to γ (central stalk)", x = "Residue", y = "Distance (Å)") +
  theme_classic()

(p_dist_nuc + p_dist_alpha) / p_dist_gamma

# ---- Classify structural context ----
chain_B_df[, context := "bulk"]
chain_B_df[dist_to_nucleotide < 15, context := "near_nucleotide"]
chain_B_df[dist_to_alpha < 8, context := "alpha_interface"]
chain_B_df[dist_to_gamma < 8, context := "gamma_interface"]
chain_B_df[dist_to_nucleotide < 15 & dist_to_alpha < 8, context := "catalytic_interface"]

message("\n=== atpB structural context ===")
print(table(chain_B_df$context))

# Visualize context on structure
p_context_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = x, y = y), 
             color = "gray85", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df, aes(x = x, y = y, color = context), size = 2) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")], 
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  scale_color_manual(values = c("bulk" = "gray50",
                                "near_nucleotide" = "orange",
                                "alpha_interface" = "blue",
                                "gamma_interface" = "green",
                                "catalytic_interface" = "red")) +
  coord_fixed() +
  labs(title = "atpB (chain B) structural context",
       subtitle = "Black diamonds = ATP/ADP binding sites",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_context_xy

# Context along sequence
p_context_seq <- ggplot(chain_B_df, aes(x = resno, y = 0)) +
  geom_point(aes(color = context), size = 2, shape = 15) +
  scale_color_manual(values = c("bulk" = "gray50",
                                "near_nucleotide" = "orange",
                                "alpha_interface" = "blue",
                                "gamma_interface" = "green",
                                "catalytic_interface" = "red")) +
  labs(title = "atpB structural context along sequence", x = "Residue number") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), legend.position = "top")

p_context_seq

# ---- Check nucleotide binding per chain ----

message("Nucleotides by chain:")
print(ligand_centers_atp[order(chain)])

# Which beta chains have nucleotides?
beta_chains <- c("B", "D", "F")
message("\nBeta subunit nucleotide status:")
for (ch in beta_chains) {
  nucs <- ligand_centers_atp[chain == ch & ligand %in% c("ATP", "ADP")]
  if (nrow(nucs) > 0) {
    message("  Chain ", ch, ": ", paste(nucs$ligand, collapse = ", "))
  } else {
    message("  Chain ", ch, ": EMPTY (βE conformation)")
  }
}

# So for chain B, we should calculate distance to the closest nucleotide
# which would be in the adjacent alpha or another beta subunit
# This represents the "nucleotide binding pocket" region even if empty

# Let's find the closest nucleotide binding site to chain B
# (likely in an adjacent alpha subunit)
closest_nuc_to_B <- nuc_coords[, .(
  mean_dist_to_B = mean(sapply(1:nrow(chain_B_df), function(i) {
    sqrt((chain_B_df$x[i] - x)^2 + (chain_B_df$y[i] - y)^2 + (chain_B_df$z[i] - z)^2)
  }))
), by = .(chain, resid)]
print(closest_nuc_to_B[order(mean_dist_to_B)])

# For functional analysis, we should look at:
# 1. The Walker A motif (P-loop, binds phosphate) - conserved GxxxxGKT
# 2. The Walker B motif (binds Mg) - conserved hhhhDE
# 3. The catalytic site between β and α

# Let's identify the nucleotide binding pocket by finding where
# the OTHER beta subunits bind ATP - that pocket location in chain B
# would be the "empty" catalytic site

# Get ATP position from chain D or F (the ones with ATP)
atp_position <- ligand_centers_atp[ligand == "ATP"][1]  # take first ATP
message("\nUsing ATP position from chain ", atp_position$chain, " as reference")
message("ATP coordinates: ", round(atp_position$lig_x, 1), ", ", 
        round(atp_position$lig_y, 1), ", ", round(atp_position$lig_z, 1))

# Now recalculate distance to this "reference" catalytic site position
# But actually - each β subunit is rotated ~120° around the central axis
# So we need to find the EQUIVALENT position in chain B

# Alternative approach: use distance to closest nucleotide from ANY chain
# This captures "proximity to a nucleotide binding region"
chain_B_df[, dist_to_any_nucleotide := sapply(1:.N, function(i) {
  min(sqrt((x[i] - nuc_coords$x)^2 + (y[i] - nuc_coords$y)^2 + (z[i] - nuc_coords$z)^2))
})]

message("\nDistance to any nucleotide - summary for chain B:")
print(summary(chain_B_df$dist_to_any_nucleotide))

# Plot comparison
p_nuc_dist <- ggplot(chain_B_df, aes(x = resno)) +
  geom_line(aes(y = dist_to_any_nucleotide), color = "red") +
  geom_hline(yintercept = 15, linetype = "dashed") +
  labs(title = "Chain B distance to nearest nucleotide (any chain)",
       subtitle = "Chain B is empty (βE) - closest nucleotides are in adjacent subunits",
       x = "Residue", y = "Distance (Å)") +
  theme_classic()

p_nuc_dist

# Update context classification using distance to any nucleotide
chain_B_df[, context := "bulk"]
chain_B_df[dist_to_any_nucleotide < 15, context := "near_nucleotide"]
chain_B_df[dist_to_alpha < 8, context := "alpha_interface"]
chain_B_df[dist_to_gamma < 8, context := "gamma_interface"]
chain_B_df[dist_to_any_nucleotide < 15 & dist_to_alpha < 8, context := "catalytic_interface"]

message("\n=== Updated atpB structural context ===")
print(table(chain_B_df$context))

# Visualize updated context
p_context_updated <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = x, y = y), 
             color = "gray85", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df, aes(x = x, y = y, color = context), size = 2) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")], 
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  scale_color_manual(values = c("bulk" = "gray50",
                                "near_nucleotide" = "orange",
                                "alpha_interface" = "blue",
                                "gamma_interface" = "green",
                                "catalytic_interface" = "red")) +
  coord_fixed() +
  labs(title = "atpB (chain B = βE empty state) structural context",
       subtitle = "Black diamonds = ATP/ADP in adjacent subunits",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_context_updated

# ---- Merge GWAS sites with structural context ----

# Add PDB resno to chain_B_df for merging
chain_B_df[, pdb_resno := resno]

# Merge structural context with GWAS sites
atpB_3d <- merge(atpB_3d, 
                 chain_B_df[, .(pdb_resno, context, dist_to_any_nucleotide, 
                                dist_to_alpha, dist_to_gamma, dist_to_other_beta)],
                 by = "pdb_resno", all.x = TRUE)

message("atpB sites with structural context: ", sum(!is.na(atpB_3d$context)), " of ", nrow(atpB_3d))

message("\n=== Site class by structural context ===")
print(table(atpB_3d$site_class, atpB_3d$context, useNA = "ifany"))

# ---- Enrichment tests ----

contexts <- unique(atpB_3d$context[!is.na(atpB_3d$context)])
site_classes <- c("sig_no_control", "sig_with_control", "sig_both")

enrichment_atpB <- rbindlist(lapply(site_classes, function(sc) {
  rbindlist(lapply(contexts, function(ctx) {
    sig_in_ctx <- sum(atpB_3d$site_class == sc & atpB_3d$context == ctx, na.rm = TRUE)
    sig_not_ctx <- sum(atpB_3d$site_class == sc & atpB_3d$context != ctx, na.rm = TRUE)
    notsig_in_ctx <- sum(atpB_3d$site_class == "not_sig" & atpB_3d$context == ctx, na.rm = TRUE)
    notsig_not_ctx <- sum(atpB_3d$site_class == "not_sig" & atpB_3d$context != ctx, na.rm = TRUE)
    
    mat <- matrix(c(sig_in_ctx, sig_not_ctx, notsig_in_ctx, notsig_not_ctx), nrow = 2)
    ft <- fisher.test(mat)
    
    prop_sig <- sig_in_ctx / max((sig_in_ctx + sig_not_ctx), 1)
    prop_notsig <- notsig_in_ctx / max((notsig_in_ctx + notsig_not_ctx), 1)
    
    data.table(
      site_class = sc,
      context = ctx,
      n_sig_in_ctx = sig_in_ctx,
      n_sig_total = sig_in_ctx + sig_not_ctx,
      prop_sig = prop_sig,
      prop_background = prop_notsig,
      fold_enrichment = ifelse(prop_notsig > 0, prop_sig / prop_notsig, NA),
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  }))
}))

enrichment_atpB[, p_adj := p.adjust(p_value, method = "BH")]
enrichment_atpB <- enrichment_atpB[order(p_adj)]

message("\n=== Context Enrichment Results ===")
print(enrichment_atpB[, .(site_class, context, n_sig_in_ctx, n_sig_total, fold_enrichment, p_value, p_adj)])

# ---- Continuous distance analysis ----

atpB_3d[, neglog_p_only := -log10(P_aa_only)]
atpB_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

# Correlations
message("\n=== Spearman correlations ===")
cor_atpB <- data.table(
  comparison = c("nucleotide vs P_only", "nucleotide vs P_pcs",
                 "alpha_interface vs P_only", "alpha_interface vs P_pcs",
                 "gamma_interface vs P_only", "gamma_interface vs P_pcs"),
  rho = c(
    cor(atpB_3d$dist_to_any_nucleotide, atpB_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpB_3d$dist_to_any_nucleotide, atpB_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpB_3d$dist_to_alpha, atpB_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpB_3d$dist_to_alpha, atpB_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpB_3d$dist_to_gamma, atpB_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpB_3d$dist_to_gamma, atpB_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpB_3d$dist_to_any_nucleotide, atpB_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpB_3d$dist_to_any_nucleotide, atpB_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpB_3d$dist_to_alpha, atpB_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpB_3d$dist_to_alpha, atpB_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpB_3d$dist_to_gamma, atpB_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpB_3d$dist_to_gamma, atpB_3d$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_atpB)

# Scatterplots with loess
p_nuc_only <- ggplot(atpB_3d[!is.na(dist_to_any_nucleotide)], 
                     aes(x = dist_to_any_nucleotide, y = neglog_p_only)) +
  geom_point(aes(color = site_class), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to nucleotide vs P_aa_only",
       x = "Distance to nearest ATP/ADP (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_nuc_pcs <- ggplot(atpB_3d[!is.na(dist_to_any_nucleotide)], 
                    aes(x = dist_to_any_nucleotide, y = neglog_p_pcs)) +
  geom_point(aes(color = site_class), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to nucleotide vs P_aa_with_pcs",
       x = "Distance to nearest ATP/ADP (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_alpha_only <- ggplot(atpB_3d[!is.na(dist_to_alpha)], 
                       aes(x = dist_to_alpha, y = neglog_p_only)) +
  geom_point(aes(color = site_class), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to α subunit vs P_aa_only",
       x = "Distance to α (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_gamma_only <- ggplot(atpB_3d[!is.na(dist_to_gamma)], 
                       aes(x = dist_to_gamma, y = neglog_p_only)) +
  geom_point(aes(color = site_class), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                                "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to γ (central stalk) vs P_aa_only",
       x = "Distance to γ (Å)", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "top")

(p_nuc_only + p_nuc_pcs) / (p_alpha_only + p_gamma_only)

# ---- Wilcoxon tests: distances by site class ----

message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpB_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_any_nucleotide", "dist_to_alpha", "dist_to_gamma")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpB_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Boxplots ----

p_box_nuc <- ggplot(atpB_3d[!is.na(context)], aes(x = site_class, y = dist_to_any_nucleotide, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to nucleotide", y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_alpha <- ggplot(atpB_3d[!is.na(context)], aes(x = site_class, y = dist_to_alpha, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to α interface", y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_gamma <- ggplot(atpB_3d[!is.na(context)], aes(x = site_class, y = dist_to_gamma, fill = site_class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to γ interface", y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_nuc + p_box_alpha + p_box_gamma

# ---- Visualize sig sites on structure ----

# Merge site class onto chain_B_df for plotting
chain_B_df <- merge(chain_B_df, pdb_to_aln_atpB[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_B_df <- merge(chain_B_df, atpB_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_B_df[is.na(site_class), site_class := "no_data"]

p_gwas_atpB <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_B_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_B_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "atpB GWAS hits on ATP synthase structure",
       subtitle = "Gold = sig_no_control, Blue = sig_with_control, Red = sig_both; Diamonds = ATP/ADP",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_gwas_atpB

# ---- Boxplots with jitter and significance annotations ----

# Store p-values for annotation
pvals_nuc <- c(not_sig = NA, sig_no_control = 0.0003, sig_with_control = 0.0862, sig_both = 0.3717)
pvals_alpha <- c(not_sig = NA, sig_no_control = 0.0031, sig_with_control = 0.6535, sig_both = 0.0364)
pvals_gamma <- c(not_sig = NA, sig_no_control = 0.5913, sig_with_control = 0.9538, sig_both = 0.0286)

# Helper function for significance stars
pval_to_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# Create annotation data frames
annot_nuc <- data.table(
  site_class = factor(c("sig_no_control", "sig_with_control", "sig_both"),
                      levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both")),
  label = pval_to_stars(pvals_nuc[c("sig_no_control", "sig_with_control", "sig_both")]),
  y = max(atpB_3d$dist_to_any_nucleotide, na.rm = TRUE) + 5
)

annot_alpha <- data.table(
  site_class = factor(c("sig_no_control", "sig_with_control", "sig_both"),
                      levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both")),
  label = pval_to_stars(pvals_alpha[c("sig_no_control", "sig_with_control", "sig_both")]),
  y = max(atpB_3d$dist_to_alpha, na.rm = TRUE) + 3
)

annot_gamma <- data.table(
  site_class = factor(c("sig_no_control", "sig_with_control", "sig_both"),
                      levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both")),
  label = pval_to_stars(pvals_gamma[c("sig_no_control", "sig_with_control", "sig_both")]),
  y = max(atpB_3d$dist_to_gamma, na.rm = TRUE) + 3
)

# Factor site_class for proper ordering
atpB_3d[, site_class := factor(site_class, 
                               levels = c("not_sig", "sig_no_control", "sig_with_control", "sig_both"))]

p_box_nuc <- ggplot(atpB_3d[!is.na(dist_to_any_nucleotide)], 
                    aes(x = site_class, y = dist_to_any_nucleotide, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_text(data = annot_nuc, aes(x = site_class, y = y, label = label), size = 5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to nucleotide", 
       subtitle = "sig_no_control farther from active site",
       y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_alpha <- ggplot(atpB_3d[!is.na(dist_to_alpha)], 
                      aes(x = site_class, y = dist_to_alpha, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_text(data = annot_alpha, aes(x = site_class, y = y, label = label), size = 5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to α interface",
       subtitle = "sig_both CLOSER to α subunits",
       y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_gamma <- ggplot(atpB_3d[!is.na(dist_to_gamma)], 
                      aes(x = site_class, y = dist_to_gamma, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_text(data = annot_gamma, aes(x = site_class, y = y, label = label), size = 5) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "Distance to γ (central stalk)",
       subtitle = "sig_both CLOSER to rotary mechanism",
       y = "Distance (Å)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_box_nuc + p_box_alpha + p_box_gamma


# ---- Three projections of atpB GWAS sites ----

p_gwas_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_B_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_B_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "XY projection (top view)",
       subtitle = "Looking down central axis",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_gwas_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_B_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_B_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_z), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "XZ projection (side view)",
       subtitle = "F1 head top, c-ring bottom",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_gwas_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "B"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_B_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_B_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_B_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_B_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_y, y = lig_z), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "YZ projection (side view 2)",
       subtitle = "Orthogonal side view",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_gwas_xy / (p_gwas_xz + p_gwas_yz)

p_gwas_xy
p_gwas_x

p_gwas_yz

# ---- Merge secondary structure with atpB 3D data ----

# Get atpB secondary structure from high_cons
atpB_struct <- struct_var[Gene == "atpB"]
message("atpB positions with secondary structure: ", nrow(atpB_struct))

# Merge with 3D data
atpB_full <- merge(atpB_3d, atpB_struct[, .(Position, consensus, consensus_freq, 
                                            prop_H, prop_E, prop_C)],
                   by = "Position", all.x = TRUE)

message("Sites with both 3D and secondary structure: ", sum(!is.na(atpB_full$consensus)))

# Cross-tabulate
message("\n=== Secondary structure by site class ===")
print(table(atpB_full$consensus, atpB_full$site_class, useNA = "ifany"))

# Proportions
message("\n=== Proportion in each secondary structure by site class ===")
ss_by_class_atpB <- atpB_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpB[, total := sum(N), by = site_class]
ss_by_class_atpB[, prop := N / total]
print(dcast(ss_by_class_atpB, site_class ~ consensus, value.var = "prop"))

# Fisher's test for coil enrichment
message("\n=== Coil (C) enrichment in atpB ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpB_full$site_class == sc & atpB_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpB_full$site_class == sc & atpB_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpB_full$site_class == "not_sig" & atpB_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpB_full$site_class == "not_sig" & atpB_full$consensus != "C", na.rm = TRUE)
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# Is coil enrichment related to interface proximity?
message("\n=== Secondary structure by 3D context ===")
print(table(atpB_full$consensus, atpB_full$context, useNA = "ifany"))

# Are coils at the α interface?
message("\n=== Proportion coil by context ===")
context_ss <- atpB_full[!is.na(consensus) & !is.na(context), .(
  n = .N,
  prop_C = mean(consensus == "C"),
  prop_H = mean(consensus == "H"),
  prop_E = mean(consensus == "E")
), by = context]
print(context_ss)

# Plot: secondary structure composition by context
p_ss_context_atpB <- ggplot(atpB_full[!is.na(consensus) & !is.na(context)],
                            aes(x = context, fill = consensus)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A")) +
  labs(title = "atpB: Secondary structure by 3D context",
       x = "Structural context", y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot: P-value by secondary structure
p_ss_pval_atpB <- ggplot(atpB_full[!is.na(consensus)],
                         aes(x = consensus, y = -log10(P_aa_with_pcs), fill = consensus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c("H" = "#E41A1C", "E" = "#377EB8", "C" = "#4DAF4A")) +
  labs(title = "atpB: P_aa_with_pcs by secondary structure",
       x = "Secondary structure", y = "-log10(P)") +
  theme_classic() +
  theme(legend.position = "none")

p_ss_context_atpB + p_ss_pval_atpB

#

# ---- atpA (α subunit) - Chain A ----

# Extract chain A coordinates
chain_A_atp <- trim.pdb(pdb_atp, chain = "A")
ca_A <- atom.select(chain_A_atp, elety = "CA")
chain_A_ca <- trim.pdb(chain_A_atp, ca_A)

coords_atpA <- as.data.table(chain_A_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpA, c("pdb_resno", "aa", "x", "y", "z"))
message("atpA (chain A) residues: ", nrow(coords_atpA))

# Get sequence
seq_A <- paste(aa_map[coords_atpA$aa], collapse = "")
message("Chain A sequence length: ", nchar(seq_A))

# Load atpA alignment and get Arabidopsis sequence
atpA_aln <- readAAStringSet("data/tmp/alignedGenes/atpA_AA_aligned.fasta")
message("atpA alignment: ", length(atpA_aln), " sequences, length ", width(atpA_aln)[1])

at_idx_atpA <- grep(arabidopsis_th_id, names(atpA_aln))
stopifnot(length(at_idx_atpA) == 1)
at_atpA_seq <- gsub("-", "", as.character(atpA_aln[[at_idx_atpA]]))
message("Arabidopsis atpA length: ", nchar(at_atpA_seq))

# Align PDB to Arabidopsis
pw_A_atpA <- pairwiseAlignment(seq_A, at_atpA_seq, type = "global")
message("Alignment score: ", round(score(pw_A_atpA), 1), ", identity: ", round(pid(pw_A_atpA), 1), "%")

# Build position mapping
aln_pattern_A <- as.character(pattern(pw_A_atpA))
aln_subject_A <- as.character(subject(pw_A_atpA))

pdb_pos <- 0L
at_pos <- 0L
pdb_to_at_atpA <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_A)) {
  p_char <- substr(aln_pattern_A, i, i)
  s_char <- substr(aln_subject_A, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpA <- rbind(pdb_to_at_atpA,
                            data.table(pdb_resno = coords_atpA$pdb_resno[pdb_pos],
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpA), " positions between PDB and Arabidopsis")

# Map to alignment positions
at_atpA_aln_seq <- as.character(atpA_aln[[at_idx_atpA]])
aln_chars_atpA <- strsplit(at_atpA_aln_seq, "")[[1]]

at_to_aln_atpA <- data.table(
  aln_pos = seq_along(aln_chars_atpA),
  aln_char = aln_chars_atpA,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpA)) {
  if (aln_chars_atpA[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpA$at_ungapped[i] <- ungapped
  }
}

pdb_to_aln_atpA <- merge(pdb_to_at_atpA, at_to_aln_atpA[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
pdb_to_aln_atpA <- merge(pdb_to_aln_atpA, coords_atpA[, .(pdb_resno, x, y, z)], by = "pdb_resno")

message("Mapped ", nrow(pdb_to_aln_atpA), " PDB positions to alignment positions")

# ---- Get atpA GWAS sites ----
atpA_sites <- sites_df[Gene == "atpA"]
atpA_sites <- merge(atpA_sites, pdb_to_aln_atpA[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpA GWAS sites: ", nrow(atpA_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpA_sites$x)))

# Classify sites
atpA_sites[, site_class := "not_sig"]
atpA_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpA_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpA_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpA_sites$site_class))

atpA_3d <- atpA_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpA_3d$site_class))

# ---- Calculate structural context for atpA (chain A) ----

chain_A_df <- as.data.table(chain_A_ca$atom[, c("resno", "x", "y", "z")])
setnames(chain_A_df, c("pdb_resno", "x", "y", "z"))

# Distance to nucleotide (α subunits have ADP bound)
nuc_coords <- as.data.table(pdb_atp$atom[pdb_atp$atom$resid %in% c("ATP", "ADP"), ])[, 
                                                                                     .(x = mean(x), y = mean(y), z = mean(z)), by = .(chain, resid)]

chain_A_df[, dist_to_any_nucleotide := sapply(1:.N, function(i) {
  min(sqrt((x[i] - nuc_coords$x)^2 + (y[i] - nuc_coords$y)^2 + (z[i] - nuc_coords$z)^2))
})]

# Distance to β subunits (α-β interface)
beta_coords <- ca_df_atp[subunit == "beta"]
chain_A_df[, dist_to_beta := sapply(1:.N, function(i) {
  min(sqrt((x[i] - beta_coords$x)^2 + (y[i] - beta_coords$y)^2 + (z[i] - beta_coords$z)^2))
})]

# Distance to γ (central stalk)
gamma_coords <- ca_df_atp[subunit == "gamma"]
chain_A_df[, dist_to_gamma := sapply(1:.N, function(i) {
  min(sqrt((x[i] - gamma_coords$x)^2 + (y[i] - gamma_coords$y)^2 + (z[i] - gamma_coords$z)^2))
})]

# Distance to other α subunits
other_alpha <- ca_df_atp[subunit == "alpha" & chain != "A"]
chain_A_df[, dist_to_other_alpha := sapply(1:.N, function(i) {
  min(sqrt((x[i] - other_alpha$x)^2 + (y[i] - other_alpha$y)^2 + (z[i] - other_alpha$z)^2))
})]

# Classify context
chain_A_df[, context := "bulk"]
chain_A_df[dist_to_any_nucleotide < 15, context := "near_nucleotide"]
chain_A_df[dist_to_beta < 8, context := "beta_interface"]
chain_A_df[dist_to_gamma < 8, context := "gamma_interface"]
chain_A_df[dist_to_any_nucleotide < 15 & dist_to_beta < 8, context := "catalytic_interface"]

message("\n=== atpA structural context ===")
print(table(chain_A_df$context))

# ---- Merge with GWAS sites ----
atpA_3d <- merge(atpA_3d,
                 chain_A_df[, .(pdb_resno, context, dist_to_any_nucleotide,
                                dist_to_beta, dist_to_gamma, dist_to_other_alpha)],
                 by = "pdb_resno", all.x = TRUE)

message("\n=== Site class by structural context ===")
print(table(atpA_3d$site_class, atpA_3d$context, useNA = "ifany"))

# ---- Wilcoxon tests ----
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpA_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_any_nucleotide", "dist_to_beta", "dist_to_gamma")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpA_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Correlations ----
atpA_3d[, neglog_p_only := -log10(P_aa_only)]
atpA_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Spearman correlations ===")
cor_atpA <- data.table(
  comparison = c("nucleotide vs P_only", "nucleotide vs P_pcs",
                 "beta_interface vs P_only", "beta_interface vs P_pcs",
                 "gamma_interface vs P_only", "gamma_interface vs P_pcs"),
  rho = c(
    cor(atpA_3d$dist_to_any_nucleotide, atpA_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpA_3d$dist_to_any_nucleotide, atpA_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpA_3d$dist_to_beta, atpA_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpA_3d$dist_to_beta, atpA_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpA_3d$dist_to_gamma, atpA_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpA_3d$dist_to_gamma, atpA_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpA_3d$dist_to_any_nucleotide, atpA_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpA_3d$dist_to_any_nucleotide, atpA_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpA_3d$dist_to_beta, atpA_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpA_3d$dist_to_beta, atpA_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpA_3d$dist_to_gamma, atpA_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpA_3d$dist_to_gamma, atpA_3d$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_atpA)

# ---- Secondary structure ----
atpA_struct <- struct_var[Gene == "atpA"]
atpA_full <- merge(atpA_3d, atpA_struct[, .(Position, consensus, consensus_freq, prop_C)],
                   by = "Position", all.x = TRUE)

message("\n=== Secondary structure by site class ===")
print(table(atpA_full$consensus, atpA_full$site_class, useNA = "ifany"))

ss_by_class_atpA <- atpA_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpA[, total := sum(N), by = site_class]
ss_by_class_atpA[, prop := N / total]
print(dcast(ss_by_class_atpA, site_class ~ consensus, value.var = "prop"))

message("\n=== Coil (C) enrichment in atpA ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpA_full$site_class == sc & atpA_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpA_full$site_class == sc & atpA_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpA_full$site_class == "not_sig" & atpA_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpA_full$site_class == "not_sig" & atpA_full$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Quick visualization ----
chain_A_df <- merge(chain_A_df, pdb_to_aln_atpA[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_A_df <- merge(chain_A_df, atpA_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_A_df[is.na(site_class), site_class := "no_data"]

p_atpA_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "A"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_A_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_A_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_A_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "atpA (α subunit) GWAS hits",
       subtitle = "Gold = sig_no_control, Blue = sig_with_control, Red = sig_both",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_atpA_xy

# ---- atpA XZ and YZ projections ----

p_atpA_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "A"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_A_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_A_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_A_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_z), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "atpA - XZ projection (side view)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atpA_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "A"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_A_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_A_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_A_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_A_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_y, y = lig_z), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "atpA - YZ projection (side view 2)",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_atpA_xy / (p_atpA_xz + p_atpA_yz)

# ---- atpF (b subunit - peripheral stalk) ----

# atpF is the "b" subunit - part of the peripheral stalk
# Check which chain is atpF in the structure
# Chain "b" (lowercase) should be it based on naming

chain_b_test <- trim.pdb(pdb_atp, chain = "b")
ca_b <- atom.select(chain_b_test, elety = "CA")
chain_b_ca <- trim.pdb(chain_b_test, ca_b)
message("Chain 'b' has ", nrow(chain_b_ca$atom), " residues")

# Get sequence
coords_atpF <- as.data.table(chain_b_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpF, c("pdb_resno", "aa", "x", "y", "z"))
seq_b <- paste(aa_map[coords_atpF$aa], collapse = "")
message("Chain b sequence length: ", nchar(seq_b))

# Load atpF alignment
atpF_aln <- readAAStringSet("data/tmp/alignedGenes/atpF_AA_aligned.fasta")
message("atpF alignment: ", length(atpF_aln), " sequences, length ", width(atpF_aln)[1])

at_idx_atpF <- grep(arabidopsis_th_id, names(atpF_aln))
stopifnot(length(at_idx_atpF) == 1)
at_atpF_seq <- gsub("-", "", as.character(atpF_aln[[at_idx_atpF]]))
message("Arabidopsis atpF length: ", nchar(at_atpF_seq))

# Align PDB to Arabidopsis
pw_atpF <- pairwiseAlignment(seq_b, at_atpF_seq, type = "local")
message("Alignment score: ", round(score(pw_atpF), 1), ", identity: ", round(pid(pw_atpF), 1), "%")
print(pw_atpF)

# ---- atpF (b subunit - peripheral stalk) continued ----

# Build position mapping
aln_pattern_F <- as.character(pattern(pw_atpF))
aln_subject_F <- as.character(subject(pw_atpF))

pdb_pos <- 0L
at_pos <- start(subject(pw_atpF)) - 1L  # Start at the alignment start position (23-1=22)
pdb_to_at_atpF <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_F)) {
  p_char <- substr(aln_pattern_F, i, i)
  s_char <- substr(aln_subject_F, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpF <- rbind(pdb_to_at_atpF,
                            data.table(pdb_resno = coords_atpF$pdb_resno[pdb_pos],
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpF), " positions between PDB and Arabidopsis")

# Map to alignment positions
at_atpF_aln_seq <- as.character(atpF_aln[[at_idx_atpF]])
aln_chars_atpF <- strsplit(at_atpF_aln_seq, "")[[1]]

at_to_aln_atpF <- data.table(
  aln_pos = seq_along(aln_chars_atpF),
  aln_char = aln_chars_atpF,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpF)) {
  if (aln_chars_atpF[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpF$at_ungapped[i] <- ungapped
  }
}

pdb_to_aln_atpF <- merge(pdb_to_at_atpF, at_to_aln_atpF[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
pdb_to_aln_atpF <- merge(pdb_to_aln_atpF, coords_atpF[, .(pdb_resno, x, y, z)], by = "pdb_resno")

message("Mapped ", nrow(pdb_to_aln_atpF), " PDB positions to alignment positions")

# ---- Get atpF GWAS sites ----
atpF_sites <- sites_df[Gene == "atpF"]
atpF_sites <- merge(atpF_sites, pdb_to_aln_atpF[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpF GWAS sites: ", nrow(atpF_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpF_sites$x)))

# Classify sites
atpF_sites[, site_class := "not_sig"]
atpF_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpF_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpF_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpF_sites$site_class))

atpF_3d <- atpF_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpF_3d$site_class))

# ---- Calculate structural context for atpF ----
# atpF (b subunit) forms the peripheral stalk connecting F1 to membrane
# Key interfaces: with F1 head (α/β), with a-subunit, with OSCP

chain_b_df <- as.data.table(chain_b_ca$atom[, c("resno", "x", "y", "z")])
setnames(chain_b_df, c("pdb_resno", "x", "y", "z"))

# Distance to F1 head (α and β subunits)
f1_coords <- ca_df_atp[subunit %in% c("alpha", "beta")]
chain_b_df[, dist_to_F1 := sapply(1:.N, function(i) {
  min(sqrt((x[i] - f1_coords$x)^2 + (y[i] - f1_coords$y)^2 + (z[i] - f1_coords$z)^2))
})]

# Distance to a-subunit (membrane interface)
a_coords <- ca_df_atp[subunit == "a"]
chain_b_df[, dist_to_a := sapply(1:.N, function(i) {
  min(sqrt((x[i] - a_coords$x)^2 + (y[i] - a_coords$y)^2 + (z[i] - a_coords$z)^2))
})]

# Distance to c-ring
c_coords <- ca_df_atp[subunit == "c-ring"]
chain_b_df[, dist_to_c_ring := sapply(1:.N, function(i) {
  min(sqrt((x[i] - c_coords$x)^2 + (y[i] - c_coords$y)^2 + (z[i] - c_coords$z)^2))
})]

# Distance to delta/OSCP (peripheral stalk top)
delta_coords <- ca_df_atp[subunit %in% c("delta", "OSCP")]
chain_b_df[, dist_to_delta := sapply(1:.N, function(i) {
  min(sqrt((x[i] - delta_coords$x)^2 + (y[i] - delta_coords$y)^2 + (z[i] - delta_coords$z)^2))
})]

# Classify context
chain_b_df[, context := "bulk"]
chain_b_df[dist_to_F1 < 8, context := "F1_interface"]
chain_b_df[dist_to_a < 8, context := "a_interface"]
chain_b_df[dist_to_delta < 8, context := "delta_interface"]

message("\n=== atpF structural context ===")
print(table(chain_b_df$context))

# ---- Merge with GWAS sites ----
atpF_3d <- merge(atpF_3d,
                 chain_b_df[, .(pdb_resno, context, dist_to_F1, dist_to_a, dist_to_c_ring, dist_to_delta)],
                 by = "pdb_resno", all.x = TRUE)

message("\n=== Site class by structural context ===")
print(table(atpF_3d$site_class, atpF_3d$context, useNA = "ifany"))

# ---- Wilcoxon tests ----
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpF_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_F1", "dist_to_a", "dist_to_delta")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpF_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Correlations ----
atpF_3d[, neglog_p_only := -log10(P_aa_only)]
atpF_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Spearman correlations ===")
cor_atpF <- data.table(
  comparison = c("F1 vs P_only", "F1 vs P_pcs",
                 "a_subunit vs P_only", "a_subunit vs P_pcs",
                 "delta vs P_only", "delta vs P_pcs"),
  rho = c(
    cor(atpF_3d$dist_to_F1, atpF_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpF_3d$dist_to_F1, atpF_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpF_3d$dist_to_a, atpF_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpF_3d$dist_to_a, atpF_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpF_3d$dist_to_delta, atpF_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpF_3d$dist_to_delta, atpF_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpF_3d$dist_to_F1, atpF_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpF_3d$dist_to_F1, atpF_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpF_3d$dist_to_a, atpF_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpF_3d$dist_to_a, atpF_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpF_3d$dist_to_delta, atpF_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpF_3d$dist_to_delta, atpF_3d$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_atpF)

# ---- Secondary structure ----
atpF_struct <- struct_var[Gene == "atpF"]
atpF_full <- merge(atpF_3d, atpF_struct[, .(Position, consensus, consensus_freq, prop_C)],
                   by = "Position", all.x = TRUE)

message("\n=== Secondary structure by site class ===")
print(table(atpF_full$consensus, atpF_full$site_class, useNA = "ifany"))

ss_by_class_atpF <- atpF_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpF[, total := sum(N), by = site_class]
ss_by_class_atpF[, prop := N / total]
print(dcast(ss_by_class_atpF, site_class ~ consensus, value.var = "prop"))

message("\n=== Coil (C) enrichment in atpF ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpF_full$site_class == sc & atpF_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpF_full$site_class == sc & atpF_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpF_full$site_class == "not_sig" & atpF_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpF_full$site_class == "not_sig" & atpF_full$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualizations ----
chain_b_df <- merge(chain_b_df, pdb_to_aln_atpF[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_b_df <- merge(chain_b_df, atpF_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_b_df[is.na(site_class), site_class := "no_data"]

# Three projections
p_atpF_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "b"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_b_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_b_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_b_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpF (b subunit, peripheral stalk) - XY",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_atpF_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "b"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_b_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_b_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_b_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpF - XZ (side view)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atpF_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "b"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_b_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_b_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_b_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_b_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpF - YZ (side view 2)",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_atpF_xy / (p_atpF_xz + p_atpF_yz)

# ---- atpE (ε subunit) ----

# atpE is the epsilon subunit - connects γ to the c-ring
# Check which chain - likely lowercase 'e'

chain_e_test <- trim.pdb(pdb_atp, chain = "e")
ca_e <- atom.select(chain_e_test, elety = "CA")
chain_e_ca <- trim.pdb(chain_e_test, ca_e)
message("Chain 'e' has ", nrow(chain_e_ca$atom), " residues")

coords_atpE <- as.data.table(chain_e_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpE, c("pdb_resno", "aa", "x", "y", "z"))
seq_e <- paste(aa_map[coords_atpE$aa], collapse = "")
message("Chain e sequence length: ", nchar(seq_e))

# Load atpE alignment
atpE_aln <- readAAStringSet("data/tmp/alignedGenes/atpE_AA_aligned.fasta")
message("atpE alignment: ", length(atpE_aln), " sequences, length ", width(atpE_aln)[1])

at_idx_atpE <- grep(arabidopsis_th_id, names(atpE_aln))
stopifnot(length(at_idx_atpE) == 1)
at_atpE_seq <- gsub("-", "", as.character(atpE_aln[[at_idx_atpE]]))
message("Arabidopsis atpE length: ", nchar(at_atpE_seq))

# Align PDB to Arabidopsis
pw_atpE <- pairwiseAlignment(seq_e, at_atpE_seq, type = "local")
message("Alignment score: ", round(score(pw_atpE), 1), ", identity: ", round(pid(pw_atpE), 1), "%")
print(pw_atpE)

# ---- atpE (ε subunit) continued ----

# Build position mapping
aln_pattern_E <- as.character(pattern(pw_atpE))
aln_subject_E <- as.character(subject(pw_atpE))

pdb_pos <- 0L
at_pos <- start(subject(pw_atpE)) - 1L  # starts at 1, so 0
pdb_to_at_atpE <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_E)) {
  p_char <- substr(aln_pattern_E, i, i)
  s_char <- substr(aln_subject_E, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpE <- rbind(pdb_to_at_atpE,
                            data.table(pdb_resno = coords_atpE$pdb_resno[pdb_pos],
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpE), " positions between PDB and Arabidopsis")

# Map to alignment positions
at_atpE_aln_seq <- as.character(atpE_aln[[at_idx_atpE]])
aln_chars_atpE <- strsplit(at_atpE_aln_seq, "")[[1]]

at_to_aln_atpE <- data.table(
  aln_pos = seq_along(aln_chars_atpE),
  aln_char = aln_chars_atpE,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpE)) {
  if (aln_chars_atpE[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpE$at_ungapped[i] <- ungapped
  }
}

pdb_to_aln_atpE <- merge(pdb_to_at_atpE, at_to_aln_atpE[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
pdb_to_aln_atpE <- merge(pdb_to_aln_atpE, coords_atpE[, .(pdb_resno, x, y, z)], by = "pdb_resno")

message("Mapped ", nrow(pdb_to_aln_atpE), " PDB positions to alignment positions")

# ---- Get atpE GWAS sites ----
atpE_sites <- sites_df[Gene == "atpE"]
atpE_sites <- merge(atpE_sites, pdb_to_aln_atpE[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpE GWAS sites: ", nrow(atpE_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpE_sites$x)))

# Classify sites
atpE_sites[, site_class := "not_sig"]
atpE_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpE_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpE_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpE_sites$site_class))

atpE_3d <- atpE_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpE_3d$site_class))

# ---- Calculate structural context for atpE ----
# atpE (ε) connects γ to c-ring, acts as a regulatory/coupling subunit

chain_e_df <- as.data.table(chain_e_ca$atom[, c("resno", "x", "y", "z")])
setnames(chain_e_df, c("pdb_resno", "x", "y", "z"))

# Distance to γ (central stalk)
gamma_coords <- ca_df_atp[subunit == "gamma"]
chain_e_df[, dist_to_gamma := sapply(1:.N, function(i) {
  min(sqrt((x[i] - gamma_coords$x)^2 + (y[i] - gamma_coords$y)^2 + (z[i] - gamma_coords$z)^2))
})]

# Distance to c-ring
c_coords <- ca_df_atp[subunit == "c-ring"]
chain_e_df[, dist_to_c_ring := sapply(1:.N, function(i) {
  min(sqrt((x[i] - c_coords$x)^2 + (y[i] - c_coords$y)^2 + (z[i] - c_coords$z)^2))
})]

# Distance to β subunits (ε can interact with β in some conformations)
beta_coords <- ca_df_atp[subunit == "beta"]
chain_e_df[, dist_to_beta := sapply(1:.N, function(i) {
  min(sqrt((x[i] - beta_coords$x)^2 + (y[i] - beta_coords$y)^2 + (z[i] - beta_coords$z)^2))
})]

# Distance to α subunits
alpha_coords <- ca_df_atp[subunit == "alpha"]
chain_e_df[, dist_to_alpha := sapply(1:.N, function(i) {
  min(sqrt((x[i] - alpha_coords$x)^2 + (y[i] - alpha_coords$y)^2 + (z[i] - alpha_coords$z)^2))
})]

# Classify context
chain_e_df[, context := "bulk"]
chain_e_df[dist_to_gamma < 8, context := "gamma_interface"]
chain_e_df[dist_to_c_ring < 8, context := "c_ring_interface"]
chain_e_df[dist_to_beta < 8, context := "beta_interface"]

message("\n=== atpE structural context ===")
print(table(chain_e_df$context))

# ---- Merge with GWAS sites ----
atpE_3d <- merge(atpE_3d,
                 chain_e_df[, .(pdb_resno, context, dist_to_gamma, dist_to_c_ring, dist_to_beta, dist_to_alpha)],
                 by = "pdb_resno", all.x = TRUE)

message("\n=== Site class by structural context ===")
print(table(atpE_3d$site_class, atpE_3d$context, useNA = "ifany"))

# ---- Wilcoxon tests ----
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpE_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_gamma", "dist_to_c_ring", "dist_to_beta")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpE_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Correlations ----
atpE_3d[, neglog_p_only := -log10(P_aa_only)]
atpE_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Spearman correlations ===")
cor_atpE <- data.table(
  comparison = c("gamma vs P_only", "gamma vs P_pcs",
                 "c_ring vs P_only", "c_ring vs P_pcs",
                 "beta vs P_only", "beta vs P_pcs"),
  rho = c(
    cor(atpE_3d$dist_to_gamma, atpE_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpE_3d$dist_to_gamma, atpE_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpE_3d$dist_to_c_ring, atpE_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpE_3d$dist_to_c_ring, atpE_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpE_3d$dist_to_beta, atpE_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpE_3d$dist_to_beta, atpE_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpE_3d$dist_to_gamma, atpE_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpE_3d$dist_to_gamma, atpE_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpE_3d$dist_to_c_ring, atpE_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpE_3d$dist_to_c_ring, atpE_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpE_3d$dist_to_beta, atpE_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpE_3d$dist_to_beta, atpE_3d$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_atpE)

# ---- Secondary structure ----
atpE_struct <- struct_var[Gene == "atpE"]
atpE_full <- merge(atpE_3d, atpE_struct[, .(Position, consensus, consensus_freq, prop_C)],
                   by = "Position", all.x = TRUE)

message("\n=== Secondary structure by site class ===")
print(table(atpE_full$consensus, atpE_full$site_class, useNA = "ifany"))

ss_by_class_atpE <- atpE_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpE[, total := sum(N), by = site_class]
ss_by_class_atpE[, prop := N / total]
print(dcast(ss_by_class_atpE, site_class ~ consensus, value.var = "prop"))

message("\n=== Coil (C) enrichment in atpE ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpE_full$site_class == sc & atpE_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpE_full$site_class == sc & atpE_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpE_full$site_class == "not_sig" & atpE_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpE_full$site_class == "not_sig" & atpE_full$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualizations ----
chain_e_df <- merge(chain_e_df, pdb_to_aln_atpE[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_e_df <- merge(chain_e_df, atpE_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_e_df[is.na(site_class), site_class := "no_data"]

# Three projections
p_atpE_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "e"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_e_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_e_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_e_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpE (ε subunit) - XY",
       subtitle = "Connects γ to c-ring",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_atpE_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "e"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_e_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_e_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_e_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpE - XZ (side view)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atpE_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "e"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_e_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_e_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_e_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_e_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpE - YZ (side view 2)",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_atpE_xy / (p_atpE_xz + p_atpE_yz)

# ---- atpH ----
# Check which chain - likely lowercase 'e'

chain_h_test <- trim.pdb(pdb_atp, chain = "T")
ca_h <- atom.select(chain_h_test, elety = "CA")
chain_h_ca <- trim.pdb(chain_h_test, ca_h)
message("Chain 'T' has ", nrow(chain_h_ca$atom), " residues")

coords_atpH <- as.data.table(chain_h_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpH, c("pdb_resno", "aa", "x", "y", "z"))
seq_h <- paste(aa_map[coords_atpH$aa], collapse = "")
message("Chain h sequence length: ", nchar(seq_h))

# Load atpE alignment
atpH_aln <- readAAStringSet("data/tmp/alignedGenes/atpH_AA_aligned.fasta")
message("atpH alignment: ", length(atpH_aln), " sequences, length ", width(atpH_aln)[1])

at_idx_atpH <- grep(arabidopsis_th_id, names(atpH_aln))
stopifnot(length(at_idx_atpH) == 1)
at_atpH_seq <- gsub("-", "", as.character(atpH_aln[[at_idx_atpH]]))
message("Arabidopsis atpH length: ", nchar(at_atpH_seq))

# Align PDB to Arabidopsis
pw_atpH <- pairwiseAlignment(seq_h, at_atpH_seq, type = "local")
message("Alignment score: ", round(score(pw_atpH), 1), ", identity: ", round(pid(pw_atpH), 1), "%")
print(pw_atpH)

# ---- atpE (ε subunit) continued ----

# Build position mapping
aln_pattern_H <- as.character(pattern(pw_atpH))
aln_subject_H <- as.character(subject(pw_atpH))

pdb_pos <- 0L
at_pos <- start(subject(pw_atpH)) - 1L  # starts at 1, so 0
pdb_to_at_atpH <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_H)) {
  p_char <- substr(aln_pattern_H, i, i)
  s_char <- substr(aln_subject_H, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpH <- rbind(pdb_to_at_atpH,
                            data.table(pdb_resno = coords_atpH$pdb_resno[pdb_pos],
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpH), " positions between PDB and Arabidopsis")

# Map to alignment positions
at_atpH_aln_seq <- as.character(atpH_aln[[at_idx_atpH]])
aln_chars_atpH <- strsplit(at_atpH_aln_seq, "")[[1]]

at_to_aln_atpH <- data.table(
  aln_pos = seq_along(aln_chars_atpH),
  aln_char = aln_chars_atpH,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpH)) {
  if (aln_chars_atpH[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpH$at_ungapped[i] <- ungapped
  }
}

pdb_to_aln_atpH <- merge(pdb_to_at_atpH, at_to_aln_atpH[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
pdb_to_aln_atpH <- merge(pdb_to_aln_atpH, coords_atpH[, .(pdb_resno, x, y, z)], by = "pdb_resno")

message("Mapped ", nrow(pdb_to_aln_atpH), " PDB positions to alignment positions")

# ---- Get atpE GWAS sites ----
atpH_sites <- sites_df[Gene == "atpH"]
atpH_sites <- merge(atpH_sites, pdb_to_aln_atpH[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpH GWAS sites: ", nrow(atpH_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpH_sites$x)))

# Classify sites
atpH_sites[, site_class := "not_sig"]
atpH_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpH_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpH_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpH_sites$site_class))

atpH_3d <- atpH_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpH_3d$site_class))

# ---- Calculate structural context for atpE ----
# atpE (ε) connects γ to c-ring, acts as a regulatory/coupling subunit

chain_h_df <- as.data.table(chain_h_ca$atom[, c("resno", "x", "y", "z")])
setnames(chain_h_df, c("pdb_resno", "x", "y", "z"))

# Distance to γ (central stalk)
gamma_coords <- ca_df_atp[subunit == "gamma"]
chain_h_df[, dist_to_gamma := sapply(1:.N, function(i) {
  min(sqrt((x[i] - gamma_coords$x)^2 + (y[i] - gamma_coords$y)^2 + (z[i] - gamma_coords$z)^2))
})]
#instead of doing this gamma thing, i want difference from membrane,
#approx by distance from stroma (Z=150) and lumen (z=100)
chain_h_df[, dist_to_stroma := sapply(1:.N, function(i) {
  abs(150 - z[i])
})]
plot(chain_h_df$dist_to_stroma)
chain_h_df[, dist_to_lumen := sapply(1:.N, function(i) {
  abs(100 - z[i])
})]
plot(chain_h_df$dist_to_lumen)
# Distance to c-ring

# Classify context
chain_h_df[, context := "bulk"]
chain_h_df[dist_to_lumen < 5, context := "lumen_interface"]
chain_h_df[dist_to_stroma < 5, context := "stroma_interface"]

message("\n=== atpH structural context ===")
print(table(chain_h_df$context))

# ---- Merge with GWAS sites ----
atpH_3d <- merge(atpH_3d,
                 chain_h_df[, .(pdb_resno, context, dist_to_stroma, dist_to_lumen)],
                 by = "pdb_resno", all.x = TRUE)

message("\n=== Site class by structural context ===")
print(table(atpH_3d$site_class, atpH_3d$context, useNA = "ifany"))

# ---- Wilcoxon tests ----
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpH_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_stroma", "dist_to_lumen")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpH_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Correlations ----
atpH_3d[, neglog_p_only := -log10(P_aa_only)]
atpH_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Spearman correlations ===")
cor_atpH <- data.table(
  comparison = c("gamma vs P_only", "gamma vs P_pcs",
                 "c_ring vs P_only", "c_ring vs P_pcs",
                 "beta vs P_only", "beta vs P_pcs"),
  rho = c(
    cor(atpH_3d$dist_to_stroma, atpH_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpH_3d$dist_to_stroma, atpH_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpH_3d$dist_to_lumen, atpH_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpH_3d$dist_to_lumen, atpH_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpH_3d$dist_to_stroma, atpH_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpH_3d$dist_to_stroma, atpH_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpH_3d$dist_to_lumen, atpH_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpH_3d$dist_to_lumen, atpH_3d$neglog_p_pcs, method = "spearman")$p.value)
)
print(cor_atpH)

# ---- Secondary structure ----
atpH_struct <- struct_var[Gene == "atpH"]
atpH_full <- merge(atpH_3d, atpE_struct[, .(Position, consensus, consensus_freq, prop_C)],
                   by = "Position", all.x = TRUE)

message("\n=== Secondary structure by site class ===")
print(table(atpH_full$consensus, atpH_full$site_class, useNA = "ifany"))

ss_by_class_atpH <- atpH_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpH[, total := sum(N), by = site_class]
ss_by_class_atpH[, prop := N / total]
print(dcast(ss_by_class_atpH, site_class ~ consensus, value.var = "prop"))

message("\n=== Coil (C) enrichment in atpE ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpH_full$site_class == sc & atpH_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpH_full$site_class == sc & atpH_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpH_full$site_class == "not_sig" & atpH_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpH_full$site_class == "not_sig" & atpH_full$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualizations ----
chain_h_df <- merge(chain_h_df, pdb_to_aln_atpH[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_h_df <- merge(chain_h_df, atpH_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_h_df[is.na(site_class), site_class := "no_data"]

# Three projections
p_atpH_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "T"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_h_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_h_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_h_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpH (c subunit) - XY",
       subtitle = "c-ring",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()
p_atpH_xy
p_atpH_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "T"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_h_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_h_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_h_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpH - XZ (side view)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atpH_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "T"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_h_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_h_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_h_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_h_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpH - YZ (side view 2)",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_atpH_xy / (p_atpH_xz + p_atpH_yz)

# ---- atpI (a subunit - proton channel) ----

# atpI is the "a" subunit - membrane embedded, forms proton channel with c-ring
# Chain "a" (lowercase) should be it

chain_a_test <- trim.pdb(pdb_atp, chain = "a")
ca_a <- atom.select(chain_a_test, elety = "CA")
chain_a_ca <- trim.pdb(chain_a_test, ca_a)
message("Chain 'a' has ", nrow(chain_a_ca$atom), " residues")

coords_atpI <- as.data.table(chain_a_ca$atom[, c("resno", "resid", "x", "y", "z")])
setnames(coords_atpI, c("pdb_resno", "aa", "x", "y", "z"))
seq_a <- paste(aa_map[coords_atpI$aa], collapse = "")
message("Chain a sequence length: ", nchar(seq_a))

# Load atpI alignment
atpI_aln <- readAAStringSet("data/tmp/alignedGenes/atpI_AA_aligned.fasta")
message("atpI alignment: ", length(atpI_aln), " sequences, length ", width(atpI_aln)[1])

at_idx_atpI <- grep(arabidopsis_th_id, names(atpI_aln))
stopifnot(length(at_idx_atpI) == 1)
at_atpI_seq <- gsub("-", "", as.character(atpI_aln[[at_idx_atpI]]))
message("Arabidopsis atpI length: ", nchar(at_atpI_seq))

# Align PDB to Arabidopsis
pw_atpI <- pairwiseAlignment(seq_a, at_atpI_seq, type = "local")
message("Alignment score: ", round(score(pw_atpI), 1), ", identity: ", round(pid(pw_atpI), 1), "%")
print(pw_atpI)

# Build position mapping
aln_pattern_I <- as.character(pattern(pw_atpI))
aln_subject_I <- as.character(subject(pw_atpI))

pdb_pos <- 0L
at_pos <- start(subject(pw_atpI)) - 1L
pdb_to_at_atpI <- data.table(pdb_resno = integer(), at_ungapped = integer())

for (i in 1:nchar(aln_pattern_I)) {
  p_char <- substr(aln_pattern_I, i, i)
  s_char <- substr(aln_subject_I, i, i)
  
  if (p_char != "-") pdb_pos <- pdb_pos + 1L
  if (s_char != "-") at_pos <- at_pos + 1L
  
  if (p_char != "-" & s_char != "-") {
    pdb_to_at_atpI <- rbind(pdb_to_at_atpI,
                            data.table(pdb_resno = coords_atpI$pdb_resno[pdb_pos],
                                       at_ungapped = at_pos))
  }
}

message("Mapped ", nrow(pdb_to_at_atpI), " positions between PDB and Arabidopsis")

# Map to alignment positions
at_atpI_aln_seq <- as.character(atpI_aln[[at_idx_atpI]])
aln_chars_atpI <- strsplit(at_atpI_aln_seq, "")[[1]]

at_to_aln_atpI <- data.table(
  aln_pos = seq_along(aln_chars_atpI),
  aln_char = aln_chars_atpI,
  at_ungapped = NA_integer_
)

ungapped <- 0L
for (i in seq_along(aln_chars_atpI)) {
  if (aln_chars_atpI[i] != "-") {
    ungapped <- ungapped + 1L
    at_to_aln_atpI$at_ungapped[i] <- ungapped
  }
}

pdb_to_aln_atpI <- merge(pdb_to_at_atpI, at_to_aln_atpI[!is.na(at_ungapped), .(aln_pos, at_ungapped)],
                         by = "at_ungapped")
pdb_to_aln_atpI <- merge(pdb_to_aln_atpI, coords_atpI[, .(pdb_resno, x, y, z)], by = "pdb_resno")

message("Mapped ", nrow(pdb_to_aln_atpI), " PDB positions to alignment positions")

# ---- Get atpI GWAS sites ----
atpI_sites <- sites_df[Gene == "atpI"]
atpI_sites <- merge(atpI_sites, pdb_to_aln_atpI[, .(aln_pos, pdb_resno, x, y, z)],
                    by.x = "Position", by.y = "aln_pos", all.x = TRUE)

message("\natpI GWAS sites: ", nrow(atpI_sites))
message("Sites with 3D coordinates: ", sum(!is.na(atpI_sites$x)))

# Classify sites
atpI_sites[, site_class := "not_sig"]
atpI_sites[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
atpI_sites[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
atpI_sites[P_aa_with_pcs < thresh_aa_pcs & P_aa_only < thresh_aa_only, site_class := "sig_both"]

message("\nSite classification:")
print(table(atpI_sites$site_class))

atpI_3d <- atpI_sites[!is.na(x)]
message("\nSites with 3D coords by class:")
print(table(atpI_3d$site_class))

# ---- Calculate structural context for atpI ----
# atpI (a subunit) is membrane-embedded, forms proton channel with c-ring
# Key features: proton half-channels, c-ring interface, membrane boundaries

chain_a_df <- as.data.table(chain_a_ca$atom[, c("resno", "x", "y", "z")])
setnames(chain_a_df, c("pdb_resno", "x", "y", "z"))

# Distance to c-ring (critical interface for proton translocation)
c_coords <- ca_df_atp[subunit == "c-ring"]
chain_a_df[, dist_to_c_ring := sapply(1:.N, function(i) {
  min(sqrt((x[i] - c_coords$x)^2 + (y[i] - c_coords$y)^2 + (z[i] - c_coords$z)^2))
})]

# Distance to b subunit (peripheral stalk anchor)
b_coords <- ca_df_atp[chain == "b"]
chain_a_df[, dist_to_b := sapply(1:.N, function(i) {
  min(sqrt((x[i] - b_coords$x)^2 + (y[i] - b_coords$y)^2 + (z[i] - b_coords$z)^2))
})]

# Membrane position - using Z coordinate like you did for atpH
# Approximate: stroma ~150, lumen ~100
chain_a_df[, dist_to_stroma := abs(150 - z)]
chain_a_df[, dist_to_lumen := abs(100 - z)]

# Check Z range
message("\natpI Z range: ", round(min(chain_a_df$z), 1), " - ", round(max(chain_a_df$z), 1))

# Classify context
chain_a_df[, context := "membrane"]
chain_a_df[dist_to_c_ring < 8, context := "c_ring_interface"]
chain_a_df[dist_to_stroma < 10, context := "stroma_exposed"]
chain_a_df[dist_to_lumen < 10, context := "lumen_exposed"]

message("\n=== atpI structural context ===")
print(table(chain_a_df$context))

# ---- Merge with GWAS sites ----
atpI_3d <- merge(atpI_3d,
                 chain_a_df[, .(pdb_resno, context, dist_to_c_ring, dist_to_b, 
                                dist_to_stroma, dist_to_lumen)],
                 by = "pdb_resno", all.x = TRUE)

message("\n=== Site class by structural context ===")
print(table(atpI_3d$site_class, atpI_3d$context, useNA = "ifany"))

# ---- Wilcoxon tests ----
message("\n=== Wilcoxon tests vs not_sig ===")
bg <- atpI_3d[site_class == "not_sig"]

for (dist_var in c("dist_to_c_ring", "dist_to_stroma", "dist_to_lumen")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atpI_3d[site_class == sc, get(dist_var)]
    bg_vals <- bg[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 3) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑ farther", "↓ closer")
    message(sprintf("  %s: median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Correlations ----
atpI_3d[, neglog_p_only := -log10(P_aa_only)]
atpI_3d[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Spearman correlations ===")
cor_atpI <- data.table(
  comparison = c("c_ring vs P_only", "c_ring vs P_pcs",
                 "stroma vs P_only", "stroma vs P_pcs",
                 "lumen vs P_only", "lumen vs P_pcs"),
  rho = c(
    cor(atpI_3d$dist_to_c_ring, atpI_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpI_3d$dist_to_c_ring, atpI_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpI_3d$dist_to_stroma, atpI_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpI_3d$dist_to_stroma, atpI_3d$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atpI_3d$dist_to_lumen, atpI_3d$neglog_p_only, method = "spearman", use = "complete"),
    cor(atpI_3d$dist_to_lumen, atpI_3d$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atpI_3d$dist_to_c_ring, atpI_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpI_3d$dist_to_c_ring, atpI_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpI_3d$dist_to_stroma, atpI_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpI_3d$dist_to_stroma, atpI_3d$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atpI_3d$dist_to_lumen, atpI_3d$neglog_p_only, method = "spearman")$p.value,
    cor.test(atpI_3d$dist_to_lumen, atpI_3d$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_atpI)

# ---- Secondary structure ----
atpI_struct <- struct_var[Gene == "atpI"]
atpI_full <- merge(atpI_3d, atpI_struct[, .(Position, consensus, consensus_freq, prop_C)],
                   by = "Position", all.x = TRUE)

message("\n=== Secondary structure by site class ===")
print(table(atpI_full$consensus, atpI_full$site_class, useNA = "ifany"))

ss_by_class_atpI <- atpI_full[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_by_class_atpI[, total := sum(N), by = site_class]
ss_by_class_atpI[, prop := N / total]
print(dcast(ss_by_class_atpI, site_class ~ consensus, value.var = "prop"))

message("\n=== Coil (C) enrichment in atpI ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atpI_full$site_class == sc & atpI_full$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atpI_full$site_class == sc & atpI_full$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atpI_full$site_class == "not_sig" & atpI_full$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atpI_full$site_class == "not_sig" & atpI_full$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / max(sig_C + sig_notC, 1)
  prop_bg <- bg_C / max(bg_C + bg_notC, 1)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualizations ----
chain_a_df <- merge(chain_a_df, pdb_to_aln_atpI[, .(pdb_resno, aln_pos)],
                    by = "pdb_resno", all.x = TRUE)
chain_a_df <- merge(chain_a_df, atpI_sites[, .(Position, site_class)],
                    by.x = "aln_pos", by.y = "Position", all.x = TRUE)
chain_a_df[is.na(site_class), site_class := "no_data"]

# Three projections
p_atpI_xy <- ggplot() +
  geom_point(data = ca_df_atp[chain != "a"], aes(x = x, y = y),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_a_df[site_class == "no_data"], aes(x = x, y = y),
             color = "gray70", size = 1) +
  geom_point(data = chain_a_df[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_a_df[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpI (a subunit, proton channel) - XY",
       subtitle = "Membrane-embedded, interfaces with c-ring",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_atpI_xz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "a"], aes(x = x, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_a_df[site_class == "no_data"], aes(x = x, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_a_df[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_a_df[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpI - XZ (side view)",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_atpI_yz <- ggplot() +
  geom_point(data = ca_df_atp[chain != "a"], aes(x = y, y = z),
             color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = chain_a_df[site_class == "no_data"], aes(x = y, y = z),
             color = "gray70", size = 1) +
  geom_point(data = chain_a_df[site_class == "not_sig"], aes(x = y, y = z),
             color = "gray50", size = 1.5) +
  geom_point(data = chain_a_df[site_class == "sig_no_control"], aes(x = y, y = z),
             color = "gold", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_with_control"], aes(x = y, y = z),
             color = "steelblue", size = 2.5) +
  geom_point(data = chain_a_df[site_class == "sig_both"], aes(x = y, y = z),
             color = "darkred", size = 3) +
  coord_fixed() +
  labs(title = "atpI - YZ (side view 2)",
       x = "Y (Å)", y = "Z (Å)") +
  theme_classic()

p_atpI_xy / (p_atpI_xz + p_atpI_yz)

# ---- COMBINED ATP SYNTHASE ANALYSIS ----

# First, let's establish membrane boundaries from the structure
message("Z-coordinate ranges by subunit type:")
ca_df_atp[, .(min_z = min(z), max_z = max(z), mean_z = mean(z)), by = subunit][order(mean_z)]

# Define compartments based on Z (adjust based on your structure)
# Looking at typical ATP synthase orientation:
# High Z = stroma (F1 head)
# Mid Z = membrane 
# Low Z = lumen

# Let's check the actual Z distribution
hist(ca_df_atp$z, breaks = 50, main = "Z distribution of ATP synthase", xlab = "Z (Å)")
abline(v = c(100, 130), col = "red", lty = 2)  # approximate membrane boundaries

# ---- Combine all subunit data ----

# Add gene labels and combine
atpA_3d[, gene := "atpA"]
atpB_3d[, gene := "atpB"]
atpE_3d[, gene := "atpE"]
atpF_3d[, gene := "atpF"]
atpH_3d[, gene := "atpH"]
atpI_3d[, gene := "atpI"]

# Get common columns
common_cols <- c("Position", "Gene", "P_aa_only", "P_aa_with_pcs", "N", 
                 "pdb_resno", "x", "y", "z", "site_class", "gene",
                 "neglog_p_only", "neglog_p_pcs")

# Standardize column names and combine
# First ensure all have x, y, z and basic columns

atp_combined <- rbindlist(list(
  atpA_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpA")],
  atpB_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpB")],
  atpE_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpE")],
  atpF_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpF")],
  atpH_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpH")],
  atpI_3d[, .(Position, Gene, P_aa_only, P_aa_with_pcs, pdb_resno, x, y, z, site_class, gene = "atpI")]
), fill = TRUE)

atp_combined[, neglog_p_only := -log10(P_aa_only)]
atp_combined[, neglog_p_pcs := -log10(P_aa_with_pcs)]

message("\n=== Combined ATP synthase dataset ===")
message("Total sites: ", nrow(atp_combined))
print(table(atp_combined$gene, atp_combined$site_class))

# ---- Define compartments based on Z ----

# Check Z ranges
message("\nZ ranges by gene:")
print(atp_combined[, .(min_z = min(z, na.rm=T), max_z = max(z, na.rm=T), mean_z = mean(z, na.rm=T)), by = gene][order(mean_z)])

# Define compartments (adjust these based on your histogram)
# Typical: lumen < 105, membrane 105-135, stroma > 135
atp_combined[, compartment := fcase(
  z < 105, "lumen",
  z >= 105 & z < 135, "membrane",
  z >= 135, "stroma"
)]

message("\n=== Sites by compartment ===")
print(table(atp_combined$compartment, atp_combined$site_class))

# ---- Distance to key features ----

# Distance to central axis (approximate as center of c-ring in XY plane)
c_ring_center <- ca_df_atp[subunit == "c-ring", .(cx = mean(x), cy = mean(y))]
atp_combined[, dist_to_axis := sqrt((x - c_ring_center$cx)^2 + (y - c_ring_center$cy)^2)]

# Distance to nearest nucleotide binding site
atp_combined[, dist_to_nucleotide := sapply(1:.N, function(i) {
  min(sqrt((x[i] - nuc_coords$x)^2 + (y[i] - nuc_coords$y)^2 + (z[i] - nuc_coords$z)^2))
})]

# Distance to c-ring
atp_combined[, dist_to_c_ring := sapply(1:.N, function(i) {
  min(sqrt((x[i] - c_coords$x)^2 + (y[i] - c_coords$y)^2 + (z[i] - c_coords$z)^2))
})]

# ---- Combined statistical tests ----

message("\n=== COMBINED ANALYSIS: Wilcoxon tests vs not_sig ===")
bg_all <- atp_combined[site_class == "not_sig"]

for (dist_var in c("z", "dist_to_axis", "dist_to_nucleotide", "dist_to_c_ring")) {
  message("\n", dist_var, ":")
  for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
    test_vals <- atp_combined[site_class == sc, get(dist_var)]
    bg_vals <- bg_all[, get(dist_var)]
    test_vals <- test_vals[!is.na(test_vals)]
    bg_vals <- bg_vals[!is.na(bg_vals)]
    if (length(test_vals) < 5) next
    
    wt <- wilcox.test(test_vals, bg_vals)
    direction <- ifelse(median(test_vals, na.rm = TRUE) > median(bg_vals, na.rm = TRUE), "↑", "↓")
    message(sprintf("  %s (n=%d): median=%.1f %s (bg=%.1f), p=%.4f",
                    sc, length(test_vals), median(test_vals, na.rm = TRUE), direction,
                    median(bg_vals, na.rm = TRUE), wt$p.value))
  }
}

# ---- Combined correlations ----

message("\n=== COMBINED Spearman correlations ===")
cor_combined <- data.table(
  comparison = c("Z (stroma-lumen) vs P_only", "Z vs P_pcs",
                 "axis_distance vs P_only", "axis_distance vs P_pcs",
                 "nucleotide vs P_only", "nucleotide vs P_pcs",
                 "c_ring vs P_only", "c_ring vs P_pcs"),
  rho = c(
    cor(atp_combined$z, atp_combined$neglog_p_only, method = "spearman", use = "complete"),
    cor(atp_combined$z, atp_combined$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_axis, atp_combined$neglog_p_only, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_axis, atp_combined$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_nucleotide, atp_combined$neglog_p_only, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_nucleotide, atp_combined$neglog_p_pcs, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_c_ring, atp_combined$neglog_p_only, method = "spearman", use = "complete"),
    cor(atp_combined$dist_to_c_ring, atp_combined$neglog_p_pcs, method = "spearman", use = "complete")
  ),
  p_value = c(
    cor.test(atp_combined$z, atp_combined$neglog_p_only, method = "spearman")$p.value,
    cor.test(atp_combined$z, atp_combined$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_axis, atp_combined$neglog_p_only, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_axis, atp_combined$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_nucleotide, atp_combined$neglog_p_only, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_nucleotide, atp_combined$neglog_p_pcs, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_c_ring, atp_combined$neglog_p_only, method = "spearman")$p.value,
    cor.test(atp_combined$dist_to_c_ring, atp_combined$neglog_p_pcs, method = "spearman")$p.value
  )
)
print(cor_combined)

# ---- Compartment enrichment ----

message("\n=== Compartment enrichment (Fisher's exact) ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  message("\n", sc, ":")
  for (comp in c("stroma", "membrane", "lumen")) {
    sig_in <- sum(atp_combined$site_class == sc & atp_combined$compartment == comp, na.rm = TRUE)
    sig_out <- sum(atp_combined$site_class == sc & atp_combined$compartment != comp, na.rm = TRUE)
    bg_in <- sum(atp_combined$site_class == "not_sig" & atp_combined$compartment == comp, na.rm = TRUE)
    bg_out <- sum(atp_combined$site_class == "not_sig" & atp_combined$compartment != comp, na.rm = TRUE)
    
    if (sig_in + sig_out == 0) next
    
    mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
    ft <- fisher.test(mat)
    
    prop_sig <- sig_in / (sig_in + sig_out)
    prop_bg <- bg_in / (bg_in + bg_out)
    
    message(sprintf("  %s: %.1f%% (vs %.1f%% bg), OR=%.2f, p=%.4f",
                    comp, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
  }
}

# ---- Combined secondary structure ----

# Merge with struct_var for all genes
atp_combined_ss <- merge(atp_combined, 
                         struct_var[Gene %in% atp_genes, .(Gene, Position, consensus)],
                         by = c("Gene", "Position"), all.x = TRUE)

message("\n=== COMBINED Secondary structure by site class ===")
print(table(atp_combined_ss$consensus, atp_combined_ss$site_class, useNA = "ifany"))

ss_combined <- atp_combined_ss[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_combined[, total := sum(N), by = site_class]
ss_combined[, prop := N / total]
print(dcast(ss_combined, site_class ~ consensus, value.var = "prop"))

message("\n=== COMBINED Coil enrichment ===")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus != "C", na.rm = TRUE)
  
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  
  prop_sig <- sig_C / (sig_C + sig_notC)
  prop_bg <- bg_C / (bg_C + bg_notC)
  
  message(sprintf("%s: %.1f%% C (vs %.1f%% bg), OR=%.2f, p=%.4f",
                  sc, prop_sig * 100, prop_bg * 100, ft$estimate, ft$p.value))
}

# ---- Visualization: all subunits on structure ----

p_combined_xy <- ggplot() +
  geom_point(data = ca_df_atp, aes(x = x, y = y), color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = atp_combined[site_class == "not_sig"], aes(x = x, y = y),
             color = "gray50", size = 1) +
  geom_point(data = atp_combined[site_class == "sig_no_control"], aes(x = x, y = y),
             color = "gold", size = 2) +
  geom_point(data = atp_combined[site_class == "sig_with_control"], aes(x = x, y = y),
             color = "steelblue", size = 2) +
  geom_point(data = atp_combined[site_class == "sig_both"], aes(x = x, y = y),
             color = "darkred", size = 3) +
  geom_point(data = ligand_centers_atp[ligand %in% c("ATP", "ADP")],
             aes(x = lig_x, y = lig_y), color = "black", size = 4, shape = 18) +
  coord_fixed() +
  labs(title = "All ATP synthase GWAS hits - XY (top view)",
       subtitle = "Gold=sig_no_ctrl, Blue=sig_with_ctrl, Red=sig_both",
       x = "X (Å)", y = "Y (Å)") +
  theme_classic()

p_combined_xz <- ggplot() +
  geom_point(data = ca_df_atp, aes(x = x, y = z), color = "gray90", size = 0.3, alpha = 0.3) +
  geom_point(data = atp_combined[site_class == "not_sig"], aes(x = x, y = z),
             color = "gray50", size = 1) +
  geom_point(data = atp_combined[site_class == "sig_no_control"], aes(x = x, y = z),
             color = "gold", size = 2) +
  geom_point(data = atp_combined[site_class == "sig_with_control"], aes(x = x, y = z),
             color = "steelblue", size = 2) +
  geom_point(data = atp_combined[site_class == "sig_both"], aes(x = x, y = z),
             color = "darkred", size = 3) +
  geom_hline(yintercept = c(105, 135), linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = min(ca_df_atp$x), y = 145, label = "Stroma", hjust = 0, size = 3) +
  annotate("text", x = min(ca_df_atp$x), y = 120, label = "Membrane", hjust = 0, size = 3) +
  annotate("text", x = min(ca_df_atp$x), y = 95, label = "Lumen", hjust = 0, size = 3) +
  coord_fixed() +
  labs(title = "All ATP synthase GWAS hits - XZ (side view)",
       subtitle = "Dashed lines = membrane boundaries",
       x = "X (Å)", y = "Z (Å)") +
  theme_classic()

p_combined_xy + p_combined_xz


# ---- HYPOTHESIS TESTING: ATP SYNTHASE GWAS SIGNAL ----

# ============================================================
# H1: GWAS hits are enriched in the STROMA (F1 head) vs membrane/lumen
#     Rationale: F1 is the catalytic portion, more likely under selection
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("H1: GWAS hits enriched in STROMA (F1 catalytic head)")
message(paste(rep("=", 60), collapse = ""))

# Test: are sig sites higher in Z (more stromal)?
message("\nZ-coordinate by site class:")
print(atp_combined[, .(n = .N, median_z = median(z, na.rm=T), mean_z = mean(z, na.rm=T)), by = site_class])

# Wilcoxon
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_z <- atp_combined[site_class == sc, z]
  bg_z <- atp_combined[site_class == "not_sig", z]
  if (length(test_z) < 3) next
  wt <- wilcox.test(test_z, bg_z)
  message(sprintf("%s: median Z = %.1f vs %.1f (bg), p = %.4f", 
                  sc, median(test_z), median(bg_z), wt$p.value))
}

# Fisher's test: stroma vs not-stroma
message("\nStroma enrichment (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_stroma <- sum(atp_combined$site_class == sc & atp_combined$compartment == "stroma", na.rm = TRUE)
  sig_other <- sum(atp_combined$site_class == sc & atp_combined$compartment != "stroma", na.rm = TRUE)
  bg_stroma <- sum(atp_combined$site_class == "not_sig" & atp_combined$compartment == "stroma", na.rm = TRUE)
  bg_other <- sum(atp_combined$site_class == "not_sig" & atp_combined$compartment != "stroma", na.rm = TRUE)
  
  if (sig_stroma + sig_other == 0) next
  mat <- matrix(c(sig_stroma, sig_other, bg_stroma, bg_other), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d in stroma (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_stroma, sig_stroma + sig_other,
                  100 * sig_stroma / (sig_stroma + sig_other),
                  100 * bg_stroma / (bg_stroma + bg_other),
                  ft$estimate, ft$p.value))
}

# ============================================================
# H2: GWAS hits are PERIPHERAL (far from central axis)
#     Rationale: Core (γ, c-ring axis) is ultra-conserved for rotation
#                Periphery (α-β interfaces, peripheral stalk) may tolerate variation
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("H2: GWAS hits are PERIPHERAL (far from rotation axis)")
message(paste(rep("=", 60), collapse = ""))

message("\nDistance to central axis by site class:")
print(atp_combined[, .(n = .N, median_dist = median(dist_to_axis, na.rm=T)), by = site_class])

for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_dist <- atp_combined[site_class == sc, dist_to_axis]
  bg_dist <- atp_combined[site_class == "not_sig", dist_to_axis]
  if (length(test_dist) < 3) next
  wt <- wilcox.test(test_dist, bg_dist)
  direction <- ifelse(median(test_dist, na.rm=T) > median(bg_dist, na.rm=T), "MORE peripheral", "MORE central")
  message(sprintf("%s: median = %.1f Å vs %.1f Å (bg), %s, p = %.4f",
                  sc, median(test_dist, na.rm=T), median(bg_dist, na.rm=T), direction, wt$p.value))
}

# Correlation
cor_axis <- cor.test(atp_combined$dist_to_axis, atp_combined$neglog_p_only, method = "spearman")
message(sprintf("\nCorrelation (axis dist vs -log10 P_only): rho = %.3f, p = %.2e", 
                cor_axis$estimate, cor_axis$p.value))

# ============================================================
# H3: GWAS hits AVOID the nucleotide binding sites
#     Rationale: ATP/ADP binding pockets are ultra-conserved catalytic residues
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("H3: GWAS hits AVOID nucleotide binding sites")
message(paste(rep("=", 60), collapse = ""))

message("\nDistance to nearest nucleotide by site class:")
print(atp_combined[, .(n = .N, median_dist = median(dist_to_nucleotide, na.rm=T)), by = site_class])

for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_dist <- atp_combined[site_class == sc, dist_to_nucleotide]
  bg_dist <- atp_combined[site_class == "not_sig", dist_to_nucleotide]
  if (length(test_dist) < 3) next
  wt <- wilcox.test(test_dist, bg_dist)
  direction <- ifelse(median(test_dist, na.rm=T) > median(bg_dist, na.rm=T), "FARTHER (avoiding)", "CLOSER")
  message(sprintf("%s: median = %.1f Å vs %.1f Å (bg), %s, p = %.4f",
                  sc, median(test_dist, na.rm=T), median(bg_dist, na.rm=T), direction, wt$p.value))
}

# Define "near nucleotide" as < 15 Å
atp_combined[, near_nucleotide := dist_to_nucleotide < 15]
message("\nNear nucleotide (<15Å) depletion:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_near <- sum(atp_combined$site_class == sc & atp_combined$near_nucleotide == TRUE, na.rm = TRUE)
  sig_far <- sum(atp_combined$site_class == sc & atp_combined$near_nucleotide == FALSE, na.rm = TRUE)
  bg_near <- sum(atp_combined$site_class == "not_sig" & atp_combined$near_nucleotide == TRUE, na.rm = TRUE)
  bg_far <- sum(atp_combined$site_class == "not_sig" & atp_combined$near_nucleotide == FALSE, na.rm = TRUE)
  
  if (sig_near + sig_far == 0) next
  mat <- matrix(c(sig_near, sig_far, bg_near, bg_far), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d near nucleotide (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_near, sig_near + sig_far,
                  100 * sig_near / (sig_near + sig_far),
                  100 * bg_near / (bg_near + bg_far),
                  ft$estimate, ft$p.value))
}

# ============================================================
# H4: GWAS hits are enriched in COILS (loops) vs structured regions
#     Rationale: Loops are flexible, can tolerate variation
#                α-helices and β-sheets are structurally constrained
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("H4: GWAS hits enriched in COILS (flexible loops)")
message(paste(rep("=", 60), collapse = ""))

message("\nSecondary structure proportions by site class:")
ss_props <- atp_combined_ss[!is.na(consensus), .N, by = .(site_class, consensus)]
ss_props[, total := sum(N), by = site_class]
ss_props[, prop := round(N / total * 100, 1)]
print(dcast(ss_props, site_class ~ consensus, value.var = "prop"))

message("\nCoil enrichment (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_C <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus == "C", na.rm = TRUE)
  sig_notC <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus != "C", na.rm = TRUE)
  bg_C <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus == "C", na.rm = TRUE)
  bg_notC <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus != "C", na.rm = TRUE)
  
  if (sig_C + sig_notC == 0) next
  mat <- matrix(c(sig_C, sig_notC, bg_C, bg_notC), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d coil (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_C, sig_C + sig_notC,
                  100 * sig_C / (sig_C + sig_notC),
                  100 * bg_C / (bg_C + bg_notC),
                  ft$estimate, ft$p.value))
}

# Also test helix depletion
message("\nHelix depletion:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_H <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus == "H", na.rm = TRUE)
  sig_notH <- sum(atp_combined_ss$site_class == sc & atp_combined_ss$consensus != "H", na.rm = TRUE)
  bg_H <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus == "H", na.rm = TRUE)
  bg_notH <- sum(atp_combined_ss$site_class == "not_sig" & atp_combined_ss$consensus != "H", na.rm = TRUE)
  
  if (sig_H + sig_notH == 0) next
  mat <- matrix(c(sig_H, sig_notH, bg_H, bg_notH), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d helix (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_H, sig_H + sig_notH,
                  100 * sig_H / (sig_H + sig_notH),
                  100 * bg_H / (bg_H + bg_notH),
                  ft$estimate, ft$p.value))
}

# ============================================================
# H5: GWAS hits cluster at SUBUNIT INTERFACES
#     Rationale: Interfaces mediate allostery and assembly
#                Variation here could affect enzyme kinetics
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("H5: GWAS hits enriched at SUBUNIT INTERFACES")
message(paste(rep("=", 60), collapse = ""))

# Calculate distance to nearest OTHER subunit for each site
# Use the full complex to find nearest non-self chain

atp_combined[, dist_to_interface := sapply(1:.N, function(i) {
  my_gene <- gene[i]
  # Get all other chains' coordinates
  if (my_gene == "atpA") other_coords <- ca_df_atp[!subunit %in% c("alpha")]
  else if (my_gene == "atpB") other_coords <- ca_df_atp[!subunit %in% c("beta")]
  else if (my_gene == "atpE") other_coords <- ca_df_atp[chain != "e"]
  else if (my_gene == "atpF") other_coords <- ca_df_atp[chain != "b"]
  else if (my_gene == "atpH") other_coords <- ca_df_atp[!subunit %in% c("c-ring")]
  else if (my_gene == "atpI") other_coords <- ca_df_atp[chain != "a"]
  else other_coords <- ca_df_atp
  
  min(sqrt((x[i] - other_coords$x)^2 + (y[i] - other_coords$y)^2 + (z[i] - other_coords$z)^2))
})]

# Define interface as < 8 Å from another subunit
atp_combined[, at_interface := dist_to_interface < 8]

message("\nInterface proximity by site class:")
print(atp_combined[, .(n = .N, 
                       median_interface_dist = median(dist_to_interface, na.rm=T),
                       pct_at_interface = 100 * mean(at_interface, na.rm=T)), 
                   by = site_class])

message("\nInterface enrichment (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_interface <- sum(atp_combined$site_class == sc & atp_combined$at_interface == TRUE, na.rm = TRUE)
  sig_bulk <- sum(atp_combined$site_class == sc & atp_combined$at_interface == FALSE, na.rm = TRUE)
  bg_interface <- sum(atp_combined$site_class == "not_sig" & atp_combined$at_interface == TRUE, na.rm = TRUE)
  bg_bulk <- sum(atp_combined$site_class == "not_sig" & atp_combined$at_interface == FALSE, na.rm = TRUE)
  
  if (sig_interface + sig_bulk == 0) next
  mat <- matrix(c(sig_interface, sig_bulk, bg_interface, bg_bulk), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d at interface (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_interface, sig_interface + sig_bulk,
                  100 * sig_interface / (sig_interface + sig_bulk),
                  100 * bg_interface / (bg_interface + bg_bulk),
                  ft$estimate, ft$p.value))
}

# Correlation: closer to interface = more significant?
cor_interface <- cor.test(atp_combined$dist_to_interface, atp_combined$neglog_p_only, method = "spearman")
message(sprintf("\nCorrelation (interface dist vs -log10 P_only): rho = %.3f, p = %.2e",
                cor_interface$estimate, cor_interface$p.value))

# ============================================================
# SUMMARY PLOT
# ============================================================

p_h1 <- ggplot(atp_combined, aes(x = site_class, y = z, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  geom_hline(yintercept = c(105, 135), linetype = "dashed", alpha = 0.5) +
  labs(title = "H1: Stromal enrichment?", y = "Z coordinate (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h2 <- ggplot(atp_combined, aes(x = site_class, y = dist_to_axis, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H2: Peripheral?", y = "Distance to axis (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h3 <- ggplot(atp_combined, aes(x = site_class, y = dist_to_nucleotide, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H3: Avoid nucleotide?", y = "Distance to ATP/ADP (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_h5 <- ggplot(atp_combined, aes(x = site_class, y = dist_to_interface, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H5: At interfaces?", y = "Distance to interface (Å)", x = "") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

(p_h1 + p_h2) / (p_h3 + p_h5)

# ---- ligand binding sites ----
# H6: GWAS hits AVOID ALL LIGAND BINDING SITES
#     Rationale: Ligand-binding residues (substrates, Mg, cofactors)
#                are functionally critical and highly conserved

message("\n", paste(rep("=", 60), collapse = ""))
message("H6: GWAS hits AVOID all ligand binding sites")
message(paste(rep("=", 60), collapse = ""))

# Get all non-protein residues (ligands, cofactors, ions)
protein_aa <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
                "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
ligand_atoms <- pdb_atp$atom[!pdb_atp$atom$resid %in% protein_aa, ]

message("\nLigands found in structure:")
print(table(ligand_atoms$resid))

# Get ligand centers (exclude water/HOH if present)
ligand_dt <- as.data.table(ligand_atoms)
ligand_dt <- ligand_dt[!resid %in% c("HOH", "WAT")]  # exclude water

ligand_centers_all <- ligand_dt[, .(
  lig_x = mean(x), lig_y = mean(y), lig_z = mean(z)
), by = .(chain, resid)]

message("\nLigand centers (excluding water):")
print(ligand_centers_all)

# Calculate distance to nearest ligand for each site
atp_combined[, dist_to_any_ligand := sapply(1:.N, function(i) {
  min(sqrt((x[i] - ligand_centers_all$lig_x)^2 + 
             (y[i] - ligand_centers_all$lig_y)^2 + 
             (z[i] - ligand_centers_all$lig_z)^2))
})]

message("\nDistance to nearest ligand by site class:")
print(atp_combined[, .(n = .N, 
                       median_dist = median(dist_to_any_ligand, na.rm=T),
                       mean_dist = mean(dist_to_any_ligand, na.rm=T)), 
                   by = site_class])

# Wilcoxon tests
message("\nWilcoxon tests vs not_sig:")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  test_dist <- atp_combined[site_class == sc, dist_to_any_ligand]
  bg_dist <- atp_combined[site_class == "not_sig", dist_to_any_ligand]
  if (length(test_dist) < 3) next
  wt <- wilcox.test(test_dist, bg_dist)
  direction <- ifelse(median(test_dist, na.rm=T) > median(bg_dist, na.rm=T), 
                      "FARTHER (avoiding)", "CLOSER")
  message(sprintf("%s: median = %.1f Å vs %.1f Å (bg), %s, p = %.4f",
                  sc, median(test_dist, na.rm=T), median(bg_dist, na.rm=T), 
                  direction, wt$p.value))
}

# Define "near ligand" as < 10 Å
atp_combined[, near_any_ligand := dist_to_any_ligand < 10]

message("\nNear any ligand (<10Å) depletion (Fisher's exact):")
for (sc in c("sig_no_control", "sig_with_control", "sig_both")) {
  sig_near <- sum(atp_combined$site_class == sc & atp_combined$near_any_ligand == TRUE, na.rm = TRUE)
  sig_far <- sum(atp_combined$site_class == sc & atp_combined$near_any_ligand == FALSE, na.rm = TRUE)
  bg_near <- sum(atp_combined$site_class == "not_sig" & atp_combined$near_any_ligand == TRUE, na.rm = TRUE)
  bg_far <- sum(atp_combined$site_class == "not_sig" & atp_combined$near_any_ligand == FALSE, na.rm = TRUE)
  
  if (sig_near + sig_far == 0) next
  mat <- matrix(c(sig_near, sig_far, bg_near, bg_far), nrow = 2)
  ft <- fisher.test(mat)
  message(sprintf("%s: %d/%d near ligand (%.1f%%) vs bg %.1f%%, OR=%.2f, p=%.4f",
                  sc, sig_near, sig_near + sig_far,
                  100 * sig_near / (sig_near + sig_far),
                  100 * bg_near / (bg_near + bg_far),
                  ft$estimate, ft$p.value))
}

# Correlation
cor_ligand <- cor.test(atp_combined$dist_to_any_ligand, atp_combined$neglog_p_only, 
                       method = "spearman")
message(sprintf("\nCorrelation (ligand dist vs -log10 P_only): rho = %.3f, p = %.2e",
                cor_ligand$estimate, cor_ligand$p.value))
plot(atp_combined$dist_to_any_ligand[atp_combined$site_class!="not_sig"],
     atp_combined$neglog_p_pcs[atp_combined$site_class!="not_sig"],
     xlab="distance from any ligand (Angstrom)",ylab="-log(p_pcs)")
plot(atp_combined$dist_to_any_ligand[atp_combined$site_class!="not_sig"],
     atp_combined$neglog_p_only[atp_combined$site_class!="not_sig"],
     xlab="distance from any ligand (Angstrom)",ylab="-log(p)")
# Add to summary plot
p_h6 <- ggplot(atp_combined, aes(x = site_class, y = dist_to_any_ligand, fill = site_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("not_sig" = "gray60", "sig_no_control" = "gold",
                               "sig_with_control" = "steelblue", "sig_both" = "darkred")) +
  labs(title = "H6: Avoid ligands?", y = "Distance to any ligand (Å)", x = "") +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Updated summary plot with H6
(p_h1 + p_h2 + p_h3) / (p_h5 + p_h6)
