#find pairs with minimal AA differences but reasonable temp differences 
library(data.table)
library(arrow)
data <- as.data.table(read_parquet("data/processed_data.parquet"))
head(data$Organism)

aa_codes <- c("A", "R", "N", "D", "C", 
              "Q", "E", "G", "H", "I", 
              "L", "K", "M", "F", "P", 
              "S", "T", "W", "Y", "V")

y <- data$pheno_Topt_site_p50
aa_dat <- data[,..aa_codes]

#find pairwise distances between y and AA proportions, and plot the one against another 

aa_dist <- dist(aa_dat, method = "euclidean")
y_dist <- dist(y, method = "euclidean")

aa_dist_vec <- as.vector(aa_dist)
y_dist_vec <- as.vector(y_dist) * 0.01
par(mfrow=c(2,1))
hist(aa_dist_vec, main="AA prop distances")
abline(v=quantile(aa_dist_vec,0.01),col="red")
hist(y_dist_vec, main="Topt_site_p50 distances")
abline(v=quantile(y_dist_vec,0.5),col="blue")

aa_cutoff <- quantile(aa_dist_vec,0.01)
y_cutoff <- quantile(y_dist_vec,0.5)


plot_df <- data.frame(
  AA_distance = aa_dist_vec,
  Phenotype_distance = y_dist_vec
)

candidate_pairs <- plot_df[which(
  plot_df$AA_distance < aa_cutoff &
    plot_df$Phenotype_distance > y_cutoff
), ]

par(mfrow=c(1,1))
plot(candidate_pairs$AA_distance, candidate_pairs$Phenotype_distance)

# Extract genus from binomial name
data[, genus := sapply(strsplit(Organism, " "), `[`, 1)]
sum(is.na(data$genus))

# Get unique genera with sufficient data
genus_counts <- data[, .N, by = genus]
sister_genera <- genus_counts[N >= 2, genus]  # At least 2 species per genus

print(paste("Found", length(sister_genera), "genera with multiple species"))

aa_codes <- c("A", "R", "N", "D", "C", 
              "Q", "E", "G", "H", "I", 
              "L", "K", "M", "F", "P", 
              "S", "T", "W", "Y", "V")

# Filter data to sister genera only
sister_data <- data[genus %in% sister_genera]

y <- sister_data$pheno_Topt_site_p50
aa_dat <- sister_data[, ..aa_codes]

# Calculate pairwise distances
aa_dist <- dist(aa_dat, method = "euclidean")
y_dist <- dist(y, method = "euclidean")
aa_dist_vec <- as.vector(aa_dist)
y_dist_vec <- as.vector(y_dist) * 0.01

# Create pair indices for distance matrices
n <- nrow(sister_data)
pairs_idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)

# Create dataframe with genus pairs and distances
genus_pairs <- data.frame(
  genus1 = sister_data$genus[pairs_idx[,1]],
  genus2 = sister_data$genus[pairs_idx[,2]],
  organism1 = sister_data$Organism[pairs_idx[,1]],
  organism2 = sister_data$Organism[pairs_idx[,2]],
  AA_distance = aa_dist_vec,
  Phenotype_distance = y_dist_vec,
  stringsAsFactors = FALSE
)

# Filter for sister genera pairs (same genus)
sister_pairs <- genus_pairs[genus_pairs$genus1 == genus_pairs$genus2, ]

print(paste("Found", nrow(sister_pairs), "sister pairs"))
hist(genus_pairs$AA_distance, main="AA prop distances (sister pairs)")
hist(sister_pairs$AA_distance, main="AA prop distances (sister pairs)", xlim=range(genus_pairs$AA_distance))

# Visualize distributions
par(mfrow=c(2,1))
hist(sister_pairs$AA_distance, main="AA prop distances (sister pairs)")
hist(sister_pairs$Phenotype_distance, main="Topt_site_p50 distances (sister pairs)")

# Plot relationship
par(mfrow=c(1,1))
plot(sister_pairs$AA_distance, sister_pairs$Phenotype_distance,
     xlab="AA Distance", ylab="Phenotype Distance",
     main="Sister Pairs: AA vs Phenotype Distance")

aa_cutoff <- quantile(sister_pairs$AA_distance,0.02)
y_cutoff <- quantile(sister_pairs$Phenotype_distance,0.5)
points(sister_pairs[sister_pairs$AA_distance < aa_cutoff & sister_pairs$Phenotype_distance > y_cutoff, "AA_distance"],
       sister_pairs[sister_pairs$AA_distance < aa_cutoff & sister_pairs$Phenotype_distance > y_cutoff, "Phenotype_distance"],
       col="red")

candidates <- sister_pairs[sister_pairs$AA_distance < aa_cutoff & sister_pairs$Phenotype_distance > y_cutoff,]
length(unique(candidates$organism1))
length(unique(candidates$organism2))

length(c(unique(candidates$organism2), unique(candidates$organism1)))

# Summary statistics
cat("\nSister pairs summary:\n")
print(summary(sister_pairs[, c("AA_distance", "Phenotype_distance")]))

# Show some example pairs
cat("\nExample sister pairs:\n")
print(head(sister_pairs[, c("organism1", "organism2", "AA_distance", "Phenotype_distance")]))



