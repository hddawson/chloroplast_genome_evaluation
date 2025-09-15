library(ape)
library(data.table)
library(arrow)
library(dplyr)
library(ape)
library(adephylo)
library(geiger)
library(phylobase)

data <- as.data.frame(read_parquet("data/processed_data.parquet"))
tree <- read.tree("data/Fasttree/Fasttree.nwk")
orders <- unique(data$Order)

# Storage for overlay plot
all_results <- list()

for (ord in orders) {
  cat("\n=== Processing order:", ord, "===\n")
  
  # Fix the filtering - should use ord, not hardcoded "Spermatophyta"
  order_data <- data[data$Order == ord, ]
  
  
  # Match tree to data
  order_tree <- keep.tip(tree, order_data$ID)
  
  # Ensure data order matches tree
  order_data <- order_data[match(order_tree$tip.label, order_data$ID), ]
  
  # Remove NAs
  complete_cases <- complete.cases(order_data$pheno_Topt_site_p50)
  order_data <- order_data[complete_cases, ]
  
  if (sum(complete_cases) < 10) {
    cat("Skipping", ord, "- insufficient complete cases\n")
    next
  }
  
  order_tree <- keep.tip(order_tree, order_data$ID)
  
  order_tree_fixed <- order_tree
  zero_branches <- order_tree_fixed$edge.length == 0
  very_small <- order_tree_fixed$edge.length < 1e-6
  
  cat("Fixing", sum(zero_branches), "zero branches and", 
      sum(very_small), "very small branches\n")
  
  # Set minimum branch length
  min_branch <- 1e-6
  order_tree_fixed$edge.length[order_tree_fixed$edge.length < min_branch] <- min_branch
  order_tree_fixed$node.label <- NULL
  
  rownames(order_data) <- order_data$ID
  
  # convert to phylo4d object
  #p4d <- phylo4d(order_tree_fixed, order_data$pheno_Topt_site_p50)
  
  tree_clean <- order_tree_fixed
  trait <- order_data$pheno_Topt_site_p50
  names(trait) <- order_data$ID
  
  # 1. compute patristic distances
  pat <- cophenetic.phylo(tree_clean)

  # 2. set up distance classes
  n.classes <- 15
  d.vec <- as.numeric(pat[upper.tri(pat)])
  breaks <- quantile(d.vec, probs = seq(0, 1, length.out = n.classes + 1))
  breaks <- unique(breaks)

  # 3. compute Moran's I for each class
  moran.results <- vector("list", length = length(breaks)-1)
  names(moran.results) <- paste0("class_", seq_len(length(breaks)-1))
  
  for(i in seq_len(length(breaks)-1)){
    cat("Processing distance class", i, "/", length(breaks)-1, "\n")
    lo <- breaks[i]
    hi <- breaks[i+1]
    
    # W: proximity matrix
    W <- matrix(0, nrow = nrow(pat), ncol = ncol(pat), dimnames = dimnames(pat))
    idx <- which(pat >= lo & pat < hi, arr.ind = TRUE)

    for(r in seq_len(nrow(idx))){
      W[idx[r,1], idx[r,2]] <- 1
      W[idx[r,2], idx[r,1]] <- 1
    }
    
    # Skip if all zeros
    if(all(W == 0)){
      moran.results[[i]] <- list(obs = NA, pvalue = NA)
      next
    }
    
    # Prepare data
    df <- as.data.frame(trait)
    rownames(df) <- names(trait)
    
    # Run test
    moran.results[[i]] <- abouheif.moran(df, W = W, nrepet = 99, alter = "two-sided")
  }
  
  # 4. Extract results
  I <- sapply(moran.results, function(x) {
    if (length(x) == 1 && is.na(x)) return(NA)
    if (is.list(x) && "obs" %in% names(x)) return(x$obs[1])
    return(NA)
  })
  
  pval <- sapply(moran.results, function(x) {
    if (length(x) == 1 && is.na(x)) return(NA)
    if (is.list(x) && "pvalue" %in% names(x)) return(x$pvalue[1])
    return(NA)
  })
  
  # Midpoints for plotting
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  
  # Save individual results
  results <- list(
    order = ord,
    moran_I = I,
    pvalues = pval,
    distances = mid,
    breaks = breaks,
    n_tips = length(trait),
    raw_results = moran.results
  )
  
  saveRDS(results, paste0("results/", ord, "_correlogram_results.rds"))
  
  # Store for overlay plot
  all_results[[ord]] <- results
  
  # Individual plot
  png(paste0("plots/", ord, "_correlogram.png"), width = 10, height = 8, units = "in", res = 300)
  plot(mid, I, type = "b", pch = 16, main = paste0(ord, " Phylogenetic Correlogram (n=", length(trait), ")"),
       xlab = "Patristic Distance", ylab = "Moran's I", cex.main = 1.2)
  abline(h = 0, lty = 2, col = "gray")
  
  # Color significant points
  sig_points <- !is.na(pval) & pval < 0.05
  if(any(sig_points)) {
    points(mid[sig_points], I[sig_points], col = "red", pch = 16, cex = 1.2)
  }
  
  legend("topright", legend = c("p < 0.05", "non-significant"), 
         pch = 16, col = c("red", "black"), bty = "n")
  dev.off()
  
  cat("Completed", ord, "\n")
}

# Create overlay correlogram
png("plots/all_orders_correlogram_overlay.png", width = 12, height = 10, units = "in", res = 300)

# Calculate appropriate ylim based on actual data
all_I_values <- unlist(lapply(all_results, function(x) x$moran_I[!is.na(x$moran_I)]))
ylim_range <- range(all_I_values, na.rm = TRUE)
ylim_padding <- diff(ylim_range) * 0.1
ylim_final <- c(ylim_range[1] - ylim_padding, ylim_range[2] + ylim_padding)

# Set up plot with normalized x-axis (0-1)
plot(NULL, xlim = c(0, 1), ylim = ylim_final, 
     xlab = "Normalized Patristic Distance", ylab = "Moran's I",
     main = "Phylogenetic Correlograms - All Orders", cex.main = 1.4)

abline(h = 0, lty = 2, col = "gray", lwd = 2)

# Single color with alpha transparency
line_color <- rgb(0.2, 0.4, 0.8, alpha = 0.6)  # Blue with 60% transparency
# Collect all normalized data for mean trend
all_norm_data <- data.frame(dist_norm = numeric(0), moran_I = numeric(0))

# Plot each order with normalized distances
for(i in seq_along(all_results)) {
  ord <- names(all_results)[i]
  res <- all_results[[i]]
  
  # Remove NAs for plotting
  valid <- !is.na(res$moran_I) & !is.na(res$distances)
  if(sum(valid) == 0) next
  
  # Normalize distances to 0-1 scale for this order
  dist_norm <- res$distances[valid] / max(res$distances, na.rm = TRUE)
  
  lines(dist_norm, res$moran_I[valid], col = line_color, lwd = 1.5, type = "l")
  
  # Collect data for mean trend
  all_norm_data <- rbind(all_norm_data, data.frame(dist_norm = dist_norm, moran_I = res$moran_I[valid]))
}

# Calculate and plot mean trend line
dist_bins <- seq(0, 1, length.out = 50)
mean_trend <- numeric(length(dist_bins) - 1)

for(i in 1:(length(dist_bins) - 1)) {
  bin_data <- all_norm_data$moran_I[all_norm_data$dist_norm >= dist_bins[i] & all_norm_data$dist_norm < dist_bins[i+1]]
  mean_trend[i] <- ifelse(length(bin_data) > 0, mean(bin_data, na.rm = TRUE), NA)
}

bin_centers <- (dist_bins[-1] + dist_bins[-length(dist_bins)]) / 2
valid_trend <- !is.na(mean_trend)

lines(bin_centers[valid_trend], mean_trend[valid_trend], col = "red", lwd = 3, type = "l")

# Add legend
legend("topright", legend = c("Individual orders", "Mean trend"), 
       col = c(line_color, "red"), lwd = c(1.5, 3), bty = "n")

# Add text showing number of orders
text(x = 0.1, y = par("usr")[4] * 0.9, 
     labels = paste("n =", length(all_results), "orders"), cex = 1.2, font = 2)

dev.off()

# Save combined results
saveRDS(all_results, "results/all_orders_correlogram_results.rds")

cat("\n=== Analysis Complete ===\n")
cat("Individual results saved in results/\n")
cat("Individual plots saved in plots/\n")
cat("Overlay plot: plots/all_orders_correlogram_overlay.png\n")
cat("Combined results: results/all_orders_correlogram_results.rds\n")


