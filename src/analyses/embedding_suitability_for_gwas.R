library(arrow)
library(Biostrings)
library(tidyverse)
library(lsa)
library(stringr)
library(dplyr)

# Load shared data once
data <- read_parquet("data/processed_data.parquet")
# Find matching embedding and alignment files
emb <- read_parquet("data/embeddings/psbA_residue_embeddings.parquet")

# Read alignment
aln <- readAAStringSet("data/tmp/alignedGenes/psbA_AA_aligned.fasta")

# Get names and sequences
seq_names <- names(aln)
seqs <- as.character(aln)

# Store in a tibble
aln_df <- tibble(ID = seq_names, seq = seqs)

aln_mat <- seqs %>%
  strsplit(split = "") %>%
  do.call(rbind, .)

rownames(aln_mat) <- seq_names
aln_mat <- as_tibble(aln_mat, rownames = "ID")

ids <- aln_mat$ID
mat <- as.matrix(aln_mat[,-1])

# Compute true (ungapped) sequence length for each alignment row
seq_lengths <- apply(mat, 1, function(x) sum(x != "-"))
# Map: alignment position → residue index for each sequence
align_to_res_idx <- lapply(1:nrow(mat), function(i) {
  aln_seq <- mat[i, ]
  idx_map <- cumsum(aln_seq != "-")  # running count of residues
  idx_map[aln_seq == "-"] <- NA      # gaps have no residue index
  idx_map
})

# Function to count differences (ignoring gaps)
count_diffs <- function(seq1, seq2) {
  sum(seq1 != seq2 & seq1 != "-" & seq2 != "-")
}

# Pick a random reference sequence
for (i in 1:5) {
  set.seed(i*10)
  
  # Random reference sequence
  ref_idx <- sample(seq_len(nrow(mat)), 1)
  ref_id <- ids[ref_idx]
  ref_seq <- mat[ref_idx, ]
  
  # Restrict to same ungapped length (simplify coordinate logic)
  ref_len <- seq_lengths[ref_idx]
  valid_candidates <- which(seq_lengths == ref_len)
  
  # Find first sequence differing at exactly one amino acid
  target_idx <- NA
  for (j in valid_candidates) {
    if (j == ref_idx) next
    diffs <- count_diffs(ref_seq, mat[j, ])
    if (diffs == 2) {
      target_idx <- j
      break
    }
  }
  
  if (is.na(target_idx)) {
    cat("No sequence found that differs by exactly one residue.\n")
    next
  }
  
  cat("Reference sequence:", ref_id, "\n")
  cat("Matching sequence:", ids[target_idx], "\n")
  
  # Alignment position of difference
  diff_pos <- which(ref_seq != mat[target_idx, ] & ref_seq != "-" & mat[target_idx, ] != "-")
  cat("Differing alignment position:", diff_pos, "\n")
  cat("Residues:", ref_seq[diff_pos], "vs", mat[target_idx, diff_pos], "\n")
  
  # Convert alignment coordinate → residue index
  ref_res_idx <- align_to_res_idx[[ref_idx]][diff_pos]
  target_res_idx <- align_to_res_idx[[target_idx]][diff_pos]
  
  cat("Residue indices in embeddings:", ref_res_idx, "(ref) vs", target_res_idx, "(target)\n")
  
  # Extract the IDs before the first "|"
  aln_ids <- tibble(
    full_id = ids,
    ID = str_extract(ids, "^[^|]+")
  )
  
  ref_clean <- str_extract(ref_id, "^[^|]+")
  target_clean <- str_extract(ids[target_idx], "^[^|]+")
  
  cat("Ref ID for embeddings:", ref_clean, "\n")
  cat("Target ID for embeddings:", target_clean, "\n")
  
  # Embedding data
  emb_ref <- emb %>% filter(ID == ref_clean)
  emb_target <- emb %>% filter(ID == target_clean)
  
  if (nrow(emb_ref) != nrow(emb_target)) {
    cat("Warning: unequal embedding lengths — skipping pair.\n")
    next
  }
  
  emb_cols <- grep("embedding", colnames(emb_ref))
  
  # Compute cosine similarity residue-wise
  emb_diff <- inner_join(
    emb_ref %>% select(Residue_Index, all_of(emb_cols)),
    emb_target %>% select(Residue_Index, all_of(emb_cols)),
    by = "Residue_Index",
    suffix = c("_ref", "_target")
  ) %>%
    rowwise() %>%
    mutate(cos_sim = cosine(
      as.numeric(c_across(ends_with("_ref"))),
      as.numeric(c_across(ends_with("_target")))
    )) %>%
    ungroup()
  
  # Plot cosine similarity with true residue index marked
  p <- ggplot(emb_diff, aes(x = Residue_Index, y = cos_sim)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_vline(xintercept = ref_res_idx, color = "red", linetype = "dashed") +
    labs(title = paste0("psbA embedding cosine similarity"),#, ref_clean, " vs ", target_clean),
         subtitle = paste("Red dashed line = site of amino acid differences"),
         x = "Residue index",
         y = "Cosine similarity between embeddings") +
    theme_minimal(base_size = 13)
  
  print(p)
}

emb_cols <- grep("embedding", colnames(emb))
pos_var <- emb %>%
  group_by(Residue_Index) %>%
  summarize(
    mean_cosine_dist = 1 - mean(combn(seq_len(n()), 2, function(p)
      cosine(as.numeric(emb_cols[p[1]]), as.numeric(emb_cols[p[2]]))
    )),
    .groups = "drop"
  )

ggplot(pos_var, aes(x = Residue_Index, y = mean_cosine_dist)) +
  geom_line(color = "darkorange") +
  labs(y = "Mean pairwise cosine distance", x = "Residue index")

### PCs

library(lsa)
library(ggplot2)
library(dplyr)
library(stringr)

genes_to_process <- c("ndhG")
for (gene in genes_to_process) {
  message("\n=== PC cosine similarity analysis for gene: ", gene, " ===")
  
  emb_file <- file.path("data/embeddings/", paste0(gene, "_residue_embeddings.parquet"))
  aln_file <- file.path("data/tmp/alignedGenes/", paste0(gene, "_AA_aligned.fasta"))
  
  # --- Read data ---
  df_emb <- tryCatch(as.data.table(read_parquet(emb_file)), error = function(e) NULL)
  aln <- tryCatch(readAAStringSet(aln_file), error = function(e) NULL)
  if (is.null(df_emb) || is.null(aln) || length(aln) == 0) next
  
  # --- Clean IDs ---
  names(aln) <- sub("\\|.*", "", names(aln))
  df_emb[, ID := sub("\\|.*", "", ID)]
  
  # --- Build alignment matrix ---
  aln_mat <- as.matrix(aln)
  seq_ids <- names(aln)
  
  # --- Per-sequence ungapped length and position mapping ---
  seq_lengths <- apply(aln_mat, 1, function(x) sum(x != "-"))
  align_to_res_idx <- lapply(1:nrow(aln_mat), function(i) {
    aln_seq <- aln_mat[i, ]
    idx_map <- cumsum(aln_seq != "-")
    idx_map[aln_seq == "-"] <- NA
    idx_map
  })
  
  # --- Helper: count differences ignoring gaps ---
  count_diffs <- function(seq1, seq2) sum(seq1 != seq2 & seq1 != "-" & seq2 != "-")
  
  # --- Try 5 random references per gene ---
  for (i in 1:5) {
    set.seed(i)
    ref_idx <- sample(seq_len(nrow(aln_mat)), 1)
    ref_seq <- aln_mat[ref_idx, ]
    ref_id <- seq_ids[ref_idx]
    
    ref_len <- seq_lengths[ref_idx]
    valid <- which(seq_lengths == ref_len)
    
    # Find first sequence differing by 3 residues
    target_idx <- NA
    for (j in valid) {
      if (j == ref_idx) next
      if (count_diffs(ref_seq, aln_mat[j, ]) == 3) {
        target_idx <- j
        break
      }
    }
    if (is.na(target_idx)) next
    
    target_id <- seq_ids[target_idx]
    diff_pos <- which(ref_seq != aln_mat[target_idx, ] & ref_seq != "-" & aln_mat[target_idx, ] != "-")
    
    cat("Ref:", ref_id, "Target:", target_id, "Diff positions:", diff_pos, "\n")
    
    ref_res_idx <- align_to_res_idx[[ref_idx]][diff_pos]
    target_res_idx <- align_to_res_idx[[target_idx]][diff_pos]
    
    # --- Get embeddings for each sequence ---
    emb_ref <- df_emb %>% filter(ID == ref_id)
    emb_target <- df_emb %>% filter(ID == target_id)
    if (nrow(emb_ref) != nrow(emb_target)) next
    
    emb_cols <- grep("^embedding_", names(emb_ref), value = TRUE)
    emb_mat_ref <- as.matrix(emb_ref[, ..emb_cols])
    emb_mat_target <- as.matrix(emb_target[, ..emb_cols])
    
    # --- PCA on combined embeddings to extract PCs ---
    all_emb <- scale(rbind(emb_mat_ref, emb_mat_target))
    pca <- prcomp(all_emb, rank. = 10, scale. = FALSE)
    pcs <- pca$x
    pcs_ref <- pcs[1:nrow(emb_mat_ref), ]
    pcs_target <- pcs[(nrow(emb_mat_ref) + 1):nrow(pcs), ]
    
    # --- Compute cosine similarity residue-wise across PCs ---
    cos_sim_df <- tibble(
      Residue_Index = emb_ref$Residue_Index,
      cos_sim = purrr::map_dbl(seq_len(nrow(pcs_ref)), function(k)
        cosine(as.numeric(pcs_ref[k, ]), as.numeric(pcs_target[k, ]))
      )
    )
    
    # --- Plot similarity profile ---
    p <- ggplot(cos_sim_df, aes(x = Residue_Index, y = cos_sim)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_vline(xintercept = ref_res_idx, color = "red", linetype = "dashed") +
      labs(
        title = paste0("Cosine similarity across top 10 PCs (", gene, ")"),
        subtitle = paste(ref_id, "vs", target_id, "| Diff residue:", paste(ref_res_idx, collapse = ",")),
        x = "Residue index", y = "Cosine similarity (10D PC space)"
      ) +
      theme_minimal(base_size = 13)
    
    print(p)
  }
}











