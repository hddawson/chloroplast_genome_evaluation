# ============================================================================
# CHLOROPLAST PROTEIN 3D STRUCTURE ANALYSIS FRAMEWORK
# For mapping GWAS temperature adaptation signals onto protein structures
# ============================================================================

library(bio3d)
library(data.table)
library(Biostrings)
library(ggplot2)
library(patchwork)

# ============================================================================
# CONFIGURATION AND DEFAULTS
# ============================================================================

#' Default complex definitions
#' Each complex specifies: PDB ID, gene-to-chain mappings, compartment boundaries
#' Chain mappings can be "auto" for sequence-based detection or explicit
DEFAULT_COMPLEXES <- list(
  
  atp_synthase = list(
    pdb = "6FQF",
    organism = "Spinacia oleracea",
    description = "Chloroplast F1Fo ATP synthase",
    # gene -> chain(s) mapping; multiple chains = homomeric subunits
    chain_map = list(
      atpA = c("A", "C", "E"),      # α subunits (F1 head)
      atpB = c("B", "D", "F"),      # β subunits (F1 head, catalytic)
      atpE = "e",                    # ε subunit (central stalk base)
      atpF = "b",                    # b subunit (peripheral stalk)
      atpH = c("G","H","I","J","K","L","M","N","O","P","Q","R","S","T"), # c-ring
      atpI = "a"                     # a subunit (proton channel)
    ),
    # Z-coordinate boundaries for compartment assignment
    compartment_z = c(lumen_stroma = 105, stroma_top = 135),
    compartment_labels = c("lumen", "membrane", "stroma")
  ),
  
  
  rubisco = list(
    pdb = "1RCX",
    organism = "Spinacia oleracea", 
    description = "Ribulose-1,5-bisphosphate carboxylase/oxygenase L8S8",
    chain_map = list(
      rbcL = c("L", "B", "E", "H", "K", "O", "R", "V")  # 8 large subunits (467 aa each)
      # rbcS (small subunits: S, C, F, I, M, P, T, W) is nuclear-encoded
    ),
    compartment_z = NULL,  # soluble enzyme, no membrane
    compartment_labels = NULL
  ),
  
  cytochrome_b6f = list(
    pdb = "6RQF",  # Spinach, 1.5Å resolution
    organism = "Spinacia oleracea",
    description = "Cytochrome b6f complex",
    chain_map = list(
      petA = c("C", "K"),   # cytochrome f ~ 285 
      petB = c("A", "I"),   # cytochrome b6 ~ 215
      # petC - Rieske Fe-S, D & L 
      petD = c("B", "J"),   # subunit IV ~ 170
      petG = c("G", "O"),   # subunit V
      petL = c("M", "E"),   # subunit VI
      petM = c("F", "N"),   # subunit VII
      petN = c("P", "H")    # subunit VIII
    ),
    compartment_z = c(60, 100),  # adjust based on structure
    compartment_labels = c("lumen", "membrane", "stroma")
  ),
  
  photosystem_ii = list(
    pdb = "5XNL",  # Thermosynechococcus, but good resolution
    organism = "Thermosynechococcus elongatus",
    description = "Photosystem II",
    chain_map = list(
      psbA = c("A", "a"),   # D1 protein
      psbB = c("B", "b"),   # CP47
      psbC = c("C", "c"),   # CP43
      psbD = c("D", "d"),    # D2 protein
      psbE = c("E", "e"),
      psbF = c("F", "f"),
      psbH = c("H", "h"),
      psbI = c("I", "i"),
      psbJ = c("J", "j"),
      psbK = c("K", "k"),
      psbM = c("M", "m"),
      #psbN = c("C", "c"), not in structure
      psbT = c("T", "t"),
      psbZ = c("Z", "z")
      # many more small subunits...
    ),
    compartment_z = c(-10, -48),
    compartment_labels = c("lumen", "membrane", "stroma")
  ),
  
  photosystem_i = list(
    pdb = "5ZJI",  # Pisum sativum (pea) 2.8A
    organism = "Pisum sativum",
    description = "Photosystem I",
    chain_map = list(
      psaA = c("A"),
      psaB = c("B"),
      psaC = c("C"),
      psaJ = c("J")
    ),
    compartment_z = c(0, 40),
    compartment_labels = c("lumen", "membrane", "stroma")
  ),
  
  ndh_complex = list( # 
    pdb = "7EU3",  # Recent cryo-EM structure
    organism = "Hordeum vulgare",
    description = "NDH-1 complex (chloroplast NADH dehydrogenase-like)",
    chain_map = list(
      ndhA = "A",
      ndhB = "B",
      ndhC = "C",
      ndhD = "D",
      ndhG = "G",
      ndhH = "H",
      ndhI = "I",
      ndhJ = "J",
      ndhK = "K"
    ),
    compartment_z = NULL,
    compartment_labels = NULL
  ),
  rps_complex = list( # 
    pdb = "5MMJ",  # Recent cryo-EM structure https://www.nature.com/articles/srep35793
    organism = "Spinacia oleracea",
    description = "small subunit of plastid ribosome",
    chain_map = list( #usually follow a nice SX pattern in the molecule name, map to alphabet
      rps3 = "c", #E auth c, #P09595
      rps4 = "d", #F auth d P13788
      rps7 = "g", # I auth G P82129, a and B?
      rps8 = "h", # J auth H P09597
      rps11 = "k", #M auth k P06506
      rps14 = "n", # P auth N P0656 # two different forms? 
      rps18 = "r", # T auth r Q9M3K7 - alsopo two forms>?
      rps19 = "s" #U auth s, alpha subunit, P06508
    ),
    compartment_z = NULL,
    compartment_labels = NULL
  ),
  rpl_complex = list( # 
    pdb = "5H1S",  # Recent cryo-EM structure https://www.nature.com/articles/srep35793
    organism = "Spinacia oleracea",
    description = "large subunit of plastid ribosome",
    chain_map = list( #usually follow a nice SX pattern in the molecule name, map to alphabet
      rpl2 = "E", # S auth E P06509
      rpl14 = "M", #E auth m, #P09596
      rpl16 = "O", #G auth O P17353
      rpl20 = "S", # K auth S P28803
      rpl23 = "V", #N auth V Q9LWB5
      rpl33 = "c", # U auth c P28805 
      rpl36 = "f" # X auth F P12230 - alsopo two forms>?
    ),
    compartment_z = NULL,
    compartment_labels = NULL
  ),
  rp_complex = list( # 
    pdb = "5X8P",  # Recent cryo-EM structure https://www.nature.com/articles/srep35793
    organism = "Spinacia oleracea",
    description = "plastid ribosome",
    chain_map = list( #usually follow a nice SX pattern in the molecule name, map to alphabet
      rpl2 = "C", #  P06509
      rpl14 = "L", # P09596
      rpl16 = "N", # P17353
      rpl20 = "R", #  P28803
      rpl23 = "U", #T auth V Q9LWB5
      rpl33 = "2", #  P28805 
      rpl36 = "5", #  P12230 - 
      rps3 = "c", #E auth c, #P09595
      rps4 = "d", #F auth d P13788 # 
      rps7 = "g", # I auth G P82129, a and B?
      rps8 = "h", # J auth H P09597
      rps11 = "k", #M auth k P06506
      rps14 = "n", # P auth N P0656 # 
      rps18 = "r", # T auth r Q9M3K7 - alsopo two forms>?
      rps19 = "s" #U auth s, alpha subunit, P06508
    ),
    compartment_z = NULL,
    compartment_labels = NULL
  ),
  rnaPol_complex = list( # 
    pdb = "8WA1",  # Recent cryo-EM structure https://www.nature.com/articles/srep35793
    organism = "Spinacia oleracea",
    description = "plastid-encoded rna polymerase",
    chain_map = list( #usually follow a nice SX pattern in the molecule name, map to alphabet
      rpoA = c("A","a"), #  A, T, auth a P06269
      rpoB = "B", # P06271
      rpoC1 = "N", # A0A140G1Q3
    ),
    compartment_z = NULL,
    compartment_labels = NULL
  )
)

#' Standard amino acid 3-letter to 1-letter conversion
AA_MAP <- c(ALA="A", CYS="C", ASP="D", GLU="E", PHE="F", GLY="G", HIS="H", 
            ILE="I", LYS="K", LEU="L", MET="M", ASN="N", PRO="P", GLN="Q", 
            ARG="R", SER="S", THR="T", VAL="V", TRP="W", TYR="Y")


# ============================================================================
# STRUCTURE LOADING AND BASIC PROCESSING
# ============================================================================

#' Load PDB structure and extract key information
#' 
#' @param pdb_id PDB identifier (e.g., "6FKF") or path to local PDB file
#' @param cache_dir Directory to cache downloaded structures
#' @return List with: pdb (raw), ca_df (CA atoms), ligands (ligand centers), metadata
load_pdb_structure <- function(pdb_id, cache_dir = "data/pdb_cache") {
  
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  
  message("Loading PDB: ", pdb_id)
  
  if (file.exists(pdb_id)) {
    pdb <- read.pdb(pdb_id)
  } else {
    
    pdb_file <- file.path(cache_dir, paste0(pdb_id, ".pdb"))
    cif_file <- file.path(cache_dir, paste0(pdb_id, ".cif"))
    
    if (file.exists(pdb_file)) {
      pdb <- read.pdb(pdb_file)
    } else if (file.exists(cif_file)) {
      pdb <- read.cif(cif_file)
    } else {
      pdb <- tryCatch({
        p <- read.pdb(pdb_id)
        write.pdb(p, pdb_file)
        p
      }, error = function(e) {
        cif <- get.pdb(pdb_id, format = "cif")
        file.rename(cif, cif_file)
        read.pdb(cif_file)
      })
    }
  }
  
  ca_sel <- atom.select(pdb, elety = "CA")
  ca_atoms <- pdb$atom[ca_sel$atom, ]
  
  ca_df <- as.data.table(ca_atoms[, c("chain", "resno", "resid", "x", "y", "z")])
  setnames(ca_df, c("chain", "resno", "aa_3letter", "x", "y", "z"))
  ca_df[, aa := AA_MAP[aa_3letter]]
  
  if (pdb_id == "4XK8") ca_df <- ca_df[chain == toupper(chain)]
  
  protein_resids <- names(AA_MAP)
  water_resids <- c("HOH", "WAT", "H2O")
  ligand_atoms <- pdb$atom[!pdb$atom$resid %in% c(protein_resids, water_resids), ]
  
  ligands <- NULL
  if (nrow(ligand_atoms) > 0) {
    ligands <- as.data.table(ligand_atoms)[, .(
      lig_x = mean(x), lig_y = mean(y), lig_z = mean(z),
      n_atoms = .N
    ), by = .(chain, resid, resno)]
    setnames(ligands, "resid", "ligand_type")
  }
  
  chain_summary <- ca_df[, .(
    n_residues = .N,
    first_resno = min(resno),
    last_resno = max(resno)
  ), by = chain][order(-n_residues)]
  
  list(
    pdb = pdb,
    ca_df = ca_df,
    ligands = ligands,
    chain_summary = chain_summary,
    pdb_id = pdb_id
  )
}



#' Get sequence from a chain in the structure
#' 
#' @param ca_df CA atom data.table
#' @param chain_id Chain identifier
#' @return Character string of 1-letter amino acid sequence
get_chain_sequence <- function(ca_df, chain_id) {
  chain_data <- ca_df[chain == chain_id][order(resno)]
  paste(chain_data$aa, collapse = "")
}


# ============================================================================
# ALIGNMENT AND CHAIN MAPPING
# ============================================================================

#' Map PDB chains to genes via sequence alignment
#' 
#' @param structure Output from load_pdb_structure()
#' @param alignment_path Path to aligned FASTA file for the gene
#' @param reference_pattern Grep pattern to find reference sequence in alignment
#' @param gene_name Name of the gene (for output)
#' @param candidate_chains Character vector of chains to try (NULL = all chains)
#' @param min_identity Minimum percent identity to accept mapping
#' @param alignment_type "local" or "global" - local better for partial overlaps
#' @return List with: mapping (data.table), alignment diagnostics, best_chain
map_chain_to_alignment <- function(structure, 
                                   alignment_path,
                                   reference_pattern,
                                   gene_name,
                                   candidate_chains = NULL,
                                   min_identity = 50,
                                   alignment_type = "local") {
  
  # Load alignment
  aln <- readAAStringSet(alignment_path)
  
  # Find reference sequence
  ref_idx <- grep(reference_pattern, names(aln))
  stopifnot("Reference sequence not found in alignment" = length(ref_idx) >= 1)
  if (length(ref_idx) > 1) {
    warning("Multiple matches for reference pattern, using first")
    ref_idx <- ref_idx[1]
  }
  
  ref_aln_seq <- as.character(aln[[ref_idx]])
  ref_ungapped <- gsub("-", "", ref_aln_seq)
  
  message("Aligning ", gene_name, " to structure ", structure$pdb_id)
  message("  Reference: ", names(aln)[ref_idx])
  message("  Reference ungapped length: ", nchar(ref_ungapped))
  
  # Get candidate chains
  if (is.null(candidate_chains)) {
    candidate_chains <- unique(structure$ca_df$chain)
  }
  
  # Try aligning each candidate chain
  alignment_results <- rbindlist(lapply(candidate_chains, function(ch) {
    chain_seq <- get_chain_sequence(structure$ca_df, ch)
    if (nchar(chain_seq) < 20) return(NULL)  
    
    pw <- tryCatch(
      pairwiseAlignment(chain_seq, ref_ungapped, type = alignment_type),
      error = function(e) NULL
    )
    
    if (is.null(pw)) return(NULL)
    
    data.table(
      chain = ch,
      chain_length = nchar(chain_seq),
      ref_length = nchar(ref_ungapped),
      score = score(pw),
      pct_identity = pid(pw),
      n_matches = nmatch(pw),
      n_mismatches = nmismatch(pw)
    )
  }))
  
  if (nrow(alignment_results) == 0) {
    warning("No chains aligned successfully to ", gene_name)
    return(NULL)
  }
  
  alignment_results <- alignment_results[order(-pct_identity)]
  
  # Filter by minimum identity
  good_chains <- alignment_results[pct_identity >= min_identity]
  
  if (nrow(good_chains) == 0) {
    warning("No chains meet minimum identity threshold for ", gene_name)
    print(alignment_results)
    return(NULL)
  }
  
  message("  Best chains (>", min_identity, "% identity):")
  print(good_chains)
  
  # Build position mapping for the best chain
  best_chain <- good_chains$chain[1]
  chain_seq <- get_chain_sequence(structure$ca_df, best_chain)
  chain_ca <- structure$ca_df[chain == best_chain]
  
  # Build residue-level mapping using new function
  position_map <- build_position_mapping(
    chain_seq, ref_ungapped, ref_aln_seq, chain_ca, alignment_type
  )
  
  position_map[, gene := gene_name]

  list(
    gene = gene_name,
    alignment_results = alignment_results,
    good_chains = good_chains$chain,
    best_chain = best_chain,
    position_map = position_map,
    ref_name = names(aln)[ref_idx],
    ref_ungapped_length = nchar(ref_ungapped),
    alignment_path = alignment_path
  )
}


#' Build position mapping from pairwise alignment
#' 
#' @param chain_seq 1-letter sequence of PDB chain
#' @param ref_ungapped Reference sequence WITHOUT gaps
#' @param ref_aln_seq Reference sequence WITH gaps (from alignment)
#' @param chain_ca CA atoms for the chain (data.table with resno, x, y, z)
#' @param alignment_type "local" or "global" - local handles partial overlaps better
#' @return data.table mapping alignment positions to PDB residues
build_position_mapping <- function(chain_seq, ref_ungapped, ref_aln_seq, chain_ca, 
                                   alignment_type = "local") {
  
  # Perform pairwise alignment
  pw <- pairwiseAlignment(chain_seq, ref_ungapped, type = alignment_type)
  
  aln_pattern <- as.character(pattern(pw))  # PDB sequence (with gaps)
  aln_subject <- as.character(subject(pw))  # Reference (with gaps from pairwise)
  
  # For local alignment, we need the start positions
  pdb_start <- start(pattern(pw))
  ref_start <- start(subject(pw))
  
  # Track positions - start from the alignment start positions
  pdb_pos <- pdb_start - 1L
  ref_ungapped_pos <- ref_start - 1L
  
  chain_ca <- chain_ca[order(resno)]
  pdb_resnos <- chain_ca$resno
  
  pdb_to_ref <- data.table(
    pdb_resno = integer(),
    ref_ungapped_pos = integer()
  )
  
  for (i in seq_len(nchar(aln_pattern))) {
    p_char <- substr(aln_pattern, i, i)
    s_char <- substr(aln_subject, i, i)
    
    if (p_char != "-") pdb_pos <- pdb_pos + 1L
    if (s_char != "-") ref_ungapped_pos <- ref_ungapped_pos + 1L
    
    if (p_char != "-" && s_char != "-" && pdb_pos <= length(pdb_resnos)) {
      pdb_to_ref <- rbind(pdb_to_ref, data.table(
        pdb_resno = pdb_resnos[pdb_pos],
        ref_ungapped_pos = ref_ungapped_pos
      ))
    }
  }
  
  # Now map ref_ungapped_pos to alignment position
  aln_chars <- strsplit(ref_aln_seq, "")[[1]]
  ungapped_to_aln <- data.table(
    aln_pos = seq_along(aln_chars),
    aln_char = aln_chars,
    ref_ungapped_pos = NA_integer_
  )
  
  ungapped <- 0L
  for (i in seq_along(aln_chars)) {
    if (aln_chars[i] != "-") {
      ungapped <- ungapped + 1L
      ungapped_to_aln$ref_ungapped_pos[i] <- ungapped
    }
  }
  
  # Merge: PDB resno -> ref ungapped -> alignment position
  position_map <- merge(pdb_to_ref, 
                        ungapped_to_aln[!is.na(ref_ungapped_pos), .(aln_pos, ref_ungapped_pos)],
                        by = "ref_ungapped_pos", all.x = TRUE)
  
  # Add coordinates from chain_ca
  position_map <- merge(position_map, 
                        chain_ca[, .(pdb_resno = resno, x, y, z, aa)],
                        by = "pdb_resno", all.x = TRUE)
  
  # Return with alignment info as attributes
  result <- position_map[, .(aln_pos, pdb_resno, aa, x, y, z)]
  attr(result, "alignment") <- pw
  attr(result, "pdb_start") <- pdb_start
  attr(result, "ref_start") <- ref_start
  attr(result, "pct_identity") <- pid(pw)
  
  result
}


#' Map multiple chains for a gene (homomeric subunits)
#' 
#' @param structure Output from load_pdb_structure()
#' @param alignment_path Path to aligned FASTA
#' @param reference_pattern Pattern to find reference in alignment
#' @param gene_name Gene name
#' @param chains Character vector of chains for this gene
#' @param alignment_type "local" or "global" - local better for partial overlaps
#' @return data.table with all chains mapped, one row per (aln_pos, chain)
map_multiple_chains <- function(structure,
                                alignment_path, 
                                reference_pattern,
                                gene_name,
                                chains,
                                alignment_type = "local") {
  
  # Load alignment
  aln <- readAAStringSet(alignment_path)
  
  # Find reference sequence
  ref_idx <- grep(reference_pattern, names(aln))
  stopifnot("Reference sequence not found in alignment" = length(ref_idx) >= 1)
  if (length(ref_idx) > 1) {
    warning("Multiple matches for reference pattern, using first")
    ref_idx <- ref_idx[1]
  }
  
  ref_aln_seq <- as.character(aln[[ref_idx]])
  ref_ungapped <- gsub("-", "", ref_aln_seq)
  
  message("Mapping ", gene_name, " (", length(chains), " chains) to alignment")
  message("  Reference: ", substr(names(aln)[ref_idx], 1, 50))
  message("  Reference length: ", nchar(ref_ungapped), " aa")
  
  # Map each chain
  all_maps <- rbindlist(lapply(chains, function(ch) {
    chain_seq <- get_chain_sequence(structure$ca_df, ch)
    if (nchar(chain_seq) < 20) {
      message("  Chain ", ch, ": too short (", nchar(chain_seq), " aa), skipping")
      return(NULL)
    }
    
    chain_ca <- structure$ca_df[chain == ch]
    
    pos_map <- tryCatch({
      pm <- build_position_mapping(chain_seq, ref_ungapped, ref_aln_seq, 
                                   chain_ca, alignment_type)
      pct_id <- attr(pm, "pct_identity")
      message("  Chain ", ch, ": ", nrow(pm), " positions mapped, ", 
              round(pct_id, 1), "% identity")
      pm
    }, error = function(e) {
      message("  Chain ", ch, ": mapping failed - ", e$message)
      NULL
    })
    
    if (is.null(pos_map) || nrow(pos_map) == 0) return(NULL)
    
    pos_map[, chain := ch]
    pos_map[, gene := gene_name]
    pos_map
  }))
  
  if (is.null(all_maps) || nrow(all_maps) == 0) {
    warning("No chains mapped successfully for ", gene_name)
    return(NULL)
  }
  
  # Add attributes
  attr(all_maps, "gene") <- gene_name
  attr(all_maps, "chains") <- chains
  attr(all_maps, "ref_name") <- names(aln)[ref_idx]
  
  all_maps
}


# ============================================================================
# STRUCTURAL FEATURE COMPUTATION
# ============================================================================

#' Compute structural features for mapped residues
#' 
#' @param position_maps data.table from map_multiple_chains (all chains)
#' @param structure Output from load_pdb_structure()
#' @param complex_def Complex definition with compartment info
#' @param other_chains Character vector of chains from OTHER subunits (for interface calc)
#' @return data.table with one row per alignment position, features aggregated across chains
compute_residue_features <- function(position_maps,
                                     structure,
                                     complex_def = NULL,
                                     other_chains = NULL) {
  
  gene_name <- attr(position_maps, "gene")
  if (is.null(gene_name)) gene_name <- position_maps$gene[1]
  
  ca_df <- structure$ca_df
  ligands <- structure$ligands
  
  # Get coordinates of other subunits for interface calculation
  if (is.null(other_chains)) {
    # Default: all chains not belonging to this gene
    my_chains <- unique(position_maps$chain)
    other_chains <- setdiff(unique(ca_df$chain), my_chains)
  }
  other_coords <- ca_df[chain %in% other_chains, .(x, y, z)]
  
  # Ligand coordinates (if any)
  lig_coords <- NULL
  if (!is.null(ligands) && nrow(ligands) > 0) {
    lig_coords <- ligands[, .(x = lig_x, y = lig_y, z = lig_z)]
  }
  
  # Compute features per (aln_pos, chain)
  position_maps[, `:=`(
    # Distance to nearest ligand
    dist_to_ligand = if (!is.null(lig_coords) && nrow(lig_coords) > 0) {
      sapply(seq_len(.N), function(i) {
        min(sqrt((x[i] - lig_coords$x)^2 + 
                   (y[i] - lig_coords$y)^2 + 
                   (z[i] - lig_coords$z)^2))
      })
    } else NA_real_,
    
    # Distance to nearest other-subunit residue (interface)
    dist_to_interface = if (nrow(other_coords) > 0) {
      sapply(seq_len(.N), function(i) {
        min(sqrt((x[i] - other_coords$x)^2 + 
                   (y[i] - other_coords$y)^2 + 
                   (z[i] - other_coords$z)^2))
      })
    } else NA_real_
  )]
  
  # Aggregate across chains (for homomeric subunits)
  features <- position_maps[, .(
    n_chains = .N,
    
    # Coordinates: mean position
    x_mean = mean(x, na.rm = TRUE),
    y_mean = mean(y, na.rm = TRUE),
    z_mean = mean(z, na.rm = TRUE),
    
    # Conformational flexibility: how much does this residue move across chains?
    conformational_range = sqrt(var(x, na.rm = TRUE) + var(y, na.rm = TRUE) + var(z, na.rm = TRUE)),
    x_range = diff(range(x, na.rm = TRUE)),
    y_range = diff(range(y, na.rm = TRUE)),
    z_range = diff(range(z, na.rm = TRUE)),
    
    # Ligand distance: use MINIMUM across chains (captures "can this ever be near ligand?")
    dist_ligand_min = min(dist_to_ligand, na.rm = TRUE),
    dist_ligand_mean = mean(dist_to_ligand, na.rm = TRUE),
    dist_ligand_max = max(dist_to_ligand, na.rm = TRUE),
    
    # Interface distance: use MINIMUM (captures "is this at any interface?")
    dist_interface_min = min(dist_to_interface, na.rm = TRUE),
    dist_interface_mean = mean(dist_to_interface, na.rm = TRUE)
    
  ), by = .(gene, aln_pos)]
  
  # Handle Inf values from empty ranges
  for (col in names(features)) {
    if (is.numeric(features[[col]])) {
      features[is.infinite(get(col)), (col) := NA_real_]
    }
  }
  
  # Add compartment if defined
  if (!is.null(complex_def$compartment_z)) {
    z_bounds <- complex_def$compartment_z
    z_labels <- complex_def$compartment_labels
    
    features[, compartment := fcase(
      z_mean < z_bounds[1], z_labels[1],
      z_mean >= z_bounds[1] & z_mean < z_bounds[2], z_labels[2],
      z_mean >= z_bounds[2], z_labels[3],
      default = NA_character_
    )]
  }
  
  # Add "at interface" binary
  features[, at_interface := dist_interface_min < 8]
  
  # Add "near ligand" binary
  features[, near_ligand := dist_ligand_min < 15]
  
  features
}


#' Compute distance from central axis (for rotary enzymes like ATP synthase)
#' 
#' @param features Output from compute_residue_features()
#' @param structure Output from load_pdb_structure()
#' @param axis_chains Chains defining the axis (e.g., c-ring for ATP synthase)
#' @return features with dist_to_axis added
add_axis_distance <- function(features, structure, axis_chains) {
  
  # Calculate center of axis-defining chains in XY plane
  axis_center <- structure$ca_df[chain %in% axis_chains, 
                                 .(cx = mean(x), cy = mean(y))]
  
  features[, dist_to_axis := sqrt((x_mean - axis_center$cx)^2 + 
                                    (y_mean - axis_center$cy)^2)]
  features
}


# ============================================================================
# GWAS INTEGRATION AND STATISTICAL TESTS
# ============================================================================

#' Merge GWAS results with structural features
#' 
#' @param gwas_sites data.table with Gene, Position (=aln_pos), P_aa_only, P_aa_with_pcs
#' @param features Output from compute_residue_features()
#' @param thresh_aa_only P-value threshold for sig_no_control classification
#' @param thresh_aa_pcs P-value threshold for sig_with_control classification
#' @return Merged data.table with site_class column
merge_gwas_with_structure <- function(gwas_sites, 
                                      features,
                                      thresh_aa_only,
                                      thresh_aa_pcs) {
  
  gene_name <- unique(features$gene)
  stopifnot(length(gene_name) == 1)
  
  # Filter GWAS to this gene
  gwas_gene <- gwas_sites[Gene == gene_name]
  
  # Merge
  merged <- merge(gwas_gene, features, 
                  by.x = c("Gene", "Position"), 
                  by.y = c("gene", "aln_pos"),
                  all.x = TRUE)
  
  # Classify sites
  merged[, site_class := "not_sig"]
  merged[P_aa_only < thresh_aa_only, site_class := "sig_no_control"]
  merged[P_aa_with_pcs < thresh_aa_pcs, site_class := "sig_with_control"]
  merged[P_aa_only < thresh_aa_only & P_aa_with_pcs < thresh_aa_pcs, 
         site_class := "sig_both"]
  
  merged[, site_class := factor(site_class, 
                                levels = c("not_sig", "sig_no_control", 
                                           "sig_with_control", "sig_both"))]
  
  # Add -log10 p-values
  merged[, neglog_p_only := -log10(P_aa_only)]
  merged[, neglog_p_pcs := -log10(P_aa_with_pcs)]
  
  message("Merged ", gene_name, ": ", nrow(merged), " sites, ",
          sum(!is.na(merged$x_mean)), " with structure")
  message("  Site classes: ", paste(names(table(merged$site_class)), 
                                    table(merged$site_class), 
                                    sep = "=", collapse = ", "))
  
  merged
}


#' Run enrichment tests for a single continuous feature
#' 
#' @param data Merged GWAS + structure data
#' @param feature_col Name of feature column
#' @param direction "greater" if sig sites expected to have larger values, "less" otherwise
#' @return data.table with test results
test_feature_enrichment <- function(data, 
                                    feature_col, 
                                    direction = "two.sided") {
  
  # Filter to sites with this feature
  data_filt <- data[!is.na(get(feature_col))]
  bg <- data_filt[site_class == "not_sig"]
  
  if (nrow(bg) < 5) {
    warning("Insufficient background sites for ", feature_col)
    return(NULL)
  }
  
  results <- rbindlist(lapply(c("sig_no_control", "sig_with_control", "sig_both"), 
                              function(sc) {
                                test_data <- data_filt[site_class == sc]
                                
                                if (nrow(test_data) < 3) {
                                  return(data.table(
                                    feature = feature_col,
                                    site_class = sc,
                                    n_sig = nrow(test_data),
                                    n_bg = nrow(bg),
                                    median_sig = NA_real_,
                                    median_bg = median(bg[[feature_col]], na.rm = TRUE),
                                    wilcox_p = NA_real_,
                                    direction = NA_character_
                                  ))
                                }
                                
                                # Wilcoxon test
                                wt <- wilcox.test(test_data[[feature_col]], bg[[feature_col]], 
                                                  alternative = direction)
                                
                                med_sig <- median(test_data[[feature_col]], na.rm = TRUE)
                                med_bg <- median(bg[[feature_col]], na.rm = TRUE)
                                
                                data.table(
                                  feature = feature_col,
                                  site_class = sc,
                                  n_sig = nrow(test_data),
                                  n_bg = nrow(bg),
                                  median_sig = med_sig,
                                  median_bg = med_bg,
                                  wilcox_p = wt$p.value,
                                  direction = ifelse(med_sig > med_bg, "higher", "lower")
                                )
                              }))
  
  # Add correlation (across all sites, not just sig)
  cor_result <- cor.test(data_filt[[feature_col]], data_filt$neglog_p_only, 
                         method = "spearman")
  
  results[, `:=`(
    spearman_rho = cor_result$estimate,
    spearman_p = cor_result$p.value
  )]
  
  results
}


#' Run enrichment test for a binary feature (e.g., at_interface)
#' 
#' @param data Merged GWAS + structure data  
#' @param feature_col Name of logical/binary feature column
#' @param feature_label Human-readable label
#' @return data.table with Fisher's exact test results
test_binary_enrichment <- function(data,
                                   feature_col,
                                   feature_label = NULL) {
  
  if (is.null(feature_label)) feature_label <- feature_col
  
  data_filt <- data[!is.na(get(feature_col))]
  
  results <- rbindlist(lapply(c("sig_no_control", "sig_with_control", "sig_both"),
                              function(sc) {
                                sig_in <- sum(data_filt$site_class == sc & data_filt[[feature_col]] == TRUE, na.rm = TRUE)
                                sig_out <- sum(data_filt$site_class == sc & data_filt[[feature_col]] == FALSE, na.rm = TRUE)
                                bg_in <- sum(data_filt$site_class == "not_sig" & data_filt[[feature_col]] == TRUE, na.rm = TRUE)
                                bg_out <- sum(data_filt$site_class == "not_sig" & data_filt[[feature_col]] == FALSE, na.rm = TRUE)
                                
                                if (sig_in + sig_out == 0) {
                                  return(data.table(
                                    feature = feature_label,
                                    site_class = sc,
                                    n_sig_in = sig_in,
                                    n_sig_total = 0,
                                    pct_sig = NA_real_,
                                    pct_bg = 100 * bg_in / max(bg_in + bg_out, 1),
                                    odds_ratio = NA_real_,
                                    fisher_p = NA_real_
                                  ))
                                }
                                
                                mat <- matrix(c(sig_in, sig_out, bg_in, bg_out), nrow = 2)
                                ft <- fisher.test(mat)
                                
                                data.table(
                                  feature = feature_label,
                                  site_class = sc,
                                  n_sig_in = sig_in,
                                  n_sig_total = sig_in + sig_out,
                                  pct_sig = 100 * sig_in / (sig_in + sig_out),
                                  pct_bg = 100 * bg_in / (bg_in + bg_out),
                                  odds_ratio = ft$estimate,
                                  fisher_p = ft$p.value
                                )
                              }))
  
  results
}


#' Run all standard structural enrichment tests
#' 
#' @param merged_data Output from merge_gwas_with_structure()
#' @return List with continuous_tests and binary_tests data.tables
run_all_enrichment_tests <- function(merged_data) {
  
  gene_name <- merged_data$Gene[1]
  message("\n=== Running enrichment tests for ", gene_name, " ===\n")
  
  # Continuous features to test
  cont_features <- c("z_mean", "dist_to_axis", "dist_ligand_min", 
                     "dist_interface_min", "conformational_range")
  cont_features <- cont_features[cont_features %in% names(merged_data)]
  
  continuous_results <- rbindlist(lapply(cont_features, function(f) {
    test_feature_enrichment(merged_data, f)
  }))
  
  # Binary features
  bin_features <- c("at_interface", "near_ligand")
  bin_features <- bin_features[bin_features %in% names(merged_data)]
  
  binary_results <- rbindlist(lapply(bin_features, function(f) {
    test_binary_enrichment(merged_data, f)
  }))
  
  # Compartment enrichment (if present)
  compartment_results <- NULL
  if ("compartment" %in% names(merged_data)) {
    compartments <- unique(merged_data$compartment[!is.na(merged_data$compartment)])
    
    compartment_results <- rbindlist(lapply(compartments, function(comp) {
      merged_data[, temp_comp := compartment == comp]
      res <- test_binary_enrichment(merged_data, "temp_comp", paste0("compartment:", comp))
      merged_data[, temp_comp := NULL]
      res
    }))
  }
  
  list(
    gene = gene_name,
    continuous = continuous_results,
    binary = rbind(binary_results, compartment_results, fill = TRUE)
  )
}


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

#' Plot structure with GWAS significance overlay
#' 
#' @param merged_data Output from merge_gwas_with_structure()
#' @param structure Output from load_pdb_structure()
#' @param position_maps Output from map_multiple_chains() - needed for per-chain coords
#' @param projection "xy", "xz", or "yz"
#' @param title Plot title
#' @param show_ligands Show ligand positions as diamonds
#' @param highlight_chains Character vector of chains to highlight (NULL = use all from position_maps)
#' @return ggplot object
plot_gwas_on_structure <- function(merged_data, 
                                   structure,
                                   position_maps,
                                   projection = "xy",
                                   title = NULL,
                                   show_ligands = TRUE,
                                   highlight_chains = NULL) {
  
  ca_df <- structure$ca_df
  ligands <- structure$ligands
  
  # Get projection coordinates
  if (projection == "xy") {
    coord1 <- "x"; coord2 <- "y"
  } else if (projection == "xz") {
    coord1 <- "x"; coord2 <- "z"
  } else {
    coord1 <- "y"; coord2 <- "z"
  }
  
  if (is.null(title)) {
    title <- paste0(merged_data$Gene[1], " GWAS hits - ", toupper(projection), " projection")
  }
  
  # Use position_maps (per-chain coordinates) for plotting
  # Merge site_class from merged_data onto position_maps
  plot_coords <- merge(
    position_maps[, .(gene, aln_pos, chain, x, y, z)],
    merged_data[, .(Gene, Position, site_class)],
    by.x = c("gene", "aln_pos"),
    by.y = c("Gene", "Position"),
    all.x = TRUE
  )
  plot_coords[is.na(site_class), site_class := "no_gwas_data"]
  
  # Filter chains if specified
  if (!is.null(highlight_chains)) {
    plot_coords <- plot_coords[chain %in% highlight_chains]
  }
  
  # Get chains being plotted for subtitle
  chains_plotted <- unique(plot_coords$chain)
  
  p <- ggplot() +
    # Background: all CA atoms (gray)
    geom_point(data = ca_df, 
               aes_string(x = coord1, y = coord2),
               color = "gray70", size = 0.3, alpha = 0.3) +
    # Not significant (on actual chain positions)
    geom_point(data = plot_coords[site_class == "not_sig"],
               aes_string(x = coord1, y = coord2),
               color = "gray50", size = 1.5, alpha = 0.5) +
    # Significant classes
    geom_point(data = plot_coords[site_class == "sig_no_control"],
               aes_string(x = coord1, y = coord2),
               color = "gold", size = 2.5) +
    geom_point(data = plot_coords[site_class == "sig_with_control"],
               aes_string(x = coord1, y = coord2),
               color = "steelblue", size = 2.5) +
    geom_point(data = plot_coords[site_class == "sig_both"],
               aes_string(x = coord1, y = coord2),
               color = "darkred", size = 3) +
    coord_fixed() +
    labs(title = title,
         subtitle = paste0("Chains: ", paste(chains_plotted, collapse = ","), 
                           " | Gold=sig_no_ctrl, Blue=sig_with_ctrl, Red=sig_both"),
         x = paste0(toupper(coord1), " (Å)"),
         y = paste0(toupper(coord2), " (Å)")) +
    theme_classic()
  
  # Add ligands if requested
  if (show_ligands && !is.null(ligands) && nrow(ligands) > 0) {
    lig_coord1 <- ifelse(coord1 == "x", "lig_x", ifelse(coord1 == "y", "lig_y", "lig_z"))
    lig_coord2 <- ifelse(coord2 == "x", "lig_x", ifelse(coord2 == "y", "lig_y", "lig_z"))
    
    p <- p + geom_point(data = ligands,
                        aes_string(x = lig_coord1, y = lig_coord2),
                        color = "black", size = 4, shape = 18)
  }
  
  p
}


#' Plot boxplots of feature by site class
#' 
#' @param merged_data Output from merge_gwas_with_structure()
#' @param feature_col Feature column name
#' @param ylabel Y-axis label
#' @param title Plot title
#' @return ggplot object
plot_feature_boxplot <- function(merged_data,
                                 feature_col,
                                 ylabel = NULL,
                                 title = NULL) {
  
  if (is.null(ylabel)) ylabel <- feature_col
  if (is.null(title)) title <- paste(merged_data$Gene[1], ":", feature_col)
  
  plot_data <- merged_data[!is.na(get(feature_col))]
  
  ggplot(plot_data, aes_string(x = "site_class", y = feature_col, fill = "site_class")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("not_sig" = "gray60", 
                                 "sig_no_control" = "gold",
                                 "sig_with_control" = "steelblue", 
                                 "sig_both" = "darkred")) +
    labs(title = title, y = ylabel, x = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}


#' Plot combined GWAS hits for entire complex (all genes, all chains)
#' 
#' @param complex_results Output from analyze_complex()
#' @param projection "xy", "xz", or "yz"
#' @param title Plot title
#' @param show_ligands Show ligand positions
#' @param compartment_lines Y-values for horizontal compartment boundary lines (for xz/yz projections)
#' @return ggplot object
plot_complex_gwas <- function(complex_results,
                              projection = "xz",
                              title = NULL,
                              show_ligands = TRUE,
                              compartment_lines = NULL) {
  
  structure <- complex_results$structure
  ca_df <- structure$ca_df
  ligands <- structure$ligands
  
  # Get projection coordinates
  if (projection == "xy") {
    coord1 <- "x"; coord2 <- "y"
  } else if (projection == "xz") {
    coord1 <- "x"; coord2 <- "z"
  } else {
    coord1 <- "y"; coord2 <- "z"
  }
  
  if (is.null(title)) {
    title <- paste0(complex_results$complex_name, " GWAS hits - ", toupper(projection), " projection")
  }
  
  # Combine position_maps from all genes, merge with site_class
  all_coords <- rbindlist(lapply(complex_results$gene_results, function(r) {
    if (is.null(r$position_maps)) return(NULL)
    
    # Merge site_class onto position_maps
    coords <- merge(
      r$position_maps[, .(gene, aln_pos, chain, x, y, z)],
      r$merged[, .(Gene, Position, site_class)],
      by.x = c("gene", "aln_pos"),
      by.y = c("Gene", "Position"),
      all.x = TRUE
    )
    coords[is.na(site_class), site_class := "no_gwas_data"]
    coords
  }), fill = TRUE)
  
  if (is.null(all_coords) || nrow(all_coords) == 0) {
    stop("No coordinate data found in complex_results")
  }
  
  # Use compartment lines from complex_def if not provided
  if (is.null(compartment_lines) && !is.null(complex_results$complex_def$compartment_z)) {
    compartment_lines <- complex_results$complex_def$compartment_z
  }
  
  p <- ggplot() +
    # Background: all CA atoms
    geom_point(data = ca_df, 
               aes_string(x = coord1, y = coord2),
               color = "gray70", size = 0.3, alpha = 0.3) +
    # Not significant
    geom_point(data = all_coords[site_class == "not_sig"],
               aes_string(x = coord1, y = coord2),
               color = "gray50", size = 1, alpha = 0.5) +
    # Significant classes
    geom_point(data = all_coords[site_class == "sig_no_control"],
               aes_string(x = coord1, y = coord2),
               color = "gold", size = 2) +
    geom_point(data = all_coords[site_class == "sig_with_control"],
               aes_string(x = coord1, y = coord2),
               color = "steelblue", size = 2) +
    geom_point(data = all_coords[site_class == "sig_both"],
               aes_string(x = coord1, y = coord2),
               color = "darkred", size = 2.5) +
    coord_fixed() +
    labs(title = title,
         subtitle = "Gold=sig_no_ctrl, Blue=sig_with_ctrl, Red=sig_both",
         x = paste0(toupper(coord1), " (Å)"),
         y = paste0(toupper(coord2), " (Å)")) +
    theme_classic()
  
  # Add compartment boundary lines (for xz or yz projections)
  if (!is.null(compartment_lines) && projection %in% c("xz", "yz")) {
    p <- p + geom_hline(yintercept = compartment_lines, linetype = "dashed", alpha = 0.5)
  }
  
  # Add ligands if requested
  if (show_ligands && !is.null(ligands) && nrow(ligands) > 0) {
    lig_coord1 <- ifelse(coord1 == "x", "lig_x", ifelse(coord1 == "y", "lig_y", "lig_z"))
    lig_coord2 <- ifelse(coord2 == "x", "lig_x", ifelse(coord2 == "y", "lig_y", "lig_z"))
    
    p <- p + geom_point(data = ligands,
                        aes_string(x = lig_coord1, y = lig_coord2),
                        color = "black", size = 4, shape = 18)
  }
  
  p
}


#' Generate diagnostic alignment plot
#' 
#' @param chain_map_result Output from map_chain_to_alignment()
#' @return ggplot showing alignment quality across chains
plot_alignment_diagnostics <- function(chain_map_result) {
  
  aln_results <- chain_map_result$alignment_results
  
  p <- ggplot(aln_results, aes(x = reorder(chain, -pct_identity), y = pct_identity)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    labs(title = paste(chain_map_result$gene, "- Chain alignment quality"),
         subtitle = paste("Reference:", chain_map_result$ref_name),
         x = "Chain", y = "Percent Identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p
}


# ============================================================================
# HIGH-LEVEL WORKFLOW FUNCTIONS
# ============================================================================

#' Analyze a single gene within a complex
#' 
#' @param gene_name Gene name (e.g., "atpB")
#' @param structure Output from load_pdb_structure()
#' @param complex_def Complex definition from DEFAULT_COMPLEXES
#' @param alignment_dir Directory containing gene alignments
#' @param reference_pattern Grep pattern for reference species
#' @param gwas_sites GWAS results data.table
#' @param thresh_aa_only Threshold for P_aa_only
#' @param thresh_aa_pcs Threshold for P_aa_with_pcs
#' @return List with all analysis results
analyze_gene <- function(gene_name,
                         structure,
                         complex_def,
                         alignment_dir,
                         reference_pattern,
                         gwas_sites,
                         thresh_aa_only,
                         thresh_aa_pcs) {
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Analyzing: ", gene_name)
  message(paste(rep("=", 60), collapse = ""))
  
  # Get chains for this gene
  chains <- complex_def$chain_map[[gene_name]]
  if (is.null(chains)) {
    warning("No chain mapping for ", gene_name)
    return(NULL)
  }
  
  # Alignment path
  aln_path <- file.path(alignment_dir, paste0(gene_name, "_AA_aligned.fasta"))
  if (!file.exists(aln_path)) {
    warning("Alignment not found: ", aln_path)
    return(NULL)
  }
  
  # Map chains to alignment
  position_maps <- map_multiple_chains(
    structure, aln_path, reference_pattern, gene_name, chains
  )
  
  if (is.null(position_maps) || nrow(position_maps) == 0) {
    warning("Failed to map chains for ", gene_name)
    return(NULL)
  }
  
  # Compute structural features
  other_chains <- setdiff(unique(structure$ca_df$chain), chains)
  features <- compute_residue_features(
    position_maps, structure, complex_def, other_chains
  )
  
  # Add axis distance if relevant (rotary enzymes)
  if (!is.null(complex_def$chain_map$atpH)) {
    # ATP synthase: use c-ring as axis
    features <- add_axis_distance(features, structure, complex_def$chain_map$atpH)
  }
  
  # Merge with GWAS
  merged <- merge_gwas_with_structure(
    gwas_sites, features, thresh_aa_only, thresh_aa_pcs
  )
  
  # Run tests
  tests <- run_all_enrichment_tests(merged)
  
  list(
    gene = gene_name,
    chains = chains,
    position_maps = position_maps,
    features = features,
    merged = merged,
    tests = tests
  )
}


#' Analyze entire complex
#' 
#' @param complex_name Name in DEFAULT_COMPLEXES
#' @param alignment_dir Directory with gene alignments
#' @param reference_pattern Reference species pattern
#' @param gwas_sites GWAS data
#' @param thresh_aa_only P-value threshold
#' @param thresh_aa_pcs P-value threshold with PC control
#' @param genes_to_analyze Subset of genes (NULL = all)
#' @return List of results per gene
analyze_complex <- function(complex_name,
                            alignment_dir,
                            reference_pattern,
                            gwas_sites,
                            thresh_aa_only,
                            thresh_aa_pcs,
                            genes_to_analyze = NULL) {
  
  complex_def <- DEFAULT_COMPLEXES[[complex_name]]
  if (is.null(complex_def)) {
    stop("Unknown complex: ", complex_name)
  }
  
  message("\n", paste(rep("#", 70), collapse = ""))
  message("ANALYZING COMPLEX: ", complex_name)
  message("  PDB: ", complex_def$pdb)
  message("  Organism: ", complex_def$organism)
  message(paste(rep("#", 70), collapse = ""))
  
  # Load structure
  structure <- load_pdb_structure(complex_def$pdb)
  
  # Determine genes to analyze
  if (is.null(genes_to_analyze)) {
    genes_to_analyze <- names(complex_def$chain_map)
  }
  
  # Analyze each gene
  results <- lapply(genes_to_analyze, function(g) {
    tryCatch(
      analyze_gene(g, structure, complex_def, alignment_dir,
                   reference_pattern, gwas_sites, 
                   thresh_aa_only, thresh_aa_pcs),
      error = function(e) {
        warning("Error analyzing ", g, ": ", e$message)
        NULL
      }
    )
  })
  names(results) <- genes_to_analyze
  
  # Remove failed analyses
  results <- results[!sapply(results, is.null)]
  
  # Combine merged data across genes
  combined_merged <- rbindlist(lapply(results, function(r) r$merged), fill = TRUE)
  
  list(
    complex_name = complex_name,
    complex_def = complex_def,
    structure = structure,
    gene_results = results,
    combined_merged = combined_merged
  )
}


# ============================================================================
# SUMMARY AND REPORTING
# ============================================================================

#' Print summary of enrichment tests
#' 
#' @param test_results Output from run_all_enrichment_tests()
print_test_summary <- function(test_results) {
  
  message("\n=== ENRICHMENT TEST SUMMARY: ", test_results$gene, " ===\n")
  
  if (!is.null(test_results$continuous) && nrow(test_results$continuous) > 0) {
    message("CONTINUOUS FEATURES (Wilcoxon vs not_sig):")
    print(test_results$continuous[, .(feature, site_class, n_sig, 
                                      median_sig, median_bg, 
                                      direction, wilcox_p)])
  }
  
  if (!is.null(test_results$binary) && nrow(test_results$binary) > 0) {
    message("\nBINARY FEATURES (Fisher's exact):")
    print(test_results$binary[, .(feature, site_class, n_sig_in, n_sig_total,
                                  pct_sig, pct_bg, odds_ratio, fisher_p)])
  }
}


#' Combine test results across genes
#' 
#' @param complex_results Output from analyze_complex()
#' @return List with combined continuous and binary test tables
combine_test_results <- function(complex_results) {
  
  continuous <- rbindlist(lapply(complex_results$gene_results, function(r) {
    if (!is.null(r$tests$continuous)) {
      r$tests$continuous[, gene := r$gene]
    }
  }), fill = TRUE)
  
  binary <- rbindlist(lapply(complex_results$gene_results, function(r) {
    if (!is.null(r$tests$binary)) {
      r$tests$binary[, gene := r$gene]
    }
  }), fill = TRUE)
  
  list(continuous = continuous, binary = binary)
}

#' Compute membrane-relative coordinate
#' 
#' Projects coordinates onto membrane normal axis.
#' Requires specifying reference points for lumen and stroma sides.
#' 
#' @param coords data.table with x, y, z columns
#' @param lumen_ref XYZ of a known lumenal residue/ligand
#' @param stroma_ref XYZ of a known stromal residue/ligand
#' @return Vector of membrane-relative positions (negative = lumen, positive = stroma)
compute_membrane_depth <- function(coords, lumen_ref, stroma_ref) {
  
  # Membrane normal vector (lumen -> stroma)
  normal <- c(stroma_ref[1] - lumen_ref[1],
              stroma_ref[2] - lumen_ref[2],
              stroma_ref[3] - lumen_ref[3])
  normal <- normal / sqrt(sum(normal^2))  # normalize
  
  
  # Membrane center
  center <- (lumen_ref + stroma_ref) / 2
  
  # Project each point onto normal
  depths <- sapply(1:nrow(coords), function(i) {
    vec <- c(coords$x[i] - center[1],
             coords$y[i] - center[2],
             coords$z[i] - center[3])
    sum(vec * normal)
  })
  
  depths
}
