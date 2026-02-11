# Chloroplast Genome Temperature Association

**Authors:** Henry Dawson, Sheng-kai Hsu, Edward S. Buckler

**Affiliations:** Section of Plant Breeding and Genetics, Cornell University, Ithaca, NY, USA 14853,  Institute for Genomic Diversity, Cornell University, Ithaca, NY, USA 14853, Agricultural Research Service, United States Department of Agriculture, Ithaca, NY, USA 14853

**Corresponding author:** hdd29@cornell.edu

---

## Abstract
TODO
<!-- ~250 words. Background, methods, key findings, conclusion. Write last. -->

## Introduction 

Plants generally struggle to grow in the cold. This affects planting seasons for many crops, and challenges yield, food security, and agricultural efficiency. A key component of cold stress in plants is photosynthetic. Thermodependent carbon sink processes slow while light continues to irradiate photosynthetic machinery. The chloroplast electron transport chain becomes overreduced and produces reactive oxygen species. The extent of this stress is such that a plant may inhibit photosynthesis under cold, unable to grow despite abundant light. 

The evolutionary history of the chloroplast offers a potential explanation for the mismatch between photosynthesis and the cold. According to the endosymbiotic theory, the chloroplast’s ancestor was a cyanobacterium that lived in the ocean ~1.5 billion years ago. The environment at that time was aquatic, likely warm, and slightly dimmer - a contrast to the dry, cold, bright conditions common in cold stress. The conservation of the chloroplast genome from angiosperms to algae suggests that the general mechanism of photosynthesis is retained from an ancestral state likely optimized for warm conditions. 
Despite strong functional constraints on variation, and a limited pool of processes that drive genomic variation, chloroplast-encoded proteins display variation across the plant kingdom. To what degree, if any, is this variation associated with plant growth temperature? 

This question is motivated by reports linking plastid genetic variation to temperature-dependent phenotypes and adaptive evolution in plastid genomes, and more broadly by the universal temperature dependence of enzymatic reaction rates, evidence for protein and proteome-level thermal adaptation in bacteria, and the well-documented sensitivity of the photosynthetic apparatus to cold stress.

Here we test the hypothesis that variation in chloroplast-encoded proteins is associated with temperature. Previous chloroplast genomic research focuses on phylogeny resolution and population genomic models of adaptive evolution, with a particular focus on rubisco. Here we analyze a general set of proteins spanning multiple complexes and pathways, and do so specifically with regards to plant growth temperature as estimated from occurrence records and climate models. 

We detected a significant convergent association, characterized by hydrophobicity changes and conservative substitutions localized near protein active sites, consistent with bacterial temperature adaptation and evolution under high functional constraint. This research offers a plastome-temperature genotype-phenotype map to guide engineering efforts, and the basis of a genomic prediction and selection tool to increase abiotic stress tolerance in plant breeding programs. 


## Methods

### Genomic Data Preparation

Chloroplast genomes were downloaded from NCBI (see [../code/1_download.py](../code/1_download.py)). Duplicates per species were resolved by ranking on gene count, number of ambiguous bases, recency, sequencing technology, and ties were resolved by choosing the genome with 2-mer frequency closest to the overall average (see [../code/dataProcessing.ipynb](../code/dataProcessing.ipynb)). Genomes were [annotated](../code/dataProcessing.ipynb) using CpGAVAS2. Genomes with outlying lengths (<50kbp, >300kbp) or gene counts (<80) were removed. Genomes belonging to parasitic or semi-parasitic plants Cuscuta x., Thesium chinense, Cassytha filiformis, Striga asiatica, Orobanche x., were removed. Samples belonging to orders with less than 50 species were removed. 

### Phenotype: 

Occurrences for each species were found through [GBIF](../code/XXX_processGBIFData.py) and BIEN. Occurrences were subjected to a [quality control pipeline](../code/XXX_CleanOccurrences.R) that filtered marine coordinates, nonsensical coordinates (e.g. 0,0), duplicate coordinates, invalid coordinates according to r_CoordinateCleaner, outlier coordinates according to r_CoordinateCleaner::cc_outl, removing occurrences at sites with a noisy measurement according to Topt-site, and keeping only species with at least 5 remaining occurrences. Temperature estimates were found by querying remaining occurrences into [worldclim BIO8](../code/XXX_EnvDataFromTifs.r), which represents the mean temperature of the wettest quarter, and taking the [median value per species](../code/XXX_processEnvData.r).


### Population structure PCS:

The genotype matrix was made by [aligning](../code/aligner.sh) the CDS of the 95% conserved gene set using MAFFT, then [processing](../code/XXX_makeAApcs.py) to drop sites with >5% gaps, and encoding the major variant at each site as 0 and all other variants as 1. Principal component analysis was [performed](../code/XXX_makeGenealogPCs.r) on the genotype matrix after median imputing the missing data.

### Residue-Temperature Association Analysis

[Tested](../code/residueGWAS_PopStepAndaJump.r) for associations between amino acid variation and plant growth temperature (species median BIO8) using position-wise linear models across 60 chloroplast genes. For each aligned position with >=8500 non-gap sequences and polymorphic residues (>= 2 alleles), we fit three nested models: (1) BIO8 ~ 1 (null), (2) BIO8 ~ residue, and (3) BIO8 ~ residue + population structure PCs. Population structure was controlled using the first 1000 principal components derived from major/minor allele-encoded sequence alignments. Significance was assessed by ANOVA comparing model fits, yielding P-values for residue effects without population control (P_aa_only) and with control (P_aa_with_pcs). Sites were classified as significant if P_aa_with_pcs fell below the 5th percentile threshold. 

### Phylogenetic reconstruction

To construct a species tree, we built a superalignment from the 60 chloroplast genes passing quality control (≥8500 non-gap sequences per site, ≥2 allelic variants). To reduce computational burden while retaining phylogenetic diversity, we [collapsed](../code/tipCollapsingAA.r) congeneric taxa within 30 amino acid substitutions of each other, selecting one representative per cluster (n = 6,094 taxa retained from ~10,857). The resulting 11,427-site amino acid superalignment was analyzed using RAxML-NG v1.2.0 with the LG+G4 substitution model. Tree search was initiated from five parsimony starting trees with a topological constraint enforcing monophyly of major angiosperm orders. Branch lengths were estimated under maximum likelihood with proportional scaling across partitions corresponding to genes. 

Ancestral amino acid sequences were [reconstructed](../data/aa_tree_ASR.raxml.log) at all internal nodes of the species tree using marginal ancestral state reconstruction in RAxML-NG v1.2.0. Reconstruction was performed on the tip-collapsed superalignment (6,094 taxa, 11,427 sites) under the LG+G4 model (alpha = 0.637) with branch lengths re-optimized by maximum likelihood. Per-site marginal probabilities were computed for each ancestral node, and the most probable state at each site was taken as the reconstructed ancestral sequence.

### Phylogenetic analysis 

### Protein structural analysis
See: [netsurf prep](../code/prep_seqs_for_netsurf.py) and [structural analysis](../code/structure_3_cuz_2_didntautosave.R)

Per-residue structural features were predicted using the NetSurfP-3.0 webserver. To reduce redundancy while preserving sequence diversity, input sequences for each gene were clustered with CD-HIT (v4.8.1) at a gene-specific identity threshold set to the 90th percentile of sampled pairwise identities from the corresponding alignment (clamped to 0.5–0.99). To accommodate limited sequence submissions, representative sequence per cluster was selected and sequences were batched for prediction. NetSurfP-3.0 output provided per-residue estimates of relative solvent accessibility (RSA), absolute solvent accessibility (ASA), backbone dihedral angles (phi, psi), intrinsic disorder probability, and secondary structure class (helix, sheet, coil). To test whether temperature-associated residues occupy structurally distinct positions, we compared structural feature distributions between GWAS-significant and non-significant sites using Wilcoxon rank-sum tests. 

## Results

### GWAS identifies [N] genome-wide significant loci

<!-- Manhattan plot, top hits table, effect sizes -->
The [GWAS](../code/residueGWAS_PopStepAndaJump.r) identified two sites in the large subunit of RuBisCO—RbcL 328 and RbcL 309—as strongly associated with temperature (P = 10^-32, 10^-30; Figure 1B). Both sites have previously been reported as under positive selection and are evolutionarily coupled. Additional significant sites were also concentrated in rbcL; beyond these, significance was distributed across protein complexes without strong enrichment for any particular complex (Supplementary Figure 1A).
![Figure 1. GWAS overview: residue–temperature associations across chloroplast genes.](figures/fig1_gwas_overview.png)

*Panel A (locus-level):* [fig1_panel_A.png](figures/fig1_panel_A.png) · *QQ plot:* [fig1_panel_qq.png](figures/fig1_panel_qq.png)

| SNP | Chr | Position | Gene | P-value | Beta |
|-----|-----|----------|------|---------|------|
| rs… | … | … | … | … | … |

See: [full results table](data/processed/gwas_top_hits.tsv)

### Phylogenetic conservation of associated genes

<!-- Are the hit genes conserved? Any signatures of selection? -->

![Phylogenetic distribution and conservation of temperature-associated genes.](figures/figure1_panel_B.png)

*Panel C:* [figure1_panel_C.png](figures/figure1_panel_C.png)

### Structural consequences of identified variants

<!-- Do variants map to functional domains? Predicted impact? -->

![Structural enrichment of temperature-associated residues.](figures/structural_enrichment.png)

*Figure 1 panel D (effect-size comparison):* [figure1_panel_D.png](figures/figure1_panel_D.png) · [barplot](figures/figure1_panel_D_barplot.png) · [comparison](figures/figure1_panel_D_comparison.png)

## Discussion

<!-- 
- Summarize key findings
- Biological interpretation — what do the hits mean?
- How phylogenetics and structural data support or contextualize the GWAS
- Limitations (sample size, ancestry, structural prediction confidence)
- Future directions
-->

## Data and code availability

All analysis code is available in [`code/`](code/). Summary statistics are in [`data/processed/`](data/processed/).

## References

<!-- 
Use a consistent format. You can store PDFs in references/ and link them:
1. Smith et al. (2024). *Title*. Journal. [PDF](references/smith2024.pdf)
-->

## Supplementary

- [Supplementary tables](data/processed/)
- [QC plots](figures/)