---
name: gwas-embedding-orthogonality
description: Explains GWAS sites vs protein embeddings, orthogonal sites, and population structure confounding in genomic prediction. Use when combining GWAS and embedding features.
---

# GWAS vs Embedding Orthogonality

## Context

- **GWAS sites**: Binary features from amino acid states at significant positions. Interpretable, sparse.
- **Embeddings**: Dense protein sequence representations (e.g. ESM). Capture broader context.

## Orthogonal Sites

Sites are "orthogonal" to embeddings if they are uncorrelated with embedding dimensions. Use orthogonal sites to add signal without duplicating embedding information. Sites correlated with embeddings may just reflect structure captured by embeddings.

## Workflow

1. Run GWAS to get significant sites.
2. Compute correlation between GWAS site features and embedding dimensions.
3. Classify: `sig_control` = significant with control; low correlation with embeddings.
4. Use orthogonal sites as additional predictors in combined model.

## Confounding

Population structure confounds GWAS: related taxa share genotype and phenotype. Orthogonal sites + phylogenetic holdout reduce confounding. Adding all GWAS sites may increase structure capture.
