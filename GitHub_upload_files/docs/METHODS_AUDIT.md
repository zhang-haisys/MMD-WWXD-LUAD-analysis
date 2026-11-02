# Methods audit / 方法核查

These are findings in the supplied source, not newly executed analyses.

## Blocking issues

1. Preprocessing module 00: `tcga <- log2(tcga+1)` is followed by merging and
   `Expr=log2(Expr+1)`. Sequential execution can double-transform TCGA values.
   Historical execution and pre-ComBat inputs must be checked. GEOquery download
   alone does not establish the submitted microarray normalization or log scale.
2. Filtering uses `rowSums(exp)>1` on the merged matrix; platform-specific filtering
   rationale and exact input scales remain unverified.
3. `ComBat(Expr, batch=group_list$dataset)` contains no explicit `mod` preserving
   tumor/control or clinical covariates. No covariates were silently added here.
4. Source modules use `is_count=FALSE`, but the complete original non-count DEG
   helper was not supplied. The separately reviewed partial edgeR helper is not
   a valid replacement and is not bundled as the missing original helper.
5. Module 00 contains `inner_join(..., by="X")` while upstream gene identifiers
   are row names; a missing intermediate conversion must be resolved.
6. TCGA survival time is in days and GSE11969 is read in months. A shared time
   unit is not documented in supplied preprocessing. Verify before pooling.
7. Module 01 sets `case="Parental"`, `control="OR"` in all three comparisons.
   This is the reverse of an OR-minus-parental contrast if the missing helper
   implements case-minus-control. Confirm the helper before finalizing direction.
   The Venn uses all significant genes, regardless of direction. No requirement
   of direction agreement was imposed by that intersection code.
8. Module 14 contains GENE1/GENE2/GENE3, example coefficients, uncreated input
   objects, survival analysis and zero filling of missing genes. It is not a
   verified external tumor/normal classifier and must not be reported as one.
9. Module 04 includes GSE54236 references although the manuscript names GSE31210.
   Do not rename an accession to make it match without using the correct data.
10. Module 09 uses `Masked Somatic Mutation_1` and inconsistent Mutation paths;
    the GDC query/workflow and actual variant manifest need verification.
11. Module 07 computes hub-gene correlations, not a single-gene GSEA workflow.
12. Module 11 exports expression and reads a TIDE.csv; TIDE scores are not
    calculated by that R file. Record the external service/version and procedure.

## Interpretation

Packaging is not a scientific correction. Existing DEGs, target genes, scores,
cellular results and biological interpretation were not recomputed or endorsed.
Where preprocessing is corrected, affected downstream analyses must be compared.
Current sessionInfo cannot establish the environment used historically.
Do not describe this release as fully reproducible until inputs, helpers,
environment, paths and outputs have been verified together.
