# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 05 — Build subtype-specific expression signatures
# Author: Kazhal Ghiasi
# Purpose:
# Identify genes differentially expressed in each major breast cancer subtype
# relative to the equally weighted average of the other three subtypes.
# Method:
# - Use the saved TCGA-BRCA tumour VST expression matrix.
# - Include Basal, Her2, LumA, and LumB tumour subtypes.
# - Fit limma linear models and apply empirical Bayes moderation.
# - Adjust p-values for multiple testing.
# - Retain up to 100 significant genes in each direction per contrast.
# Interpretation:
# These are subtype-versus-rest signatures, not proof of subtype exclusivity.

# 1. Check dependency 

if (!requireNamespace("limma", quietly = TRUE)) {
  stop(
    "Package 'limma' is required. Install it with BiocManager::install('limma'), ",
    "then rerun this script."
  )
}

library(limma)

# 2. Define inputs and outputs 

expr_file <- "results/ml_tumour_vst_matrix_1095.rds"
meta_file <- "results/ml_tumour_metadata_1095.csv"
ann_file <- "results/ml_gene_annotation.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

for (f in c(expr_file, meta_file, ann_file)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f)
  }
}

# 3. Load expression matrix and metadata 

expr <- readRDS(expr_file)

meta <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("expression_id", "sample_id", "subtype")
missing_columns <- setdiff(required_columns, colnames(meta))

if (length(missing_columns) > 0) {
  stop(
    "Metadata is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!is.matrix(expr) || !is.numeric(expr)) {
  stop("The expression object must be a numeric matrix.")
}

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("The expression matrix must have gene and sample identifiers.")
}

if (anyDuplicated(rownames(expr))) {
  stop("Duplicate gene identifiers were found in the expression matrix.")
}

if (anyDuplicated(colnames(expr))) {
  stop("Duplicate sample identifiers were found in the expression matrix.")
}

if (anyDuplicated(meta$expression_id)) {
  stop("Duplicate expression_id values were found in the metadata.")
}

# Explicitly align metadata to expression columns.
match_idx <- match(colnames(expr), meta$expression_id)

if (anyNA(match_idx)) {
  stop("Some expression samples have no matching metadata record.")
}

meta <- meta[match_idx, , drop = FALSE]

stopifnot(identical(colnames(expr), meta$expression_id))

if (any(!is.finite(expr))) {
  stop("The expression matrix contains missing or non-finite values.")
}

cat("\n=== INPUT DATA ===\n")
cat("Genes:", nrow(expr), "\n")
cat("Expression samples:", ncol(expr), "\n")

# 4. Select the four main tumour subtypes 

main_subtypes <- c("Basal", "Her2", "LumA", "LumB")

subtype_text <- trimws(as.character(meta$subtype))
keep <- !is.na(subtype_text) & subtype_text %in% main_subtypes

expr4 <- expr[, keep, drop = FALSE]
meta4 <- meta[keep, , drop = FALSE]

meta4$subtype <- factor(
  trimws(as.character(meta4$subtype)),
  levels = main_subtypes
)

if (nrow(meta4) == 0) {
  stop("No samples matched the four requested subtype labels.")
}

if (any(table(meta4$subtype) == 0)) {
  stop("At least one requested subtype has no samples.")
}

if (any(table(meta4$subtype) < 2)) {
  stop("At least one subtype has fewer than two samples.")
}

cat("\n=== SUBTYPE SAMPLE COUNTS ===\n")
cat("Samples included:", ncol(expr4), "\n")
print(table(meta4$subtype))

# 5. Fit limma model 

design <- model.matrix(~ 0 + subtype, data = meta4)
colnames(design) <- main_subtypes

if (qr(design)$rank != ncol(design)) {
  stop("The subtype design matrix is not full rank.")
}

fit <- limma::lmFit(expr4, design)

# Each subtype is compared with the equally weighted mean of the other three.
contrast_matrix <- limma::makeContrasts(
  Basal_vs_rest = Basal - (Her2 + LumA + LumB) / 3,
  Her2_vs_rest  = Her2 - (Basal + LumA + LumB) / 3,
  LumA_vs_rest  = LumA - (Basal + Her2 + LumB) / 3,
  LumB_vs_rest  = LumB - (Basal + Her2 + LumA) / 3,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

# 6. Load and validate gene annotations 

ann <- read.csv(
  ann_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_annotation_columns <- c("ensembl_id", "gene_symbol", "gene_type")
missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(ann)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "Gene annotation is missing required columns: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

clean_id <- function(x) {
  sub("\\..*$", "", as.character(x))
}

ann$ensembl_id_clean <- clean_id(ann$ensembl_id)

# Prefer annotation rows with a gene symbol when duplicate cleaned IDs exist.
ann_has_symbol <- !is.na(ann$gene_symbol) &
  trimws(as.character(ann$gene_symbol)) != ""

ann <- ann[order(!ann_has_symbol), , drop = FALSE]
ann <- ann[!duplicated(ann$ensembl_id_clean), , drop = FALSE]

annotation_idx <- match(
  clean_id(rownames(expr4)),
  ann$ensembl_id_clean
)

cat("\n=== ANNOTATION CHECK ===\n")
cat("Expression genes:", nrow(expr4), "\n")
cat("Genes matched to annotation:", sum(!is.na(annotation_idx)), "\n")
cat("Genes without annotation match:", sum(is.na(annotation_idx)), "\n")

# 7. Extract complete results and build compact signatures 

all_results <- list()
signature_list <- list()

top_n <- 100L
fdr_threshold <- 0.05

for (contrast_name in colnames(contrast_matrix)) {
  
  res <- limma::topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    sort.by = "P"
  )
  
  res$ensembl_id <- clean_id(rownames(res))
  res$gene_symbol <- ann$gene_symbol[match(
    res$ensembl_id,
    ann$ensembl_id_clean
  )]
  res$gene_type <- ann$gene_type[match(
    res$ensembl_id,
    ann$ensembl_id_clean
  )]
  
  res$subtype_comparison <- contrast_name
  
  significant <- !is.na(res$adj.P.Val) &
    res$adj.P.Val < fdr_threshold
  
  up <- res[significant & res$logFC > 0, , drop = FALSE]
  down <- res[significant & res$logFC < 0, , drop = FALSE]
  
  # Rank significant genes by adjusted p-value, then raw p-value,
  # then absolute effect size. Retain up to top_n per direction.
  rank_signature <- function(x, n) {
    if (nrow(x) == 0) {
      return(x)
    }
    
    x <- x[
      order(x$adj.P.Val, x$P.Value, -abs(x$logFC)),
      ,
      drop = FALSE
    ]
    
    head(x, n)
  }
  
  up_top <- rank_signature(up, top_n)
  down_top <- rank_signature(down, top_n)
  
  selected <- rbind(up_top, down_top)
  
  if (nrow(selected) > 0) {
    selected$signature_direction <- ifelse(
      selected$logFC > 0,
      "Higher_in_subtype",
      "Lower_in_subtype"
    )
  } else {
    selected$signature_direction <- character(0)
  }
  
  all_results[[contrast_name]] <- res
  signature_list[[contrast_name]] <- selected
  
  cat("\n", contrast_name, "\n", sep = "")
  cat("Significant genes (FDR <", fdr_threshold, "): ",
      sum(significant), "\n", sep = "")
  cat("Significantly higher in subtype: ", nrow(up), "\n", sep = "")
  cat("Significantly lower in subtype: ", nrow(down), "\n", sep = "")
  cat("Genes retained in compact signature: ", nrow(selected), "\n", sep = "")
}

all_results_df <- do.call(rbind, all_results)
signatures_df <- do.call(rbind, signature_list)

# 8. Save results 

write.csv(
  all_results_df,
  "results/subtype_differential_expression_all.csv",
  row.names = FALSE
)

write.csv(
  signatures_df,
  "results/subtype_signatures_top100_each_direction.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    expression_id = meta4$expression_id,
    sample_id = meta4$sample_id,
    subtype = as.character(meta4$subtype),
    stringsAsFactors = FALSE
  ),
  "results/subtype_signature_sample_metadata.csv",
  row.names = FALSE
)

saveRDS(
  list(
    design = design,
    contrasts = contrast_matrix,
    fit = fit,
    fit_contrasts = fit2,
    main_subtypes = main_subtypes,
    fdr_threshold = fdr_threshold,
    top_n_per_direction = top_n,
    annotation_matched = sum(!is.na(annotation_idx)),
    annotation_unmatched = sum(is.na(annotation_idx))
  ),
  "results/subtype_signature_analysis_details.rds"
)

cat("\n=== SCRIPT 05 COMPLETE ===\n")
cat("Saved complete differential-expression results and compact signatures.\n")

