# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 06 — Refine subtype signatures by effect size
# Author: Kazhal Ghiasi
# Purpose:
# Test the sensitivity of subtype-specific expression signatures to
# increasingly stringent effect-size thresholds.
# Primary filter:
# - Benjamini-Hochberg adjusted p-value < 0.05
# - Absolute log2 fold change >= 1.0
# Sensitivity thresholds:
# - |log2FC| >= 0.5
# - |log2FC| >= 1.0
# - |log2FC| >= 1.5
# The primary filtered table retains all genes passing the filters.
# A separate compact signature retains up to 100 higher and 100 lower
# expressed genes per subtype comparison.

# 1. Define inputs and outputs 

input_file <- "results/subtype_differential_expression_all.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop(
    "Required input not found: ", input_file,
    ". Run Script 05 successfully before running Script 06."
  )
}

# 2. Load and validate differential-expression results 

results <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "logFC",
  "adj.P.Val",
  "subtype_comparison",
  "gene_symbol",
  "ensembl_id"
)

missing_columns <- setdiff(required_columns, names(results))

if (length(missing_columns) > 0) {
  stop(
    "Input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(results) == 0) {
  stop("The differential-expression results file contains no rows.")
}

results$logFC <- as.numeric(results$logFC)
results$adj.P.Val <- as.numeric(results$adj.P.Val)
results$subtype_comparison <- as.character(results$subtype_comparison)

comparisons <- sort(unique(
  results$subtype_comparison[
    !is.na(results$subtype_comparison) &
      nzchar(results$subtype_comparison)
  ]
))

if (length(comparisons) == 0) {
  stop("No valid subtype comparisons were found in the input.")
}

cat("\n=== INPUT CHECK ===\n")
cat("Genes across all comparisons:", nrow(results), "\n")
cat("Subtype comparisons:", length(comparisons), "\n")
print(comparisons)

# 3. Define thresholds and significance criteria 

thresholds <- c(0.5, 1.0, 1.5)
fdr_cutoff <- 0.05
effect_cutoff <- 1.0
top_n_per_direction <- 100L

# |log2FC| >= 0.5: approximately 1.41-fold difference
# |log2FC| >= 1.0: 2-fold difference
# |log2FC| >= 1.5: approximately 2.83-fold difference

# 4. Summarise sensitivity to effect-size thresholds 

threshold_rows <- list()

for (comparison in comparisons) {
  
  x <- results[
    !is.na(results$subtype_comparison) &
      results$subtype_comparison == comparison,
    ,
    drop = FALSE
  ]
  
  for (threshold in thresholds) {
    
    keep <- !is.na(x$adj.P.Val) &
      !is.na(x$logFC) &
      x$adj.P.Val < fdr_cutoff &
      abs(x$logFC) >= threshold
    
    threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
      subtype_comparison = comparison,
      log2FC_threshold = threshold,
      significant_genes = sum(keep),
      higher_in_subtype = sum(keep & x$logFC > 0, na.rm = TRUE),
      lower_in_subtype = sum(keep & x$logFC < 0, na.rm = TRUE)
    )
  }
}

threshold_summary <- do.call(rbind, threshold_rows)
rownames(threshold_summary) <- NULL

write.csv(
  threshold_summary,
  "results/subtype_signature_threshold_sensitivity.csv",
  row.names = FALSE
)

cat("\n=== GENE COUNTS AT DIFFERENT EFFECT-SIZE THRESHOLDS ===\n")
print(threshold_summary, row.names = FALSE)

# 5. Build the primary effect-filtered signatures 

significant <- !is.na(results$adj.P.Val) &
  !is.na(results$logFC) &
  !is.na(results$subtype_comparison) &
  results$adj.P.Val < fdr_cutoff

passes_effect <- abs(results$logFC) >= effect_cutoff

filtered <- results[
  significant & passes_effect,
  ,
  drop = FALSE
]

# Add explicit direction for downstream interpretation.
filtered$signature_direction <- ifelse(
  filtered$logFC > 0,
  "Higher_in_subtype",
  ifelse(
    filtered$logFC < 0,
    "Lower_in_subtype",
    "Zero_effect"
  )
)

# Sort deterministically within each comparison.
if (nrow(filtered) > 0) {
  filtered <- filtered[
    order(
      filtered$subtype_comparison,
      filtered$adj.P.Val,
      -abs(filtered$logFC),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}

write.csv(
  filtered,
  "results/subtype_signatures_effect_filtered_all.csv",
  row.names = FALSE
)

# 6. Create compact signatures: up to 100 genes per direction 

signature_list <- list()
compact_summary <- list()

for (comparison in comparisons) {
  
  x <- filtered[
    filtered$subtype_comparison == comparison,
    ,
    drop = FALSE
  ]
  
  up <- x[x$logFC > 0, , drop = FALSE]
  down <- x[x$logFC < 0, , drop = FALSE]
  
  # Higher-expression genes: strongest positive effects first after FDR.
  if (nrow(up) > 0) {
    up <- up[
      order(up$adj.P.Val, -up$logFC, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }
  
  # Lower-expression genes: strongest negative effects first after FDR.
  if (nrow(down) > 0) {
    down <- down[
      order(down$adj.P.Val, down$logFC, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }
  
  selected <- rbind(
    head(up, top_n_per_direction),
    head(down, top_n_per_direction)
  )
  
  signature_list[[comparison]] <- selected
  
  compact_summary[[comparison]] <- data.frame(
    subtype_comparison = comparison,
    genes_passing_primary_filter = nrow(x),
    genes_higher_before_top_n = nrow(up),
    genes_lower_before_top_n = nrow(down),
    genes_retained_in_compact_signature = nrow(selected)
  )
  
  cat(
    "\n", comparison,
    "\nGenes passing primary filters: ", nrow(x),
    "\nHigher in subtype before top-N limit: ", nrow(up),
    "\nLower in subtype before top-N limit: ", nrow(down),
    "\nGenes retained in compact signature: ", nrow(selected),
    "\n",
    sep = ""
  )
}

compact_signatures <- do.call(rbind, signature_list)
rownames(compact_signatures) <- NULL

compact_summary_df <- do.call(rbind, compact_summary)
rownames(compact_summary_df) <- NULL

write.csv(
  compact_signatures,
  "results/subtype_signatures_effect_filtered_top100.csv",
  row.names = FALSE
)

write.csv(
  compact_summary_df,
  "results/subtype_signature_effect_filter_summary.csv",
  row.names = FALSE
)

# 7. Save analysis settings 

saveRDS(
  list(
    input_file = input_file,
    subtype_comparisons = comparisons,
    fdr_cutoff = fdr_cutoff,
    effect_size_thresholds = thresholds,
    primary_effect_cutoff = effect_cutoff,
    top_n_per_direction = top_n_per_direction,
    threshold_summary = threshold_summary,
    compact_signature_summary = compact_summary_df
  ),
  "results/subtype_signature_refinement_details.rds"
)

cat("\n=== SCRIPT 06 COMPLETE ===\n")
cat("Saved threshold sensitivity, filtered genes, compact signatures, and settings.\n")

