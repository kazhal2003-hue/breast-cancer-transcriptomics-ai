# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 03 - paired differential expression
# Author: Kazhal Ghiasi
# Purpose: Restore and verify the saved paired DESeq2 analysis.
# Dataset: TCGA-BRCA; 113 matched tumour-normal patient pairs.
# This script restores saved objects and performs diagnostic checks.
# It does not refit the DESeq2 model.

# 1. Setup 

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required but is not installed.")
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

model_file <- "results/deseq2_paired_113.rds"
results_file <- "results/deseq2_paired_113_results.rds"

if (!file.exists(model_file)) {
  stop("Saved paired model not found: ", model_file)
}
if (!file.exists(results_file)) {
  stop("Saved paired results not found: ", results_file)
}

dds_paired <- readRDS(model_file)
res_paired <- readRDS(results_file)

if (!inherits(dds_paired, "DESeqDataSet")) {
  stop("The saved model is not a DESeqDataSet.")
}
if (!inherits(res_paired, "DESeqResults")) {
  stop("The saved results are not a DESeqResults object.")
}

# 2. Verify design and patient pairing 

model_design <- DESeq2::design(dds_paired)
sample_metadata <- as.data.frame(
  SummarizedExperiment::colData(dds_paired)
)

if (!all(c("patient", "condition") %in% names(sample_metadata))) {
  stop("Expected patient and condition columns are missing.")
}

pair_table <- with(
  sample_metadata,
  table(patient, condition, useNA = "ifany")
)

pairing_valid <- (
  ncol(pair_table) == 2L &&
    all(rowSums(pair_table) == 2L) &&
    all(rowSums(pair_table == 1L) == 2L) &&
    !anyNA(sample_metadata$patient) &&
    !anyNA(sample_metadata$condition)
)

expected_design <- setequal(
  attr(terms(model_design), "term.labels"),
  c("patient", "condition")
) && length(attr(terms(model_design), "term.labels")) == 2L

cat("\n=== MODEL DESIGN ===\n")
print(model_design)

cat("\n=== MODEL DIMENSIONS (GENES x SAMPLES) ===\n")
print(dim(dds_paired))

cat("\n=== SAMPLE COUNTS BY CONDITION ===\n")
print(table(sample_metadata$condition, useNA = "ifany"))

cat("\n=== PATIENT-LEVEL PAIRING ===\n")
cat("Unique patients:", length(unique(sample_metadata$patient)), "\n")
cat("Valid one-normal/one-tumour pairing:", pairing_valid, "\n")
cat("Expected model design:", expected_design, "\n")

write.csv(
  sample_metadata,
  "results/paired_sample_metadata.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(pair_table),
  "results/paired_patient_condition_table.csv",
  row.names = FALSE
)

# 3. Check convergence information, if available 

model_metadata <- S4Vectors::mcols(dds_paired)

if ("betaConv" %in% names(model_metadata)) {
  beta_conv <- model_metadata$betaConv
  convergence_failures <- sum(!beta_conv, na.rm = TRUE)
  convergence_missing <- sum(is.na(beta_conv))
} else {
  convergence_failures <- NA_integer_
  convergence_missing <- NA_integer_
}

cat("\n=== CONVERGENCE DIAGNOSTIC ===\n")
cat("Genes with coefficient convergence failures:",
    convergence_failures, "\n")
cat("Genes with unavailable convergence status:",
    convergence_missing, "\n")

# 4. Compare the saved results with results extracted from the model 

model_results <- tryCatch(
  DESeq2::results(dds_paired, name = "condition_Tumour_vs_Normal"),
  error = function(e) NULL
)

# A coefficient name may differ by factor reference level or naming convention.
if (is.null(model_results)) {
  available_coefficients <- tryCatch(
    DESeq2::resultsNames(dds_paired),
    error = function(e) character(0)
  )
  
  cat("\nAvailable model coefficients:\n")
  print(available_coefficients)
  
  comparison_status <- paste(
    "Could not extract results using the expected coefficient name;",
    "manual coefficient verification is required."
  )
  results_match <- NA
} else {
  saved_df <- as.data.frame(res_paired)
  model_df <- as.data.frame(model_results)
  
  common_genes <- intersect(rownames(saved_df), rownames(model_df))
  
  if (length(common_genes) == 0L) {
    results_match <- FALSE
  } else {
    columns_to_compare <- intersect(
      c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
      intersect(names(saved_df), names(model_df))
    )
    
    results_match <- length(columns_to_compare) > 0L &&
      all(vapply(columns_to_compare, function(column) {
        isTRUE(all.equal(
          saved_df[common_genes, column],
          model_df[common_genes, column],
          check.attributes = FALSE
        ))
      }, logical(1)))
  }
  
  comparison_status <- if (isTRUE(results_match)) {
    "Saved results match results extracted from the saved model."
  } else {
    "Saved results do not match, or could not be fully matched to, the extracted results. Review required."
  }
}

cat("\n=== SAVED RESULTS CHECK ===\n")
cat("Saved results dimensions:\n")
print(dim(as.data.frame(res_paired)))
cat(comparison_status, "\n")

# 5. Save results, summaries, and diagnostics 

res_paired_df <- as.data.frame(res_paired)
res_paired_df$gene_id <- rownames(res_paired_df)
res_paired_df <- res_paired_df[
  c("gene_id", setdiff(names(res_paired_df), "gene_id"))
]

write.csv(
  res_paired_df,
  "results/deseq2_paired_113_results.csv",
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "Number of genes in saved model",
    "Number of samples in saved model",
    "Number of unique patients",
    "Valid patient pairing",
    "Expected model design",
    "Genes with non-missing adjusted p-values",
    "Genes with adjusted p-value < 0.05",
    "Genes with coefficient convergence failures",
    "Saved results match extracted model results"
  ),
  value = c(
    nrow(dds_paired),
    ncol(dds_paired),
    length(unique(sample_metadata$patient)),
    pairing_valid,
    expected_design,
    sum(!is.na(res_paired$padj)),
    sum(res_paired$padj < 0.05, na.rm = TRUE),
    convergence_failures,
    results_match
  )
)

write.csv(
  summary_table,
  "results/paired_deseq2_summary.csv",
  row.names = FALSE
)

writeLines(
  c(
    "Paired DESeq2 diagnostic checkpoint",
    paste("Model design:", paste(deparse(model_design), collapse = "")),
    paste("Valid patient pairing:", pairing_valid),
    paste("Expected model design:", expected_design),
    paste("Results comparison:", comparison_status),
    paste("Convergence failures:", convergence_failures),
    "The model was not refitted by this script."
  ),
  "results/paired_deseq2_checkpoint_notes.txt"
)

cat("\n=== SCRIPT 03 COMPLETE ===\n")
cat("Saved model and results restored; no model refitting performed.\n")
cat("Review results/paired_deseq2_summary.csv and checkpoint notes.\n")
