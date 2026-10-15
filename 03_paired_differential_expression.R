# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 03 - paired differential-expression validation
# Author: Kazhal Ghiasi
# Purpose: Validate the existing paired DESeq2 model and saved results.
#   Expected cohort: 113 patients, one Normal and one Tumour sample each.
# This script does NOT refit DESeq2 or overwrite the original model
# or saved differential-expression results.

# 1. Setup and load existing objects

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("Package DESeq2 is required.")
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("Package SummarizedExperiment is required.")
}

model_file <- "results/deseq2_paired_113.rds"
results_file <- "results/deseq2_paired_113_results.rds"

if (!file.exists(model_file)) {
  stop("Saved model not found: ", model_file)
}

if (!file.exists(results_file)) {
  stop("Saved results not found: ", results_file)
}

dds <- readRDS(model_file)
saved_results <- readRDS(results_file)

if (!inherits(dds, "DESeqDataSet")) {
  stop("The saved model is not a DESeqDataSet.")
}

if (!inherits(saved_results, "DESeqResults")) {
  stop("The saved results are not a DESeqResults object.")
}

# 2. Validate the model design and sample pairing

model_design <- DESeq2::design(dds)

metadata <- as.data.frame(
  SummarizedExperiment::colData(dds)
)

required_metadata <- c("patient", "condition")

if (!all(required_metadata %in% names(metadata))) {
  stop(
    "Model metadata must contain patient and condition columns."
  )
}

if (anyNA(metadata$patient) || anyNA(metadata$condition)) {
  stop("Missing patient or condition labels were found.")
}

metadata$patient <- as.character(metadata$patient)
metadata$condition <- as.character(metadata$condition)

condition_levels <- sort(unique(metadata$condition))

if (!setequal(condition_levels, c("Normal", "Tumour"))) {
  stop(
    "Expected condition labels Normal and Tumour; found: ",
    paste(condition_levels, collapse = ", ")
  )
}

pair_table <- table(
  patient = metadata$patient,
  condition = factor(
    metadata$condition,
    levels = c("Normal", "Tumour")
  )
)

pairing_valid <-
  nrow(pair_table) == 113L &&
  ncol(pair_table) == 2L &&
  all(rowSums(pair_table) == 2L) &&
  all(pair_table[, "Normal"] == 1L) &&
  all(pair_table[, "Tumour"] == 1L)

design_terms <- attr(
  terms(model_design),
  "term.labels"
)

design_valid <- setequal(
  design_terms,
  c("patient", "condition")
) && length(design_terms) == 2L

cat("\n=== PAIRED MODEL VALIDATION ===\n")
cat("Model dimensions (genes x samples):\n")
print(dim(dds))

cat("\nModel design:\n")
print(model_design)

cat("\nCoefficient names:\n")
print(DESeq2::resultsNames(dds))

cat("\nSample counts:\n")
print(table(metadata$condition))

cat("\nUnique patients:", length(unique(metadata$patient)), "\n")
cat("Exactly 113 complete pairs:", pairing_valid, "\n")
cat("Expected paired design:", design_valid, "\n")

if (!pairing_valid) {
  stop(
    "Pairing validation failed. Inspect the saved sample metadata."
  )
}

if (!design_valid) {
  stop(
    "The saved model design differs from the expected paired design."
  )
}

# 3. Extract results from the saved model, without refitting

coefficient_names <- DESeq2::resultsNames(dds)

expected_coefficient <- "condition_Tumour_vs_Normal"

if (!expected_coefficient %in% coefficient_names) {
  stop(
    "Expected coefficient not found: ",
    expected_coefficient,
    "\nAvailable coefficients: ",
    paste(coefficient_names, collapse = ", ")
  )
}

model_results <- DESeq2::results(
  dds,
  name = expected_coefficient
)

saved_df <- as.data.frame(saved_results)
model_df <- as.data.frame(model_results)

saved_gene_ids <- rownames(saved_df)
model_gene_ids <- rownames(model_df)

if (is.null(saved_gene_ids) || is.null(model_gene_ids)) {
  stop("Gene identifiers are missing from one or both result objects.")
}

same_gene_set <- setequal(saved_gene_ids, model_gene_ids)

common_genes <- intersect(saved_gene_ids, model_gene_ids)

columns_to_compare <- intersect(
  c("baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"),
  intersect(names(saved_df), names(model_df))
)

if (length(columns_to_compare) == 0L) {
  stop("No comparable statistical columns were found.")
}

comparison_by_column <- vapply(
  columns_to_compare,
  function(column) {
    isTRUE(
      all.equal(
        saved_df[common_genes, column],
        model_df[common_genes, column],
        tolerance = 1e-8,
        check.attributes = FALSE
      )
    )
  },
  logical(1)
)

results_match <-
  same_gene_set &&
  length(common_genes) == nrow(saved_df) &&
  length(common_genes) == nrow(model_df) &&
  all(comparison_by_column)

cat("\n=== SAVED RESULTS COMPARISON ===\n")
cat("Saved result dimensions:\n")
print(dim(saved_df))

cat("\nExtracted result dimensions:\n")
print(dim(model_df))

cat("\nSame gene set:", same_gene_set, "\n")
cat("Genes compared:", length(common_genes), "\n")
cat("Columns compared:\n")
print(comparison_by_column)

cat("\nAll saved statistics match:", results_match, "\n")

# 4. Summarise the saved differential-expression results

significant_count <- sum(
  saved_df$padj < 0.05,
  na.rm = TRUE
)

summary_table <- data.frame(
  metric = c(
    "Genes in saved model",
    "Samples in saved model",
    "Unique patients",
    "Normal samples",
    "Tumour samples",
    "Exactly 113 complete pairs",
    "Expected paired model design",
    "Genes in saved results",
    "Genes with non-missing adjusted p-values",
    "Genes with adjusted p-value below 0.05",
    "Same gene set in saved and extracted results",
    "All comparable statistics match"
  ),
  value = c(
    nrow(dds),
    ncol(dds),
    length(unique(metadata$patient)),
    sum(metadata$condition == "Normal"),
    sum(metadata$condition == "Tumour"),
    pairing_valid,
    design_valid,
    nrow(saved_df),
    sum(!is.na(saved_df$padj)),
    significant_count,
    same_gene_set,
    results_match
  ),
  stringsAsFactors = FALSE
)

# 5. Write separate diagnostic outputs

write.csv(
  summary_table,
  "results/script03_validation_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    column = names(comparison_by_column),
    matches = unname(comparison_by_column)
  ),
  "results/script03_result_comparison.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    patient = rownames(pair_table),
    Normal = as.integer(pair_table[, "Normal"]),
    Tumour = as.integer(pair_table[, "Tumour"]),
    stringsAsFactors = FALSE
  ),
  "results/script03_pairing_audit.csv",
  row.names = FALSE
)

report_lines <- c(
  "Script 03: paired DESeq2 validation",
  paste("Model file:", model_file),
  paste("Saved results file:", results_file),
  paste("Model design:", paste(deparse(model_design), collapse = "")),
  paste("Genes:", nrow(dds)),
  paste("Samples:", ncol(dds)),
  paste("Unique patients:", length(unique(metadata$patient))),
  paste("Complete 113-pair validation:", pairing_valid),
  paste("Expected model design:", design_valid),
  paste("Same gene set:", same_gene_set),
  paste("Comparable statistics match:", results_match),
  paste("Genes with adjusted p-value < 0.05:", significant_count),
  "",
  "No DESeq2 model refitting was performed.",
  "Original model and results files were not overwritten."
)

writeLines(
  report_lines,
  "results/script03_validation_notes.txt"
)

cat("\n=== SCRIPT 03 COMPLETE ===\n")
cat("Validation summary saved to results/script03_validation_summary.csv\n")
cat("Original model and results remain unchanged.\n")

if (!results_match) {
  warning(
    "Saved and extracted results did not fully match. ",
    "Review script03_result_comparison.csv before interpretation."
  )
}
