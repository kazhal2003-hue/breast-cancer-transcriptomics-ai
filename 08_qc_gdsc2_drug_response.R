# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 08 — Quality control of GDSC2 breast cancer drug-response data
# Author: Kazhal Ghiasi
# Purpose: Inspect GDSC2 drug-response completeness, repeated records, cell-line
# coverage, drug coverage, and overlap with the expression dataset.
# Important: This is a QC script. It does not collapse repeated records or select the
# final drug-response measurement for downstream association testing.
# Coverage thresholds for downstream analyses:
# Primary: at least 20 matched cell lines per drug
# Sensitivity: at least 15 and at least 25 matched cell lines per drug

# 1. Check dependencies 

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Package 'readxl' is required. Install it with ",
    "install.packages('readxl'), then rerun this script."
  )
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop(
    "Package 'dplyr' is required. Install it with ",
    "install.packages('dplyr'), then rerun this script."
  )
}

# 2. Define input files 

project_dir <- getwd()

response_file <- file.path(
  "data", "drug_response",
  "GDSC2_fitted_dose_response_27Oct23.xlsx"
)

expression_file <- file.path(
  "results", "gdsc2_breast_cell_line_tpm_51.rds"
)

results_dir <- file.path(project_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(response_file)) {
  stop("Drug-response workbook not found: ", response_file)
}

if (!file.exists(expression_file)) {
  stop(
    "Matched expression dataset not found: ", expression_file,
    ". Run Script 07 successfully before running Script 08."
  )
}

# 3. Load response records and validate columns 

gdsc <- readxl::read_excel(response_file)

required_columns <- c(
  "DATASET", "CELL_LINE_NAME", "SANGER_MODEL_ID",
  "CANCER_TYPE", "DRUG_ID", "DRUG_NAME",
  "LN_IC50", "AUC", "RMSE", "Z_SCORE",
  "NLME_RESULT_ID", "NLME_CURVE_ID"
)

missing_columns <- setdiff(required_columns, names(gdsc))

if (length(missing_columns) > 0) {
  stop(
    "Expected columns not found in GDSC2 workbook: ",
    paste(missing_columns, collapse = ", ")
  )
}

# 4. Identify breast cancer response records 

cancer_type <- as.character(gdsc$CANCER_TYPE)

is_breast <- !is.na(cancer_type) &
  grepl("breast", cancer_type, ignore.case = TRUE)

breast_gdsc <- gdsc[is_breast, , drop = FALSE]

if (nrow(breast_gdsc) == 0) {
  stop("No breast cancer records were identified in the response workbook.")
}

# Standardise IDs for matching without modifying the original source columns.
breast_gdsc$SANGER_MODEL_ID <- trimws(
  as.character(breast_gdsc$SANGER_MODEL_ID)
)

breast_gdsc$DRUG_ID <- as.character(breast_gdsc$DRUG_ID)

cat("\n=== DATA OVERVIEW ===\n")
cat("All response records:", nrow(gdsc), "\n")
cat("Breast cancer response records:", nrow(breast_gdsc), "\n")
cat(
  "Unique breast cancer model IDs:",
  dplyr::n_distinct(
    breast_gdsc$SANGER_MODEL_ID[
      !is.na(breast_gdsc$SANGER_MODEL_ID) &
        nzchar(breast_gdsc$SANGER_MODEL_ID)
    ]
  ),
  "\n"
)
cat(
  "Unique drug IDs:",
  dplyr::n_distinct(
    breast_gdsc$DRUG_ID[!is.na(breast_gdsc$DRUG_ID)]
  ),
  "\n"
)

cat("\nRecords by dataset:\n")
print(table(breast_gdsc$DATASET, useNA = "ifany"))

# 5. Match response models to the expression dataset 

expression_data <- readRDS(expression_file)

if (!is.data.frame(expression_data) || ncol(expression_data) < 4) {
  stop("The saved expression dataset has an unexpected structure.")
}

expression_model_ids <- trimws(
  names(expression_data)[-(1:3)]
)

if (anyDuplicated(expression_model_ids)) {
  stop("Duplicate cell-line IDs found in the expression dataset.")
}

response_model_ids <- unique(
  breast_gdsc$SANGER_MODEL_ID[
    !is.na(breast_gdsc$SANGER_MODEL_ID) &
      nzchar(breast_gdsc$SANGER_MODEL_ID)
  ]
)

matched_model_ids <- intersect(
  response_model_ids,
  expression_model_ids
)

unmatched_response_ids <- setdiff(
  response_model_ids,
  expression_model_ids
)

cat("\n=== EXPRESSION/RESPONSE MODEL OVERLAP ===\n")
cat("Expression cell lines:", length(expression_model_ids), "\n")
cat("Breast response cell lines:", length(response_model_ids), "\n")
cat("Matched cell lines:", length(matched_model_ids), "\n")
cat("Response models without expression:", length(unmatched_response_ids), "\n")

# Preserve all breast response rows for the full QC output.
# Also save the subset that can be matched to expression data.
breast_gdsc_matched <- breast_gdsc[
  !is.na(breast_gdsc$SANGER_MODEL_ID) &
    breast_gdsc$SANGER_MODEL_ID %in% matched_model_ids,
  ,
  drop = FALSE
]

# 6. Missingness summary 

qc_columns <- c(
  "DATASET", "CELL_LINE_NAME", "SANGER_MODEL_ID",
  "DRUG_ID", "DRUG_NAME", "LN_IC50", "AUC",
  "RMSE", "Z_SCORE", "NLME_RESULT_ID", "NLME_CURVE_ID"
)

missing_summary <- data.frame(
  column = qc_columns,
  missing_n = vapply(
    breast_gdsc[qc_columns],
    function(x) sum(is.na(x)),
    numeric(1)
  ),
  missing_percent = vapply(
    breast_gdsc[qc_columns],
    function(x) mean(is.na(x)) * 100,
    numeric(1)
  )
)

cat("\n=== MISSING VALUES ===\n")
print(missing_summary, row.names = FALSE)

# 7. Inspect repeated cell-line/drug combinations 

pair_counts <- dplyr::count(
  breast_gdsc,
  DATASET, SANGER_MODEL_ID, DRUG_ID,
  name = "n_records"
)

pair_counts_across_datasets <- dplyr::count(
  breast_gdsc,
  SANGER_MODEL_ID, DRUG_ID,
  name = "n_records"
)

cat("\n=== REPEATED COMBINATIONS ===\n")
cat("Dataset/model/drug combinations:", nrow(pair_counts), "\n")
cat(
  "Combinations with multiple records within dataset:",
  sum(pair_counts$n_records > 1), "\n"
)
cat(
  "Model/drug combinations with multiple records across datasets:",
  sum(pair_counts_across_datasets$n_records > 1), "\n"
)

cat("\nDistribution of records per dataset/model/drug combination:\n")
print(table(pair_counts$n_records))

# Check exact repeated rows separately. Multiple records are not necessarily
# erroneous: they may represent different curve fits or response measurements.
response_columns <- intersect(
  c(
    "DATASET", "CELL_LINE_NAME", "SANGER_MODEL_ID",
    "DRUG_ID", "DRUG_NAME", "LN_IC50", "AUC", "RMSE",
    "Z_SCORE", "NLME_RESULT_ID", "NLME_CURVE_ID"
  ),
  names(breast_gdsc)
)

exact_duplicate_rows <- duplicated(
  as.data.frame(breast_gdsc[, response_columns, drop = FALSE])
)

cat(
  "\nExact duplicate rows across selected response fields:",
  sum(exact_duplicate_rows), "\n"
)

# 8. Cell-line coverage 

model_summary <- dplyr::summarise(
  dplyr::group_by(
    breast_gdsc,
    SANGER_MODEL_ID, CELL_LINE_NAME
  ),
  n_records = dplyr::n(),
  n_drugs = dplyr::n_distinct(DRUG_ID, na.rm = TRUE),
  n_missing_LN_IC50 = sum(is.na(LN_IC50)),
  n_missing_AUC = sum(is.na(AUC)),
  .groups = "drop"
)

matched_model_summary <- dplyr::summarise(
  dplyr::group_by(
    breast_gdsc_matched,
    SANGER_MODEL_ID, CELL_LINE_NAME
  ),
  n_records = dplyr::n(),
  n_drugs = dplyr::n_distinct(DRUG_ID, na.rm = TRUE),
  n_missing_LN_IC50 = sum(is.na(LN_IC50)),
  n_missing_AUC = sum(is.na(AUC)),
  .groups = "drop"
)

cat("\n=== CELL-LINE COVERAGE ===\n")
cat("Cell-line summary rows:", nrow(model_summary), "\n")
cat("Matched cell-line summary rows:", nrow(matched_model_summary), "\n")

# 9. Drug coverage 

drug_summary <- dplyr::summarise(
  dplyr::group_by(
    breast_gdsc,
    DRUG_ID, DRUG_NAME
  ),
  n_records = dplyr::n(),
  n_cell_lines = dplyr::n_distinct(SANGER_MODEL_ID, na.rm = TRUE),
  n_missing_LN_IC50 = sum(is.na(LN_IC50)),
  n_missing_AUC = sum(is.na(AUC)),
  .groups = "drop"
)

matched_drug_summary <- dplyr::summarise(
  dplyr::group_by(
    breast_gdsc_matched,
    DRUG_ID, DRUG_NAME
  ),
  n_records = dplyr::n(),
  n_cell_lines = dplyr::n_distinct(SANGER_MODEL_ID, na.rm = TRUE),
  n_missing_LN_IC50 = sum(is.na(LN_IC50)),
  n_missing_AUC = sum(is.na(AUC)),
  .groups = "drop"
)

cat("\n=== DRUG COVERAGE: ALL BREAST RESPONSE RECORDS ===\n")
cat("Drug summary rows:", nrow(drug_summary), "\n")
cat("\nDrugs tested in the fewest cell lines:\n")
print(
  head(drug_summary[order(drug_summary$n_cell_lines), ], 20),
  row.names = FALSE
)

cat("\n=== DRUG COVERAGE: MATCHED EXPRESSION MODELS ===\n")
cat("Drug summary rows:", nrow(matched_drug_summary), "\n")

# 10. Assess prespecified drug coverage thresholds 

coverage_thresholds <- c(15L, 20L, 25L)

threshold_summary <- do.call(
  rbind,
  lapply(coverage_thresholds, function(threshold) {
    data.frame(
      minimum_matched_cell_lines = threshold,
      drugs_passing_threshold = sum(
        matched_drug_summary$n_cell_lines >= threshold
      )
    )
  })
)

cat("\n=== MATCHED DRUG COVERAGE THRESHOLDS ===\n")
print(threshold_summary, row.names = FALSE)

matched_drug_summary$passes_primary_20 <- (
  matched_drug_summary$n_cell_lines >= 20
)

# 11. Save QC outputs 

write.csv(
  missing_summary,
  file.path(results_dir, "gdsc2_breast_response_missingness.csv"),
  row.names = FALSE
)

write.csv(
  pair_counts,
  file.path(results_dir, "gdsc2_breast_response_pair_counts.csv"),
  row.names = FALSE
)

write.csv(
  pair_counts_across_datasets,
  file.path(results_dir, "gdsc2_breast_response_pair_counts_across_datasets.csv"),
  row.names = FALSE
)

write.csv(
  model_summary,
  file.path(results_dir, "gdsc2_breast_response_model_coverage.csv"),
  row.names = FALSE
)

write.csv(
  matched_model_summary,
  file.path(results_dir, "gdsc2_breast_response_matched_model_coverage.csv"),
  row.names = FALSE
)

write.csv(
  drug_summary,
  file.path(results_dir, "gdsc2_breast_response_drug_coverage.csv"),
  row.names = FALSE
)

write.csv(
  matched_drug_summary,
  file.path(results_dir, "gdsc2_breast_response_matched_drug_coverage.csv"),
  row.names = FALSE
)

write.csv(
  threshold_summary,
  file.path(results_dir, "gdsc2_breast_response_coverage_thresholds.csv"),
  row.names = FALSE
)

# Preserve the original full breast-response output for compatibility.
saveRDS(
  breast_gdsc,
  file.path(results_dir, "gdsc2_breast_response_records.rds")
)

# Save the matched subset separately for downstream analyses.
saveRDS(
  breast_gdsc_matched,
  file.path(results_dir, "gdsc2_breast_response_records_matched_expression.rds")
)

saveRDS(
  list(
    response_file = response_file,
    expression_file = expression_file,
    n_breast_response_records = nrow(breast_gdsc),
    n_matched_response_records = nrow(breast_gdsc_matched),
    n_response_models = length(response_model_ids),
    n_expression_models = length(expression_model_ids),
    matched_model_ids = matched_model_ids,
    unmatched_response_model_ids = unmatched_response_ids,
    coverage_thresholds = threshold_summary,
    exact_duplicate_rows_in_selected_fields = sum(exact_duplicate_rows)
  ),
  file.path(results_dir, "gdsc2_breast_response_qc_details.rds")
)

cat("\n=== SCRIPT 08 COMPLETE ===\n")
cat("Full QC and matched-expression outputs saved in:", results_dir, "\n")

