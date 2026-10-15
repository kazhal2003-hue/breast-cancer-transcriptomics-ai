# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 07 - Prepare GDSC2 expression
# Author: Kazhal Ghiasi
# Purpose: Extract and validate GDSC2 expression profiles for breast-cancer
# cell lines while preserving established outputs.
# Missing expression values are reported, not silently imputed.
# Existing standard outputs are never overwritten by a differing candidate.

options(stringsAsFactors = FALSE)

rna_zip <- "data/drug_response/rnaseq_merged_20260323.zip"
tpm_member <- "rnaseq_merged_rsem_tpm_20260323.csv"
sources_file <- "results/gdsc2_breast_model_sources.csv"

out_rds <- "results/gdsc2_breast_cell_line_tpm_51.rds"
out_csv <- "results/gdsc2_breast_cell_line_tpm_51.csv"
candidate_rds <- "results/script07_candidate_gdsc2_breast_cell_line_tpm.rds"
candidate_csv <- "results/script07_candidate_gdsc2_breast_cell_line_tpm.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rna_zip)) stop("Missing source ZIP: ", rna_zip)
if (!file.exists(sources_file)) stop("Missing model-source metadata: ", sources_file)

# 1. Read the three metadata rows containing the source model IDs,
#    model names, and data sources.

con <- unz(rna_zip, tpm_member, open = "r")
metadata_lines <- tryCatch(
  readLines(con, n = 3, warn = FALSE),
  finally = close(con)
)

if (length(metadata_lines) != 3L) {
  stop("Could not read all three metadata rows from the source CSV.")
}

split_fields <- function(z) strsplit(z, ",", fixed = TRUE)[[1]]
metadata <- lapply(metadata_lines, split_fields)
field_counts <- vapply(metadata, length, integer(1))

if (length(unique(field_counts)) != 1L) {
  stop("Metadata rows have inconsistent numbers of fields.")
}

model_ids <- trimws(metadata[[1]][-(1:3)])
model_names <- trimws(metadata[[2]][-(1:3)])
data_sources <- trimws(metadata[[3]][-(1:3)])

if (!length(model_ids) || anyNA(model_ids) || any(model_ids == "")) {
  stop("Missing model IDs in metadata row 1.")
}
if (anyDuplicated(model_ids)) {
  stop("Duplicate model IDs detected in metadata row 1.")
}
if (anyNA(model_names) || any(model_names == "") ||
    anyNA(data_sources) || any(data_sources == "")) {
  stop("Missing model names or data sources in metadata.")
}

model_metadata <- data.frame(
  model_id = model_ids,
  model_name_full = model_names,
  data_source_full = data_sources,
  stringsAsFactors = FALSE
)

# 2. Load and validate the established breast-cancer model list.

sources <- read.csv(
  sources_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "model_id", "model_name_matched", "data_source_matched",
  "model_name_full", "data_source_full"
)

if (!all(required_cols %in% names(sources))) {
  stop("Model-source metadata is missing required columns.")
}

sources$model_id <- trimws(sources$model_id)

if (anyNA(sources$model_id) || any(sources$model_id == "") ||
    anyDuplicated(sources$model_id)) {
  stop("Missing or duplicated model IDs in breast model-source metadata.")
}

missing_ids <- setdiff(sources$model_id, model_ids)

if (length(missing_ids)) {
  stop(
    "Breast model IDs absent from source expression metadata: ",
    paste(missing_ids, collapse = ", ")
  )
}

# Keep the source-file order while reading only the selected models.
selected_source_ids <- model_ids[model_ids %in% sources$model_id]

# The source has three annotation columns followed by expression columns.
# NULL tells read.csv() to skip unwanted cell-line columns.
column_classes <- rep("NULL", length(model_ids) + 3L)
column_classes[1:3] <- "character"
selected_positions <- match(selected_source_ids, model_ids) + 3L
column_classes[selected_positions] <- "character"

# 3. Read gene annotations and selected expression columns only.

con <- unz(rna_zip, tpm_member, open = "r")
expression <- tryCatch(
  read.csv(
    con,
    skip = 3,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    colClasses = column_classes,
    na.strings = character()
  ),
  finally = close(con)
)

if (ncol(expression) != length(selected_source_ids) + 3L) {
  stop("Unexpected number of selected expression columns.")
}

if (!identical(
  names(expression)[1:3],
  c("gene_symbol", "ensembl_gene_id", "gene_id")
)) {
  stop("Unexpected gene-annotation columns in the expression table.")
}

gene_symbol <- trimws(expression$gene_symbol)
ensembl_gene_id <- trimws(expression$ensembl_gene_id)
gene_id <- trimws(expression$gene_id)

if (anyNA(gene_symbol) || any(gene_symbol == "")) {
  stop("Missing gene symbols detected.")
}

# Convert expression values to numeric and detect failed conversions.
expr_raw <- expression[-(1:3)]
expr_numeric <- lapply(
  expr_raw,
  function(z) suppressWarnings(as.numeric(z))
)

bad_numeric <- vapply(seq_along(expr_raw), function(i) {
  raw <- expr_raw[[i]]
  converted <- expr_numeric[[i]]
  any(!is.na(raw) & raw != "" & is.na(converted))
}, logical(1))

if (any(bad_numeric)) {
  stop("Non-numeric values found in expression columns.")
}

expr_all <- as.data.frame(expr_numeric, check.names = FALSE)
names(expr_all) <- selected_source_ids

# Restore the established model order from the metadata file.
breast_expression <- expr_all[, sources$model_id, drop = FALSE]

# 4. Validate values. Missing values are counted separately from infinities.
expression_matrix <- as.matrix(breast_expression)

missing_values <- sum(is.na(expression_matrix))
infinite_values <- sum(is.infinite(expression_matrix), na.rm = TRUE)
nan_values <- sum(is.nan(expression_matrix), na.rm = TRUE)
negative_values <- sum(expression_matrix < 0, na.rm = TRUE)

message("Selected breast models: ", ncol(breast_expression))
message("Expression genes: ", nrow(breast_expression))
message("Missing expression values: ", missing_values)
message("Infinite expression values: ", infinite_values)
message("NaN expression values: ", nan_values)
message("Negative TPM values: ", negative_values)

if (infinite_values > 0) {
  stop("Infinite expression values detected in selected breast models.")
}

if (nan_values > 0) {
  stop("NaN expression values detected; inspect the source data.")
}

if (negative_values > 0) {
  stop("Negative TPM values detected in selected breast models.")
}

# 5. Assemble the candidate, preserving all gene annotations.

candidate <- data.frame(
  gene_symbol = gene_symbol,
  ensembl_gene_id = ensembl_gene_id,
  gene_id = gene_id,
  breast_expression,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 6. Write validation diagnostics.

checks <- data.frame(
  check = c(
    "metadata_rows_have_equal_field_counts",
    "model_ids_are_unique",
    "expression_columns_match_selected_models",
    "all_breast_model_ids_found",
    "expression_columns_are_numeric",
    "no_infinite_expression_values",
    "no_nan_expression_values",
    "expression_values_are_nonnegative_or_missing",
    "gene_symbols_are_present",
    "missing_expression_values_reported"
  ),
  passed = c(
    length(unique(field_counts)) == 1L,
    !anyDuplicated(model_ids),
    ncol(breast_expression) == nrow(sources),
    all(sources$model_id %in% names(expr_all)),
    all(vapply(breast_expression, is.numeric, logical(1))),
    infinite_values == 0L,
    nan_values == 0L,
    all(expression_matrix >= 0 | is.na(expression_matrix)),
    !anyNA(gene_symbol) && all(gene_symbol != ""),
    !is.na(missing_values)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  checks,
  "results/script07_validation_checks.csv",
  row.names = FALSE
)

write.csv(
  model_metadata,
  "results/script07_all_model_metadata.csv",
  row.names = FALSE
)

write.csv(
  data.frame(model_id = sources$model_id),
  "results/script07_matched_model_ids.csv",
  row.names = FALSE
)

validation_summary <- data.frame(
  metadata_fields_per_row = field_counts[1],
  source_cell_lines = length(model_ids),
  expression_genes = nrow(candidate),
  matched_breast_models = ncol(candidate) - 3L,
  missing_expression_values = missing_values,
  infinite_expression_values = infinite_values,
  nan_expression_values = nan_values,
  negative_expression_values = negative_values,
  duplicated_gene_symbols = sum(duplicated(gene_symbol)),
  all_checks_passed = all(checks$passed)
)

write.csv(
  validation_summary,
  "results/script07_validation_summary.csv",
  row.names = FALSE
)

if (!all(checks$passed)) {
  stop("Validation failed. Review script07_validation_checks.csv.")
}

# 7. Compare against existing output without overwriting it.

candidate_matches_existing <- function(old, new) {
  if (!is.data.frame(old) ||
      !identical(dim(old), dim(new)) ||
      !identical(names(old), names(new))) {
    return(FALSE)
  }
  
  annotation_cols <- c(
    "gene_symbol", "ensembl_gene_id", "gene_id"
  )
  
  for (nm in annotation_cols) {
    if (!identical(as.character(old[[nm]]), as.character(new[[nm]]))) {
      return(FALSE)
    }
  }
  
  expression_cols <- setdiff(names(new), annotation_cols)
  
  for (nm in expression_cols) {
    if (!is.numeric(old[[nm]]) ||
        !isTRUE(all.equal(
          old[[nm]],
          new[[nm]],
          tolerance = 1e-10,
          check.attributes = FALSE
        ))) {
      return(FALSE)
    }
  }
  
  TRUE
}

if (file.exists(out_rds)) {
  existing <- readRDS(out_rds)
  
  if (candidate_matches_existing(existing, candidate)) {
    message(
      "PASS: candidate matches existing RDS. Original outputs preserved."
    )
  } else {
    saveRDS(candidate, candidate_rds)
    write.csv(candidate, candidate_csv, row.names = FALSE)
    
    message(
      "Candidate differs from existing output. Original files preserved. ",
      "Review ", candidate_rds, " before replacing anything."
    )
  }
} else {
  saveRDS(candidate, out_rds)
  write.csv(candidate, out_csv, row.names = FALSE)
  
  message("Created standard expression outputs.")
}

print(validation_summary)
print(checks)


