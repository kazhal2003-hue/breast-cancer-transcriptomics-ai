
# ============================================================
# 05_prepare_gdsc2_expression.R
# Prepare TPM expression data for matched GDSC2 breast models
# ============================================================

# 1. Set the project directory
# Use the active RStudio project directory.
project_dir <- getwd()

# 2. Define input files
rna_zip <- "data/drug_response/rnaseq_merged_20260323.zip"
gdsc_file <- "data/drug_response/GDSC2_fitted_dose_response_27Oct23.xlsx"
tpm_member <- "rnaseq_merged_rsem_tpm_20260323.csv"

# 3. Check that the required files exist
stopifnot(file.exists(rna_zip), file.exists(gdsc_file))

dir.create("results", showWarnings = FALSE)

# 4. Load GDSC2 and identify breast cancer models
gdsc <- readxl::read_excel(gdsc_file)

breast_gdsc <- gdsc[
  grepl("breast", gdsc$CANCER_TYPE, ignore.case = TRUE),
]

gdsc_breast_ids <- unique(
  as.character(breast_gdsc$SANGER_MODEL_ID)
)

gdsc_breast_ids <- gdsc_breast_ids[
  !is.na(gdsc_breast_ids) & nzchar(gdsc_breast_ids)
]

# 5. Read the TPM model-ID header directly from the ZIP
con <- unz(rna_zip, tpm_member, open = "r")
header_line <- readLines(con, n = 1)
close(con)

header_fields <- strsplit(
  header_line, ",", fixed = TRUE
)[[1]]

tpm_model_ids <- header_fields[-(1:3)]

# 6. Match model IDs
matched_ids <- intersect(gdsc_breast_ids, tpm_model_ids)
missing_ids <- setdiff(gdsc_breast_ids, tpm_model_ids)

# Preserve the column order in the original TPM matrix
selected_ids <- tpm_model_ids[tpm_model_ids %in% matched_ids]

cat("GDSC2 breast models:", length(gdsc_breast_ids), "\n")
cat("Matched models:", length(selected_ids), "\n")
cat("Unmatched models:\n")
print(missing_ids)

stopifnot(length(selected_ids) > 0)

# 7. Read only the three gene-annotation columns and
#    the matched expression columns.
#    Unselected columns are skipped by R.
col_classes <- rep("NULL", length(tpm_model_ids) + 3)

col_classes[1:3] <- "character"

selected_positions <- which(
  tpm_model_ids %in% selected_ids
) + 3

col_classes[selected_positions] <- "numeric"

con <- unz(rna_zip, tpm_member, open = "r")

tpm_subset <- read.csv(
  con,
  skip = 3,
  header = TRUE,
  check.names = FALSE,
  colClasses = col_classes,
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

close(con)

# 8. Give columns clear, unique names
names(tpm_subset) <- c(
  "gene_symbol",
  "ensembl_gene_id",
  "gene_id",
  selected_ids
)

# 9. Check the extracted data
cat("\nExpression matrix dimensions:\n")
print(dim(tpm_subset))

cat("\nFirst six rows:\n")
print(tpm_subset[1:min(6, nrow(tpm_subset)), 1:min(8, ncol(tpm_subset))])

cat("\nMissing expression values:\n")
print(sum(is.na(as.matrix(tpm_subset[, selected_ids, drop = FALSE]))))

cat("\nDuplicate gene identifiers:\n")
print(sum(duplicated(tpm_subset$ensembl_gene_id)))

# 10. Save the reusable expression dataset
saveRDS(
  tpm_subset,
  "results/gdsc2_breast_cell_line_tpm_51.rds"
)

write.csv(
  tpm_subset,
  "results/gdsc2_breast_cell_line_tpm_51.csv",
  row.names = FALSE,
  na = ""
)

# Save the matched and unmatched model IDs for traceability
write.csv(
  data.frame(SANGER_MODEL_ID = selected_ids),
  "results/gdsc2_breast_matched_model_ids.csv",
  row.names = FALSE
)

write.csv(
  data.frame(SANGER_MODEL_ID = missing_ids),
  "results/gdsc2_breast_unmatched_model_ids.csv",
  row.names = FALSE
)

cat("\nDONE: expression data saved in the results folder.\n")
