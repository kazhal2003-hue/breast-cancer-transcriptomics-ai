# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 02 - quality control and analysis
# Author: Kazhal Ghiasi
# Purpose: Inspect saved TCGA-BRCA metadata and analysis outputs.
# Status: Preliminary inspection; sample-level quality control remains pending.
# Notes:
#   - This script does not download data or rerun DESeq2.
#   - Missing files are reported for investigation.
#   - Repeated tumour sample identifiers require further review before
#     final interpretation of downstream analyses.

# 1. Confirm the project directory 

cat("\n=== PROJECT DIRECTORY ===\n")
cat(normalizePath(getwd()), "\n")

# 2. Check expected files 

expected_files <- c(
  "results/tcga_brca_counts.rds",
  "results/tcga_brca_sample_metadata.csv",
  "results/deseq2_results.rds",
  "results/deseq2_completed.rds",
  "results/qc_checkpoint_report.txt"
)

file_check <- data.frame(
  file = expected_files,
  exists = file.exists(expected_files),
  row.names = NULL
)

cat("\n=== EXPECTED FILE CHECK ===\n")
print(file_check)

missing_files <- file_check$file[!file_check$exists]

if (length(missing_files) > 0L) {
  cat("\nFiles not found:\n")
  cat(paste0(" - ", missing_files), sep = "\n")
  cat("\n\nCheck whether these files are expected at this stage ",
      "of the workflow.\n", sep = "")
} else {
  cat("\nAll expected files were found.\n")
}

# 3. Inspect saved sample metadata 

metadata_file <- "results/tcga_brca_sample_metadata.csv"

if (file.exists(metadata_file)) {
  
  sample_metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("\n=== SAMPLE METADATA ===\n")
  cat("Dimensions:", nrow(sample_metadata), "rows x",
      ncol(sample_metadata), "columns\n")
  cat("Columns:\n")
  print(names(sample_metadata))
  
  if ("condition" %in% names(sample_metadata)) {
    cat("\nSample counts by condition:\n")
    print(table(sample_metadata$condition, useNA = "ifany"))
  } else {
    cat(
      "\nWARNING: No 'condition' column was found. ",
      "Confirm the metadata schema before downstream analysis.\n",
      sep = ""
    )
  }
  
  if (anyDuplicated(names(sample_metadata))) {
    cat("\nWARNING: Duplicate metadata column names were detected.\n")
  }
  
} else {
  cat("\nMetadata inspection skipped: file not found.\n")
}

# 4. Inspect saved differential-expression results 

results_file <- "results/deseq2_results.rds"

if (file.exists(results_file)) {
  
  saved_results <- readRDS(results_file)
  
  cat("\n=== SAVED DIFFERENTIAL-EXPRESSION RESULTS ===\n")
  cat("Object class:\n")
  print(class(saved_results))
  
  result_dimensions <- dim(saved_results)
  
  if (!is.null(result_dimensions)) {
    cat("Dimensions:\n")
    print(result_dimensions)
  } else {
    cat("This object does not report tabular or array dimensions.\n")
    cat("Structure:\n")
    str(saved_results, max.level = 1L)
  }
  
  rm(saved_results)
  
} else {
  cat("\nDifferential-expression inspection skipped: file not found.\n")
}

# 5. Record outstanding quality-control concerns 

cat("\n=== OUTSTANDING QUALITY-CONTROL ITEMS ===\n")
cat(
  "- Investigate repeated tumour sample identifiers and their aliquots.\n",
  "- Verify sample selection and count-profile consistency.\n",
  "- Confirm count-matrix and metadata alignment before downstream analysis.\n",
  "- Perform sample-level QC before final interpretation.\n",
  sep = ""
)

# 6. Finish 

cat("\n=== DATA INSPECTION COMPLETE ===\n")
cat("No differential-expression analysis was rerun by this script.\n")
