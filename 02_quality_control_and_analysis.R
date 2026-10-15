# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 02 - repository and analysis quality control
# Author: Kazhal Ghiasi
# Purpose:
#   1. Inventory existing results without deleting or modifying them.
#   2. Check expected files across the TCGA-BRCA/GDSC2 workflow.
#   3. Inspect the saved paired DESeq2 model and sample pairing.
#   4. Record missing files and outstanding QC concerns.
# This script does NOT download data, refit DESeq2, or overwrite
# the original model or differential-expression results.

# 1. Confirm the project directory

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

cat("\n=== PROJECT DIRECTORY ===\n")
cat(project_dir, "\n")

if (!dir.exists("results")) {
  stop(
    "The results/ directory was not found. ",
    "Open the RStudio project from the repository root first."
  )
}

# 2. Inventory existing results

result_files <- list.files(
  "results",
  recursive = TRUE,
  all.files = TRUE,
  full.names = FALSE
)

result_inventory <- data.frame(
  file = result_files,
  size_bytes = file.info(
    file.path("results", result_files)
  )$size,
  stringsAsFactors = FALSE
)

cat("\n=== RESULTS INVENTORY ===\n")
cat("Files found:", nrow(result_inventory), "\n")
print(result_inventory, row.names = FALSE)

# Save a diagnostic inventory, not the underlying data.
write.csv(
  result_inventory,
  "results/script02_results_inventory.csv",
  row.names = FALSE
)

# 3. Check important workflow files

expected_files <- c(
  # TCGA-BRCA and paired differential expression
  "results/tcga_brca_counts.rds",
  "results/tcga_brca_sample_metadata.csv",
  "results/deseq2_paired_113.rds",
  "results/deseq2_paired_113_results.rds",
  "results/deseq2_paired_113_results.csv",
  
  # TCGA tumour modelling and subtype signatures
  "results/ml_tumour_vst_matrix_1095.rds",
  "results/ml_tumour_metadata_1095.csv",
  "results/subtype_signatures_effect_filtered_top100.csv",
  
  # GDSC2 expression and drug response
  "results/gdsc2_breast_cell_line_tpm_51.rds",
  "results/gdsc2_breast_model_sources.csv",
  "results/gdsc2_breast_response_records.rds",
  "results/gdsc2_breast_response_records_matched_expression.rds",
  
  # Drug association and robustness analyses
  "results/gdsc2_subtype_drug_associations.csv",
  "results/gdsc2_drug_ranking_robustness.csv",
  "results/broad_only_drug_associations.csv"
)

file_check <- data.frame(
  file = expected_files,
  exists = file.exists(expected_files),
  stringsAsFactors = FALSE
)

cat("\n=== EXPECTED FILE CHECK ===\n")
print(file_check, row.names = FALSE)

write.csv(
  file_check,
  "results/script02_expected_file_check.csv",
  row.names = FALSE
)

missing_files <- file_check$file[!file_check$exists]

if (length(missing_files) > 0L) {
  cat("\nFiles not currently found:\n")
  cat(paste0(" - ", missing_files), sep = "\n")
  cat(
    "\nA missing file is a diagnostic finding, not automatically ",
    "an error. Some files may be created by later scripts.\n",
    sep = ""
  )
} else {
  cat("\nAll checked workflow files were found.\n")
}

# 4. Inspect the saved paired DESeq2 model

model_file <- "results/deseq2_paired_113.rds"
paired_results_file <- "results/deseq2_paired_113_results.rds"

paired_qc <- list(
  model_found = file.exists(model_file),
  results_found = file.exists(paired_results_file),
  model_dimensions = NULL,
  design = NULL,
  condition_counts = NULL,
  unique_patients = NA_integer_,
  pairing_valid = NA,
  notes = character()
)

if (file.exists(model_file)) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    paired_qc$notes <- c(
      paired_qc$notes,
      "DESeq2 is not installed; model inspection was skipped."
    )
  } else {
    
    dds <- readRDS(model_file)
    
    if (!inherits(dds, "DESeqDataSet")) {
      paired_qc$notes <- c(
        paired_qc$notes,
        "Saved paired model is not a DESeqDataSet."
      )
    } else {
      
      sample_metadata <- as.data.frame(
        SummarizedExperiment::colData(dds)
      )
      
      paired_qc$model_dimensions <- dim(dds)
      paired_qc$design <- paste(
        deparse(DESeq2::design(dds)),
        collapse = ""
      )
      
      cat("\n=== SAVED PAIRED MODEL ===\n")
      cat("Dimensions (genes x samples):\n")
      print(dim(dds))
      
      cat("\nDesign:\n")
      print(DESeq2::design(dds))
      
      if (all(c("patient", "condition") %in%
              names(sample_metadata))) {
        
        cat("\nSamples by condition:\n")
        print(table(
          sample_metadata$condition,
          useNA = "ifany"
        ))
        
        paired_qc$condition_counts <- paste(
          capture.output(print(table(
            sample_metadata$condition,
            useNA = "ifany"
          ))),
          collapse = " "
        )
        
        pair_table <- table(
          patient = sample_metadata$patient,
          condition = sample_metadata$condition,
          useNA = "ifany"
        )
        
        # Exactly two observed condition levels, one sample from
        # each condition per patient, and no missing pair labels.
        pairing_valid <-
          ncol(pair_table) == 2L &&
          all(rowSums(pair_table) == 2L) &&
          all(rowSums(pair_table == 1L) == 2L) &&
          !anyNA(sample_metadata$patient) &&
          !anyNA(sample_metadata$condition)
        
        paired_qc$unique_patients <- length(
          unique(sample_metadata$patient)
        )
        paired_qc$pairing_valid <- pairing_valid
        
        cat("\nUnique patients:",
            paired_qc$unique_patients, "\n")
        cat("Valid one-normal/one-tumour pairing:",
            pairing_valid, "\n")
        
        if (!pairing_valid) {
          paired_qc$notes <- c(
            paired_qc$notes,
            "Review patient-condition pairing before interpretation."
          )
        }
        
      } else {
        paired_qc$notes <- c(
          paired_qc$notes,
          "Saved model lacks patient and/or condition metadata."
        )
      }
      
      if (!identical(
        paired_qc$design,
        "~ patient + condition"
      )) {
        cat(
          "\nNote: inspect the displayed design carefully; ",
          "it may be equivalent but formatted differently.\n",
          sep = ""
        )
      }
      
      rm(dds)
    }
  }
  
} else {
  paired_qc$notes <- c(
    paired_qc$notes,
    "Saved paired model was not found; no model was refitted."
  )
}

# 5. Check the saved paired results

if (file.exists(paired_results_file)) {
  
  saved_results <- readRDS(paired_results_file)
  
  cat("\n=== SAVED PAIRED RESULTS ===\n")
  cat("Class:\n")
  print(class(saved_results))
  
  cat("Dimensions:\n")
  print(dim(saved_results))
  
  if ("padj" %in% names(as.data.frame(saved_results))) {
    results_df <- as.data.frame(saved_results)
    
    cat(
      "Genes with non-missing adjusted p-values:",
      sum(!is.na(results_df$padj)), "\n"
    )
    
    cat(
      "Genes with adjusted p-value < 0.05:",
      sum(results_df$padj < 0.05, na.rm = TRUE), "\n"
    )
  }
  
  rm(saved_results)
}

# 6. Record QC findings

cat("\n=== OUTSTANDING QC ITEMS ===\n")
cat(
  "- Investigate repeated TCGA tumour sample identifiers.\n",
  "- Verify count-matrix and metadata alignment before refitting.\n",
  "- Confirm the paired cohort contains 113 complete patient pairs.\n",
  "- Confirm downstream scripts use the actual saved filenames.\n",
  "- Distinguish missing expected files from analyses not yet run.\n",
  sep = ""
)

summary_lines <- c(
  "Script 02: repository and analysis QC",
  paste("Project directory:", project_dir),
  paste("Files in results/:", nrow(result_inventory)),
  paste("Missing expected files:", length(missing_files)),
  paste("Paired model found:", paired_qc$model_found),
  paste("Paired results found:", paired_qc$results_found),
  paste("Paired model dimensions:",
        if (is.null(paired_qc$model_dimensions)) {
          "not inspected"
        } else {
          paste(paired_qc$model_dimensions, collapse = " x ")
        }),
  paste("Unique patients:", paired_qc$unique_patients),
  paste("Pairing valid:", paired_qc$pairing_valid),
  "No model refitting was performed."
)

if (length(paired_qc$notes) > 0L) {
  summary_lines <- c(
    summary_lines,
    "",
    "Notes:",
    paste0("- ", paired_qc$notes)
  )
}

writeLines(
  summary_lines,
  "results/script02_qc_summary.txt"
)

cat("\n=== SCRIPT 02 COMPLETE ===\n")
cat("Diagnostic files saved in results/.\n")
cat("Original analysis objects were not modified.\n")
