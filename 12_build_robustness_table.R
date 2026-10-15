# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 12 - build robustness table 
# Author: Kazhal Ghiasi
# Notes: Compares drug-response associations from the full matched
# breast cancer cell-line panel with those from the Broad-only
# sensitivity analysis.
# Primary response metric: LN_IC50
# Primary drug-coverage threshold: >= 20 unique cell lines
# Comparison measures:
#   - Spearman correlation between association rankings
#   - Overlap between the top-ranked drugs
#   - Concordance of association directions
#   - Agreement in FDR-based statistical significance
# This is a sensitivity analysis, not evidence of clinical efficacy.

# 0. Packages and input validation 

required_packages <- c("dplyr", "readr", "tibble")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(readr)
library(tibble)

results_dir <- "results"

original_file <- file.path(
  results_dir,
  "gdsc2_subtype_drug_associations.csv"
)

broad_file <- file.path(
  results_dir,
  "broad_only_drug_associations.csv"
)

if (!dir.exists(results_dir)) {
  stop("The results/ directory does not exist.")
}

missing_files <- c(original_file, broad_file)[
  !file.exists(c(original_file, broad_file))
]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) not found:\n",
    paste(missing_files, collapse = "\n"),
    "\nRun the relevant upstream scripts first."
  )
}

# 1. Load association results 

original <- read_csv(
  original_file,
  show_col_types = FALSE
)

broad <- read_csv(
  broad_file,
  show_col_types = FALSE
)

# 2. Standardise and validate column names 

standardise_associations <- function(data, dataset_name) {
  
  # Accept either response_metric or metric.
  if (!"metric" %in% names(data)) {
    if ("response_metric" %in% names(data)) {
      data <- data %>%
        rename(metric = response_metric)
    } else {
      stop(
        dataset_name,
        " results must contain either 'metric' or ",
        "'response_metric'."
      )
    }
  }
  
  required_columns <- c(
    "DRUG_ID",
    "subtype_comparison",
    "metric",
    "n_models",
    "rho",
    "FDR"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      dataset_name,
      " results are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Preserve drug names where available.
  if (!"DRUG_NAME" %in% names(data)) {
    data$DRUG_NAME <- NA_character_
  }
  
  data <- data %>%
    mutate(
      DRUG_ID = as.character(DRUG_ID),
      subtype_comparison = as.character(subtype_comparison),
      metric = as.character(metric),
      n_models = as.numeric(n_models),
      rho = as.numeric(rho),
      FDR = as.numeric(FDR),
      DRUG_NAME = as.character(DRUG_NAME)
    )
  
  if (any(is.na(data$DRUG_ID) | data$DRUG_ID == "")) {
    stop(dataset_name, " results contain missing drug IDs.")
  }
  
  if (any(
    is.na(data$subtype_comparison) |
    data$subtype_comparison == ""
  )) {
    stop(
      dataset_name,
      " results contain missing subtype comparisons."
    )
  }
  
  data
}

original <- standardise_associations(
  original,
  "Full-panel"
)

broad <- standardise_associations(
  broad,
  "Broad-only"
)

# 3. Select the primary analysis

primary_metric <- "LN_IC50"
primary_threshold <- 20

original_primary <- original %>%
  filter(
    metric == primary_metric,
    n_models >= primary_threshold,
    is.finite(rho)
  )

broad_primary <- broad %>%
  filter(
    metric == primary_metric,
    n_models >= primary_threshold,
    is.finite(rho)
  )

if (nrow(original_primary) == 0) {
  stop(
    "No eligible full-panel LN_IC50 associations remain ",
    "at the >= 20-cell-line threshold."
  )
}

if (nrow(broad_primary) == 0) {
  stop(
    "No eligible Broad-only LN_IC50 associations remain ",
    "at the >= 20-cell-line threshold."
  )
}

# A drug should appear only once per subtype and metric.
check_duplicate_associations <- function(data, dataset_name) {
  
  duplicates <- data %>%
    count(
      subtype_comparison,
      metric,
      DRUG_ID,
      name = "records"
    ) %>%
    filter(records > 1)
  
  if (nrow(duplicates) > 0) {
    write_csv(
      duplicates,
      file.path(
        results_dir,
        paste0(
          "script12_",
          gsub("[^A-Za-z0-9]+", "_", dataset_name),
          "_duplicate_association_audit.csv"
        )
      )
    )
    
    stop(
      dataset_name,
      " results contain duplicate drug/subtype/metric records. ",
      "Inspect the duplicate-association audit before proceeding."
    )
  }
}

check_duplicate_associations(
  original_primary,
  "full_panel"
)

check_duplicate_associations(
  broad_primary,
  "broad_only"
)

# 4. Identify shared subtype signatures 

signatures <- intersect(
  unique(original_primary$subtype_comparison),
  unique(broad_primary$subtype_comparison)
)

if (length(signatures) == 0) {
  stop(
    "No shared subtype comparisons were found between ",
    "the full-panel and Broad-only results."
  )
}

message(
  "Shared subtype comparisons: ",
  paste(signatures, collapse = ", ")
)

# 5. Compare ranking robustness by subtype 

robustness_results <- lapply(signatures, function(sig) {
  
  orig_sig <- original_primary %>%
    filter(subtype_comparison == sig) %>%
    select(
      DRUG_ID,
      DRUG_NAME,
      rho_original = rho,
      FDR_original = FDR
    )
  
  broad_sig <- broad_primary %>%
    filter(subtype_comparison == sig) %>%
    select(
      DRUG_ID,
      rho_broad = rho,
      FDR_broad = FDR
    )
  
  # Compare only drugs tested in both analyses.
  shared <- inner_join(
    orig_sig,
    broad_sig,
    by = "DRUG_ID"
  ) %>%
    filter(
      is.finite(rho_original),
      is.finite(rho_broad)
    ) %>%
    mutate(
      sig_original = !is.na(FDR_original) &
        FDR_original < 0.05,
      sig_broad = !is.na(FDR_broad) &
        FDR_broad < 0.05
    )
  
  n_shared <- nrow(shared)
  
  if (n_shared == 0) {
    return(
      tibble(
        signature = sig,
        shared_drugs = 0L,
        spearman_ranking_correlation = NA_real_,
        top10_overlap_count = NA_integer_,
        top10_overlap_percent = NA_real_,
        direction_comparable_drugs = 0L,
        direction_concordance_percent = NA_real_,
        significant_in_both = 0L,
        significant_original_only = 0L,
        significant_broad_only = 0L,
        significant_in_neither = 0L
      )
    )
  }
  
  # Rank agreement across all shared drugs.
  ranking_rho <- if (
    n_shared >= 3 &&
    dplyr::n_distinct(shared$rho_original) > 1 &&
    dplyr::n_distinct(shared$rho_broad) > 1
  ) {
    cor(
      shared$rho_original,
      shared$rho_broad,
      method = "spearman"
    )
  } else {
    NA_real_
  }
  
  # Top drugs ranked by absolute association strength.
  top_original <- shared %>%
    arrange(
      desc(abs(rho_original)),
      DRUG_ID
    ) %>%
    slice_head(n = 10) %>%
    pull(DRUG_ID)
  
  top_broad <- shared %>%
    arrange(
      desc(abs(rho_broad)),
      DRUG_ID
    ) %>%
    slice_head(n = 10) %>%
    pull(DRUG_ID)
  
  top_overlap <- length(
    intersect(top_original, top_broad)
  )
  
  # Percentage uses the actual number of drugs ranked
  # in each top list, which may be fewer than ten.
  top_list_size <- length(top_original)
  
  top_overlap_percent <- if (top_list_size > 0) {
    100 * top_overlap / top_list_size
  } else {
    NA_real_
  }
  
  # Direction concordance excludes exact zero correlations.
  nonzero <- shared %>%
    filter(
      rho_original != 0,
      rho_broad != 0
    )
  
  n_direction_comparable <- nrow(nonzero)
  
  direction_concordance <- if (
    n_direction_comparable > 0
  ) {
    mean(
      sign(nonzero$rho_original) ==
        sign(nonzero$rho_broad)
    )
  } else {
    NA_real_
  }
  
  # FDR significance agreement among shared drugs.
  tibble(
    signature = sig,
    shared_drugs = n_shared,
    spearman_ranking_correlation = ranking_rho,
    top10_overlap_count = as.integer(top_overlap),
    top10_overlap_percent = top_overlap_percent,
    direction_comparable_drugs = n_direction_comparable,
    direction_concordance_percent =
      100 * direction_concordance,
    significant_in_both = sum(
      shared$sig_original & shared$sig_broad
    ),
    significant_original_only = sum(
      shared$sig_original & !shared$sig_broad
    ),
    significant_broad_only = sum(
      !shared$sig_original & shared$sig_broad
    ),
    significant_in_neither = sum(
      !shared$sig_original & !shared$sig_broad
    )
  )
})

robustness_table <- bind_rows(robustness_results) %>%
  mutate(
    across(
      c(
        spearman_ranking_correlation,
        top10_overlap_percent,
        direction_concordance_percent
      ),
      ~ round(.x, 3)
    )
  )


# 6. Save the robustness summary 

output_file <- file.path(
  results_dir,
  "broad_only_robustness_summary.csv"
)

write_csv(
  robustness_table,
  output_file
)

# Save the settings alongside the summary for reproducibility.
details <- list(
  script = "12_build_robustness_table.R",
  author = "Kazhal Ghiasi",
  project = "From Cancer Transcriptomics to AI-Guided Drug Prioritization",
  comparison = "Full matched panel versus Broad-only sensitivity analysis",
  primary_metric = primary_metric,
  primary_minimum_cell_lines = primary_threshold,
  fdr_significance_threshold = 0.05,
  ranking_method = "Spearman correlation of association coefficients",
  top_drug_selection = "Top 10 by absolute Spearman rho among shared drugs",
  full_panel_input = original_file,
  broad_only_input = broad_file,
  shared_signatures = signatures,
  summary_output = output_file,
  results = robustness_table
)

saveRDS(
  details,
  file.path(
    results_dir,
    "broad_only_robustness_summary_details.rds"
  )
)

# 7. Report results 

cat("\n=== BROAD-ONLY ROBUSTNESS SUMMARY ===\n")
cat("Primary metric:", primary_metric, "\n")
cat("Minimum cell lines per drug:", primary_threshold, "\n\n")

print(
  robustness_table,
  n = Inf,
  width = Inf
)

cat(
  "\nSaved: ",
  output_file,
  "\nSaved: results/broad_only_robustness_summary_details.rds\n",
  sep = ""
)

cat("\n=== SCRIPT 12 COMPLETE ===\n")

