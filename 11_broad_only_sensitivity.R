# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 11 - broad only sensitivity
# Author: Kazhal Ghiasi
# Notes: Recalculates breast cancer subtype-signature scores using
# Broad-source cell lines only and retests associations with
# GDSC2 drug-response measurements.
# This is a sensitivity analysis. Original results are preserved.
# Minimum coverage: >= 20 unique cell lines per drug.
# Associations are exploratory and do not establish clinical efficacy.

library(dplyr)
library(readr)

options(na.print = "NA")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. Load existing data

tpm_file <- "results/gdsc2_breast_cell_line_tpm_51.rds"

matched_response_file <-
  "results/gdsc2_breast_response_records_matched_expression.rds"

full_response_file <-
  "results/gdsc2_breast_response_records.rds"

signature_file <-
  "results/subtype_signatures_effect_filtered_top100.csv"

metadata_file <-
  "results/gdsc2_breast_model_sources.csv"

if (!file.exists(tpm_file)) {
  stop("Missing expression file: ", tpm_file)
}

if (!file.exists(signature_file)) {
  stop("Missing signature file: ", signature_file)
}

if (file.exists(matched_response_file)) {
  response_file <- matched_response_file
} else if (file.exists(full_response_file)) {
  warning(
    "Matched response file not found. Using full breast response data ",
    "and restricting it to the selected Broad cell lines."
  )
  response_file <- full_response_file
} else {
  stop("No GDSC2 response RDS file was found.")
}

tpm <- readRDS(tpm_file)
response <- readRDS(response_file)

signatures <- read_csv(
  signature_file,
  show_col_types = FALSE
)

required_tpm <- c("gene_symbol", "gene_id")

required_response <- c(
  "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME",
  "LN_IC50", "AUC"
)

required_signatures <- c(
  "gene_symbol", "logFC", "subtype_comparison"
)

if (!all(required_tpm %in% names(tpm))) {
  stop("Expression data are missing required annotation columns.")
}

if (!all(required_response %in% names(response))) {
  stop("Response data are missing required columns.")
}

if (!all(required_signatures %in% names(signatures))) {
  stop("Signature data are missing required columns.")
}


# 2. Load and validate verified cell-line source metadata

matched_metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

required_metadata <- c(
  "model_id",
  "data_source_matched"
)

if (!all(required_metadata %in% names(matched_metadata))) {
  stop(
    "Source metadata must contain model_id and data_source_matched."
  )
}

# Check conflicting duplicate IDs before deduplication.
metadata_check <- matched_metadata %>%
  mutate(
    model_id = trimws(as.character(model_id)),
    data_source_matched = trimws(
      as.character(data_source_matched)
    )
  ) %>%
  filter(!is.na(model_id), model_id != "")

metadata_conflicts <- metadata_check %>%
  filter(
    !is.na(data_source_matched),
    data_source_matched != ""
  ) %>%
  distinct(model_id, data_source_matched) %>%
  count(model_id, name = "n_sources") %>%
  filter(n_sources > 1)

if (nrow(metadata_conflicts) > 0) {
  write_csv(
    metadata_conflicts,
    "results/broad_only_metadata_conflicts.csv"
  )
  stop(
    "Conflicting source labels found for duplicate model IDs. ",
    "See results/broad_only_metadata_conflicts.csv."
  )
}

matched_metadata <- metadata_check %>%
  filter(
    !is.na(data_source_matched),
    data_source_matched != ""
  ) %>%
  distinct(model_id, .keep_all = TRUE) %>%
  transmute(
    model_id,
    data_source = data_source_matched
  )

if (anyDuplicated(matched_metadata$model_id)) {
  stop("Source metadata contain duplicate model IDs.")
}

cat(
  "Verified metadata model IDs:",
  nrow(matched_metadata), "\n"
)

cat("Source labels:\n")
print(table(matched_metadata$data_source))

# 3. Identify Broad cell lines

tpm_model_ids <- setdiff(
  names(tpm),
  c("gene_symbol", "ensembl_gene_id", "gene_id")
)

broad_metadata_ids <- matched_metadata %>%
  filter(tolower(data_source) == "broad") %>%
  pull(model_id) %>%
  unique()

broad_ids <- intersect(
  broad_metadata_ids,
  tpm_model_ids
)

cat("Broad cell lines in expression matrix:", length(broad_ids), "\n")

if (length(broad_ids) < 20) {
  stop(
    "Fewer than 20 Broad cell lines are available. ",
    "Check the verified source metadata and model-ID matching."
  )
}

write_csv(
  data.frame(SANGER_MODEL_ID = broad_ids),
  "results/broad_only_model_ids.csv"
)

# 4. Prepare Broad-only expression data

broad_expr_df <- tpm[, broad_ids, drop = FALSE]

broad_expr_df[] <- lapply(
  broad_expr_df,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

broad_expr <- as.matrix(broad_expr_df)

if (any(is.infinite(broad_expr), na.rm = TRUE)) {
  stop("Infinite TPM values were found in the Broad expression data.")
}

if (any(broad_expr < 0, na.rm = TRUE)) {
  stop("Negative TPM values were found in the Broad expression data.")
}

complete_genes <- rowSums(is.na(broad_expr)) == 0

broad_expr <- broad_expr[complete_genes, , drop = FALSE]

gene_symbols <- trimws(
  as.character(tpm$gene_symbol[complete_genes])
)

valid_symbols <- !is.na(gene_symbols) & gene_symbols != ""

broad_expr <- broad_expr[valid_symbols, , drop = FALSE]
gene_symbols <- gene_symbols[valid_symbols]

# Do not silently choose an arbitrary row when a signature gene
# maps to multiple expression rows.
duplicate_symbols <- unique(
  gene_symbols[
    duplicated(gene_symbols) |
      duplicated(gene_symbols, fromLast = TRUE)
  ]
)

signature_symbols <- unique(
  trimws(as.character(signatures$gene_symbol))
)

problem_symbols <- intersect(
  signature_symbols,
  duplicate_symbols
)

if (length(problem_symbols) > 0) {
  
  write_csv(
    data.frame(gene_symbol = problem_symbols),
    "results/broad_only_duplicate_signature_symbols.csv"
  )
  
  stop(
    "Signature genes map to duplicate expression rows. See ",
    "results/broad_only_duplicate_signature_symbols.csv."
  )
}

rownames(broad_expr) <- gene_symbols

cat("Complete, unique gene symbols:", nrow(broad_expr), "\n")

# Transform and standardise within the Broad-only cell-line panel.
log_expr <- log2(broad_expr + 1)

gene_means <- rowMeans(log_expr)
gene_sds <- apply(log_expr, 1, sd)

variable_genes <- is.finite(gene_sds) & gene_sds > 0

log_expr <- log_expr[variable_genes, , drop = FALSE]
gene_means <- gene_means[variable_genes]
gene_sds <- gene_sds[variable_genes]

z_expr <- sweep(log_expr, 1, gene_means, "-")
z_expr <- sweep(z_expr, 1, gene_sds, "/")

# 5. Recalculate subtype-signature scores

signature_names <- c(
  "Basal_vs_rest",
  "Her2_vs_rest",
  "LumA_vs_rest",
  "LumB_vs_rest"
)

missing_signatures <- setdiff(
  signature_names,
  unique(signatures$subtype_comparison)
)

if (length(missing_signatures) > 0) {
  stop(
    "Missing signature comparisons: ",
    paste(missing_signatures, collapse = ", ")
  )
}

score_list <- list()
overlap_list <- list()

for (sig_name in signature_names) {
  
  sig <- signatures %>%
    filter(
      subtype_comparison == sig_name,
      !is.na(gene_symbol),
      gene_symbol != "",
      is.finite(logFC)
    )
  
  # A gene must not contribute to both directions.
  sig <- sig %>%
    distinct(gene_symbol, .keep_all = TRUE)
  
  up_genes <- intersect(
    sig$gene_symbol[sig$logFC > 0],
    rownames(z_expr)
  )
  
  down_genes <- intersect(
    sig$gene_symbol[sig$logFC < 0],
    rownames(z_expr)
  )
  
  if (length(up_genes) < 10 || length(down_genes) < 10) {
    stop("Too few mapped genes for signature: ", sig_name)
  }
  
  up_score <- colMeans(
    z_expr[up_genes, , drop = FALSE]
  )
  
  down_score <- colMeans(
    z_expr[down_genes, , drop = FALSE]
  )
  
  score_list[[sig_name]] <- up_score - down_score
  
  overlap_list[[sig_name]] <- data.frame(
    subtype_comparison = sig_name,
    up_genes_available = length(up_genes),
    down_genes_available = length(down_genes),
    total_genes_available = length(union(up_genes, down_genes))
  )
}

broad_scores <- data.frame(
  SANGER_MODEL_ID = colnames(z_expr),
  as.data.frame(score_list, check.names = FALSE),
  check.names = FALSE
)

broad_overlap <- bind_rows(overlap_list)

write_csv(
  broad_scores,
  "results/broad_only_subtype_scores.csv"
)

write_csv(
  broad_overlap,
  "results/broad_only_signature_gene_overlap.csv"
)

saveRDS(
  broad_scores,
  "results/broad_only_subtype_scores.rds"
)

cat("\nBroad-only signature scores saved.\n")
print(head(broad_scores))
print(broad_overlap)

# 6. Restrict response data to Broad cell lines

response$SANGER_MODEL_ID <-
  trimws(as.character(response$SANGER_MODEL_ID))

response$DRUG_ID <- as.character(response$DRUG_ID)

broad_response <- response %>%
  filter(
    SANGER_MODEL_ID %in% broad_ids,
    !is.na(DRUG_ID),
    !is.na(SANGER_MODEL_ID)
  ) %>%
  inner_join(
    broad_scores,
    by = "SANGER_MODEL_ID"
  )

cat(
  "\nBroad-only drug-response records:",
  nrow(broad_response), "\n"
)

cat(
  "Broad-only cell lines with response data:",
  n_distinct(broad_response$SANGER_MODEL_ID), "\n"
)

# Audit uniqueness after restricting to the selected models.
duplicate_pairs <- broad_response %>%
  count(SANGER_MODEL_ID, DRUG_ID, name = "n_records") %>%
  filter(n_records > 1)

if (nrow(duplicate_pairs) > 0) {
  
  write_csv(
    duplicate_pairs,
    "results/broad_only_duplicate_response_pairs.csv"
  )
  
  stop(
    "Repeated Broad cell-line/drug pairs were found. See ",
    "results/broad_only_duplicate_response_pairs.csv. ",
    "No response records were silently removed."
  )
}

# 7. Calculate drug-response associations
# Primary metric: LN_IC50.
# Secondary metric: AUC.
# Minimum coverage: 20 unique cell lines per drug.

metrics <- c("LN_IC50", "AUC")

association_results <- list()
result_index <- 1

for (sig_name in signature_names) {
  
  for (metric_name in metrics) {
    
    for (drug_id in unique(broad_response$DRUG_ID)) {
      
      drug_data <- broad_response %>%
        filter(DRUG_ID == drug_id) %>%
        select(
          SANGER_MODEL_ID,
          DRUG_NAME,
          all_of(metric_name),
          all_of(sig_name)
        )
      
      names(drug_data)[names(drug_data) == sig_name] <-
        "signature_score"
      
      drug_data[[metric_name]] <- suppressWarnings(
        as.numeric(drug_data[[metric_name]])
      )
      
      drug_data <- drug_data %>%
        filter(
          is.finite(.data[[metric_name]]),
          is.finite(signature_score)
        )
      
      n_models <- n_distinct(drug_data$SANGER_MODEL_ID)
      
      if (n_models < 20) next
      
      if (
        sd(drug_data[[metric_name]]) == 0 ||
        sd(drug_data$signature_score) == 0
      ) {
        next
      }
      
      test_result <- suppressWarnings(
        cor.test(
          drug_data$signature_score,
          drug_data[[metric_name]],
          method = "spearman",
          exact = FALSE
        )
      )
      
      association_results[[result_index]] <- data.frame(
        subtype_comparison = sig_name,
        DRUG_ID = drug_id,
        DRUG_NAME = as.character(drug_data$DRUG_NAME[1]),
        metric = metric_name,
        n_models = n_models,
        rho = unname(test_result$estimate),
        p_value = test_result$p.value
      )
      
      result_index <- result_index + 1
    }
  }
}

if (length(association_results) == 0) {
  
  broad_associations <- data.frame(
    subtype_comparison = character(),
    DRUG_ID = character(),
    DRUG_NAME = character(),
    metric = character(),
    n_models = integer(),
    rho = double(),
    p_value = double(),
    FDR = double(),
    rank_by_abs_rho = integer()
  )
  
  warning(
    "No Broad-only associations met the >=20-cell-line threshold. ",
    "Check coverage before interpreting this sensitivity analysis."
  )
  
} else {
  
  broad_associations <- bind_rows(association_results) %>%
    group_by(subtype_comparison, metric) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    group_by(subtype_comparison, metric) %>%
    mutate(
      rank_by_abs_rho = rank(-abs(rho), ties.method = "min")
    ) %>%
    ungroup()
}

write_csv(
  broad_associations,
  "results/broad_only_drug_associations.csv"
)

# 8. Summarise the strongest results

cat(
  "\nNumber of Broad-only association tests:",
  nrow(broad_associations), "\n"
)

cat("\nBroad-only eligible drug counts by metric and signature:\n")

if (nrow(broad_associations) > 0) {
  
  print(
    broad_associations %>%
      group_by(subtype_comparison, metric) %>%
      summarise(
        eligible_drugs = n_distinct(DRUG_ID),
        .groups = "drop"
      )
  )
  
  cat("\nTop LN_IC50 associations by signature:\n")
  
  top_results <- broad_associations %>%
    filter(metric == "LN_IC50") %>%
    group_by(subtype_comparison) %>%
    arrange(p_value, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    select(
      subtype_comparison,
      DRUG_NAME,
      n_models,
      rho,
      p_value,
      FDR
    )
  
  print(top_results, n = Inf)
  
} else {
  
  cat(
    "No eligible associations were produced. ",
    "This is a coverage limitation, not evidence of no association.\n"
  )
}

cat("\nBroad-only sensitivity analysis complete.\n")
