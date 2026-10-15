# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 09 - score signatures and test drug response 
# Author: Kazhal Ghiasi
# Score breast cancer subtype signatures and test GDSC2 drug response.
# Primary ranking coverage: >= 20 unique cell lines per drug.
# Sensitivity analyses: >= 15 and >= 25 unique cell lines.
# Repeated response records are audited, never silently averaged or removed.

library(dplyr)
library(readr)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. Load and validate inputs

tpm_file <- "results/gdsc2_breast_cell_line_tpm_51.rds"

# Prefer the expression-matched response file produced by Script 08.
matched_response_file <-
  "results/gdsc2_breast_response_records_matched_expression.rds"

full_response_file <-
  "results/gdsc2_breast_response_records.rds"

signature_file <-
  "results/subtype_signatures_effect_filtered_top100.csv"

if (!file.exists(tpm_file)) {
  stop("Missing expression input: ", tpm_file)
}

if (file.exists(matched_response_file)) {
  response_file <- matched_response_file
} else if (file.exists(full_response_file)) {
  warning(
    "Matched response file not found. Using the full breast response ",
    "dataset and matching its model IDs to the expression data."
  )
  response_file <- full_response_file
} else {
  stop("No GDSC2 response RDS file was found.")
}

if (!file.exists(signature_file)) {
  stop("Missing signature input: ", signature_file)
}

tpm <- readRDS(tpm_file)
response <- readRDS(response_file)
signatures <- read_csv(signature_file, show_col_types = FALSE)

required_tpm <- c("gene_symbol", "gene_id")

required_response <- c(
  "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME",
  "LN_IC50", "AUC", "RMSE"
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

# Normalise identifiers and remove invalid signature rows.
tpm$gene_symbol <- trimws(as.character(tpm$gene_symbol))

response$SANGER_MODEL_ID <-
  trimws(as.character(response$SANGER_MODEL_ID))

signatures$gene_symbol <-
  trimws(as.character(signatures$gene_symbol))

signatures$subtype_comparison <-
  as.character(signatures$subtype_comparison)

signatures$logFC <- suppressWarnings(
  as.numeric(signatures$logFC)
)

signatures <- signatures %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    is.finite(logFC),
    !is.na(subtype_comparison)
  )

subtype_names <- c(
  "Basal_vs_rest",
  "Her2_vs_rest",
  "LumA_vs_rest",
  "LumB_vs_rest"
)

missing_subtypes <- setdiff(
  subtype_names,
  unique(signatures$subtype_comparison)
)

if (length(missing_subtypes) > 0) {
  stop(
    "Missing signature comparisons: ",
    paste(missing_subtypes, collapse = ", ")
  )
}

# 2. Match expression samples and response records

annotation_columns <- c(
  "gene_symbol", "ensembl_gene_id", "gene_id"
)

model_ids <- setdiff(names(tpm), annotation_columns)

model_ids <- intersect(
  model_ids,
  unique(response$SANGER_MODEL_ID)
)

if (length(model_ids) < 25) {
  stop(
    "Only ", length(model_ids),
    " expression-matched cell lines were found. ",
    "Check model identifiers and Script 07/08 outputs."
  )
}

message("Expression-matched cell lines: ", length(model_ids))

# Keep only rows with complete expression measurements across
# the matched cell lines, as in the original analysis.
expression_data <- tpm[, model_ids, drop = FALSE]

expression_data[] <- lapply(
  expression_data,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

expression_values <- as.matrix(expression_data)

if (any(is.infinite(expression_values), na.rm = TRUE)) {
  stop("Infinite expression values were found.")
}

if (any(expression_values < 0, na.rm = TRUE)) {
  stop("Negative TPM values were found.")
}

complete_mask <- rowSums(is.na(expression_values)) == 0

complete_tpm <- tpm[complete_mask, , drop = FALSE]
expression_values <- expression_values[complete_mask, , drop = FALSE]

gene_symbols <- trimws(as.character(complete_tpm$gene_symbol))

valid_symbols <- !is.na(gene_symbols) & gene_symbols != ""

expression_values <- expression_values[valid_symbols, , drop = FALSE]
gene_symbols <- gene_symbols[valid_symbols]

# Log-transform TPM.
expression_matrix <- log2(expression_values + 1)
rownames(expression_matrix) <- gene_symbols

# Do not arbitrarily select one row if a gene symbol occurs more
# than once in the measured expression data.
duplicate_symbols <- unique(
  gene_symbols[duplicated(gene_symbols) |
                 duplicated(gene_symbols, fromLast = TRUE)]
)

problem_symbols <- intersect(
  unique(signatures$gene_symbol),
  duplicate_symbols
)

if (length(problem_symbols) > 0) {
  write_csv(
    data.frame(gene_symbol = problem_symbols),
    "results/gdsc2_signature_duplicate_expression_symbols.csv"
  )
  
  stop(
    "Some signature genes map to multiple expression rows. ",
    "See results/gdsc2_signature_duplicate_expression_symbols.csv"
  )
}

# 3. Standardise gene expression across matched cell lines

gene_means <- rowMeans(expression_matrix)
gene_sds <- apply(expression_matrix, 1, sd)

valid_genes <- is.finite(gene_sds) & gene_sds > 0

expression_matrix <- expression_matrix[valid_genes, , drop = FALSE]
gene_means <- gene_means[valid_genes]
gene_sds <- gene_sds[valid_genes]

z_matrix <- sweep(expression_matrix, 1, gene_means, "-")
z_matrix <- sweep(z_matrix, 1, gene_sds, "/")

# 4. Calculate subtype signature scores

score_list <- list()
overlap_list <- list()

for (subtype in subtype_names) {
  
  sig <- signatures %>%
    filter(subtype_comparison == subtype)
  
  up_genes <- unique(sig$gene_symbol[sig$logFC > 0])
  down_genes <- unique(sig$gene_symbol[sig$logFC < 0])
  
  up_genes <- intersect(up_genes, rownames(z_matrix))
  down_genes <- intersect(down_genes, rownames(z_matrix))
  
  if (length(up_genes) < 10 || length(down_genes) < 10) {
    stop(
      "Fewer than 10 available genes in one signature direction for ",
      subtype,
      ". See the signature overlap results after inspecting gene mapping."
    )
  }
  
  # Equal weighting within each direction.
  # Positive score = relatively stronger upregulated than
  # downregulated signature expression in this cell-line panel.
  up_score <- colMeans(z_matrix[up_genes, , drop = FALSE])
  down_score <- colMeans(z_matrix[down_genes, , drop = FALSE])
  
  score_list[[subtype]] <- up_score - down_score
  
  overlap_list[[subtype]] <- data.frame(
    subtype_comparison = subtype,
    up_genes_available = length(up_genes),
    down_genes_available = length(down_genes),
    total_genes_available = length(union(up_genes, down_genes))
  )
}

score_table <- data.frame(
  SANGER_MODEL_ID = colnames(z_matrix),
  as.data.frame(score_list, check.names = FALSE),
  check.names = FALSE
)

overlap_table <- bind_rows(overlap_list)

write_csv(
  score_table,
  "results/gdsc2_breast_subtype_scores.csv"
)

write_csv(
  overlap_table,
  "results/gdsc2_signature_gene_overlap.csv"
)

saveRDS(
  score_table,
  "results/gdsc2_breast_subtype_scores.rds"
)

# 5. Match response data and audit repeated cell-line/drug pairs

response_matched <- response %>%
  filter(SANGER_MODEL_ID %in% model_ids) %>%
  inner_join(score_table, by = "SANGER_MODEL_ID")

duplicate_pairs <- response_matched %>%
  count(SANGER_MODEL_ID, DRUG_ID, name = "n_records") %>%
  filter(n_records > 1) %>%
  arrange(desc(n_records), DRUG_ID, SANGER_MODEL_ID)

write_csv(
  duplicate_pairs,
  "results/gdsc2_cell_line_drug_duplicate_audit.csv"
)

if (nrow(duplicate_pairs) > 0) {
  stop(
    nrow(duplicate_pairs),
    " cell-line/drug pairs have multiple response records. ",
    "Review results/gdsc2_cell_line_drug_duplicate_audit.csv and ",
    "Script 08 QC outputs before running association tests. ",
    "No duplicate records have been averaged or removed."
  )
}

# 6. Test signature-response associations

response_metrics <- c("LN_IC50", "AUC")
coverage_thresholds <- c(15, 20, 25)

association_results <- list()
result_index <- 1

for (metric in response_metrics) {
  
  for (subtype in subtype_names) {
    
    for (drug_id in unique(response_matched$DRUG_ID)) {
      
      drug_data <- response_matched %>%
        filter(DRUG_ID == drug_id) %>%
        filter(
          is.finite(.data[[metric]]),
          is.finite(.data[[subtype]])
        )
      
      n_models <- n_distinct(drug_data$SANGER_MODEL_ID)
      
      # Include >=15 so the primary >=20 and sensitivity >=15/25
      # rankings can all be derived from the same association results.
      if (n_models < 15) next
      
      if (
        sd(drug_data[[metric]]) == 0 ||
        sd(drug_data[[subtype]]) == 0
      ) {
        next
      }
      
      test <- suppressWarnings(
        cor.test(
          drug_data[[subtype]],
          drug_data[[metric]],
          method = "spearman",
          exact = FALSE
        )
      )
      
      association_results[[result_index]] <- data.frame(
        DRUG_ID = drug_id,
        DRUG_NAME = as.character(drug_data$DRUG_NAME[1]),
        subtype_comparison = subtype,
        response_metric = metric,
        n_models = n_models,
        rho = unname(test$estimate),
        p_value = test$p.value
      )
      
      result_index <- result_index + 1
    }
  }
}

if (length(association_results) == 0) {
  stop("No eligible drug/signature associations could be tested.")
}

associations <- bind_rows(association_results) %>%
  group_by(subtype_comparison, response_metric) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Negative rho means higher signature score is associated with
# lower LN_IC50 or AUC in this cell-line dataset.
# It is not evidence of clinical efficacy.

write_csv(
  associations,
  "results/gdsc2_subtype_drug_associations.csv"
)

# 7. Ranking robustness at each coverage threshold

robustness_results <- list()
robust_index <- 1

for (metric in response_metrics) {
  
  for (subtype in subtype_names) {
    
    for (threshold in coverage_thresholds) {
      
      ranked <- associations %>%
        filter(
          response_metric == metric,
          subtype_comparison == subtype,
          n_models >= threshold
        ) %>%
        arrange(rho, p_value, DRUG_ID) %>%
        mutate(
          rank = row_number(),
          coverage_threshold = threshold,
          coverage_role = ifelse(
            threshold == 20,
            "primary",
            "sensitivity"
          )
        )
      
      robustness_results[[robust_index]] <- ranked
      robust_index <- robust_index + 1
    }
  }
}

robustness_table <- bind_rows(robustness_results)

write_csv(
  robustness_table,
  "results/gdsc2_drug_ranking_robustness.csv"
)

# 8. Summary

cat("\nSubtype scores calculated for ",
    nrow(score_table), " cell lines.\n", sep = "")

cat("Unique cell-line/drug pairs: ",
    nrow(response_matched), "\n", sep = "")

cat("Association tests: ", nrow(associations), "\n", sep = "")

cat("\nDrugs represented at each coverage threshold:\n")

for (threshold in coverage_thresholds) {
  
  eligible <- associations %>%
    filter(n_models >= threshold)
  
  cat(
    "At least ", threshold, " cell lines: ",
    n_distinct(eligible$DRUG_ID), " distinct drugs\n",
    sep = ""
  )
}

cat("\nSignature gene overlap:\n")
print(overlap_table)

cat("\nTop 10 most negative associations for each subtype, ",
    "using LN_IC50 and at least 20 cell lines:\n", sep = "")

for (subtype in subtype_names) {
  
  cat("\n", subtype, "\n", sep = "")
  
  print(
    associations %>%
      filter(
        subtype_comparison == subtype,
        response_metric == "LN_IC50",
        n_models >= 20
      ) %>%
      arrange(rho, p_value) %>%
      select(
        DRUG_NAME, n_models, rho, p_value, FDR
      ) %>%
      head(10)
  )
}

cat("\nScript 09 finished.\n")

