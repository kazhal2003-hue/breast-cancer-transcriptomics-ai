# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 10 - evaluate ranking robustness
# Author: Kazhal Ghiasi
# Notes: Evaluates the stability of drug-prioritisation rankings across
# cell-line coverage thresholds, drug-response metrics, and
# breast cancer subtype signatures. Also assesses agreement
# between LN_IC50 and AUC within individual drugs.
# Primary coverage threshold: >= 20 unique cell lines per drug.
# Sensitivity thresholds: >= 15 and >= 25 unique cell lines.
# This script evaluates ranking stability, not clinical efficacy.

library(dplyr)
library(readr)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. Load and validate saved results

assoc_file <- "results/gdsc2_subtype_drug_associations.csv"
scores_file <- "results/gdsc2_breast_subtype_scores.rds"

matched_response_file <-
  "results/gdsc2_breast_response_records_matched_expression.rds"

full_response_file <-
  "results/gdsc2_breast_response_records.rds"

if (!file.exists(assoc_file)) {
  stop("Missing association results: ", assoc_file)
}

if (!file.exists(scores_file)) {
  stop("Missing subtype scores: ", scores_file)
}

if (file.exists(matched_response_file)) {
  response_file <- matched_response_file
} else if (file.exists(full_response_file)) {
  warning(
    "Matched response file not found. Using full breast response data ",
    "and matching to the scored cell lines."
  )
  response_file <- full_response_file
} else {
  stop("No GDSC2 response RDS file was found.")
}

assoc <- read_csv(assoc_file, show_col_types = FALSE)
response <- readRDS(response_file)
scores <- readRDS(scores_file)

required_assoc <- c(
  "DRUG_ID", "DRUG_NAME", "subtype_comparison",
  "response_metric", "n_models", "rho", "p_value", "FDR"
)

required_response <- c(
  "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME",
  "LN_IC50", "AUC"
)

required_scores <- c(
  "SANGER_MODEL_ID",
  "Basal_vs_rest", "Her2_vs_rest",
  "LumA_vs_rest", "LumB_vs_rest"
)

if (!all(required_assoc %in% names(assoc))) {
  stop("Association results are missing required columns.")
}

if (!all(required_response %in% names(response))) {
  stop("Response data are missing required columns.")
}

if (!all(required_scores %in% names(scores))) {
  stop("Subtype score data are missing required columns.")
}

thresholds <- c(15, 20, 25)

subtypes <- c(
  "Basal_vs_rest",
  "Her2_vs_rest",
  "LumA_vs_rest",
  "LumB_vs_rest"
)

metrics <- c("LN_IC50", "AUC")

response$SANGER_MODEL_ID <-
  trimws(as.character(response$SANGER_MODEL_ID))

scores$SANGER_MODEL_ID <-
  trimws(as.character(scores$SANGER_MODEL_ID))

assoc$DRUG_ID <- as.character(assoc$DRUG_ID)
response$DRUG_ID <- as.character(response$DRUG_ID)

if (anyDuplicated(
  assoc[c("DRUG_ID", "subtype_comparison", "response_metric")]
)) {
  stop(
    "Duplicate drug/signature/metric rows in association results. ",
    "Check Script 09 before calculating ranking stability."
  )
}

# Safe Jaccard similarity: returns NA when both sets are empty.
jaccard_similarity <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  
  union_size <- length(union(a, b))
  
  if (union_size == 0) {
    return(NA_real_)
  }
  
  length(intersect(a, b)) / union_size
}

# 2. Rank stability across coverage thresholds

ranking_list <- list()
i <- 1

for (metric in metrics) {
  for (subtype in subtypes) {
    for (threshold in thresholds) {
      
      ranked <- assoc %>%
        filter(
          response_metric == metric,
          subtype_comparison == subtype,
          n_models >= threshold,
          is.finite(rho)
        ) %>%
        arrange(rho, DRUG_ID) %>%
        mutate(rank = row_number())
      
      ranking_list[[i]] <- ranked %>%
        mutate(coverage_threshold = threshold)
      
      i <- i + 1
    }
  }
}

rankings <- bind_rows(ranking_list)

write_csv(
  rankings,
  "results/robustness_rankings_by_threshold.csv"
)

# Compare adjacent thresholds using drugs present at both thresholds.
# Ranks remain relative to each complete eligible drug list.

threshold_comparisons <- list()
i <- 1

for (metric in metrics) {
  for (subtype in subtypes) {
    
    for (j in seq_len(length(thresholds) - 1)) {
      
      low_threshold <- thresholds[j]
      high_threshold <- thresholds[j + 1]
      
      low <- rankings %>%
        filter(
          response_metric == metric,
          subtype_comparison == subtype,
          coverage_threshold == low_threshold
        ) %>%
        select(DRUG_ID, rho, rank)
      
      high <- rankings %>%
        filter(
          response_metric == metric,
          subtype_comparison == subtype,
          coverage_threshold == high_threshold
        ) %>%
        select(DRUG_ID, rho, rank)
      
      common <- inner_join(
        low, high,
        by = "DRUG_ID",
        suffix = c("_low", "_high")
      )
      
      rank_rho <- if (
        nrow(common) >= 3 &&
        sd(common$rank_low) > 0 &&
        sd(common$rank_high) > 0
      ) {
        cor(
          common$rank_low,
          common$rank_high,
          method = "spearman"
        )
      } else {
        NA_real_
      }
      
      top_low <- low %>%
        arrange(rank) %>%
        slice_head(n = 10) %>%
        pull(DRUG_ID)
      
      top_high <- high %>%
        arrange(rank) %>%
        slice_head(n = 10) %>%
        pull(DRUG_ID)
      
      threshold_comparisons[[i]] <- data.frame(
        subtype = subtype,
        metric = metric,
        threshold_low = low_threshold,
        threshold_high = high_threshold,
        drugs_low = nrow(low),
        drugs_high = nrow(high),
        common_drugs = nrow(common),
        rank_spearman = rank_rho,
        top10_overlap = length(intersect(top_low, top_high)),
        top10_jaccard = jaccard_similarity(top_low, top_high)
      )
      
      i <- i + 1
    }
  }
}

threshold_stability <- bind_rows(threshold_comparisons)

write_csv(
  threshold_stability,
  "results/robustness_threshold_stability.csv"
)

# 3. Agreement between LN_IC50 and AUC rankings

metric_comparisons <- list()
i <- 1

for (subtype in subtypes) {
  for (threshold in thresholds) {
    
    ln <- rankings %>%
      filter(
        subtype_comparison == subtype,
        response_metric == "LN_IC50",
        coverage_threshold == threshold
      ) %>%
      select(DRUG_ID, rho, rank)
    
    auc <- rankings %>%
      filter(
        subtype_comparison == subtype,
        response_metric == "AUC",
        coverage_threshold == threshold
      ) %>%
      select(DRUG_ID, rho, rank)
    
    common <- inner_join(
      ln, auc,
      by = "DRUG_ID",
      suffix = c("_ln", "_auc")
    )
    
    rank_rho <- if (
      nrow(common) >= 3 &&
      sd(common$rank_ln) > 0 &&
      sd(common$rank_auc) > 0
    ) {
      cor(
        common$rank_ln,
        common$rank_auc,
        method = "spearman"
      )
    } else {
      NA_real_
    }
    
    metric_comparisons[[i]] <- data.frame(
      subtype = subtype,
      coverage_threshold = threshold,
      drugs_ln_ic50 = nrow(ln),
      drugs_auc = nrow(auc),
      common_drugs = nrow(common),
      rank_spearman = rank_rho
    )
    
    i <- i + 1
  }
}

metric_stability <- bind_rows(metric_comparisons)

write_csv(
  metric_stability,
  "results/robustness_LNIC50_vs_AUC_rank_agreement.csv"
)

# 4. Agreement between subtype rankings

subtype_comparisons <- list()
i <- 1

subtype_pairs <- combn(subtypes, 2, simplify = FALSE)

for (metric in metrics) {
  for (threshold in thresholds) {
    for (pair in subtype_pairs) {
      
      a <- rankings %>%
        filter(
          subtype_comparison == pair[1],
          response_metric == metric,
          coverage_threshold == threshold
        ) %>%
        select(DRUG_ID, rho, rank)
      
      b <- rankings %>%
        filter(
          subtype_comparison == pair[2],
          response_metric == metric,
          coverage_threshold == threshold
        ) %>%
        select(DRUG_ID, rho, rank)
      
      common <- inner_join(
        a, b,
        by = "DRUG_ID",
        suffix = c("_a", "_b")
      )
      
      rank_rho <- if (
        nrow(common) >= 3 &&
        sd(common$rank_a) > 0 &&
        sd(common$rank_b) > 0
      ) {
        cor(
          common$rank_a,
          common$rank_b,
          method = "spearman"
        )
      } else {
        NA_real_
      }
      
      top_a <- a %>%
        arrange(rank) %>%
        slice_head(n = 10) %>%
        pull(DRUG_ID)
      
      top_b <- b %>%
        arrange(rank) %>%
        slice_head(n = 10) %>%
        pull(DRUG_ID)
      
      subtype_comparisons[[i]] <- data.frame(
        metric = metric,
        coverage_threshold = threshold,
        subtype_a = pair[1],
        subtype_b = pair[2],
        drugs_a = nrow(a),
        drugs_b = nrow(b),
        common_drugs = nrow(common),
        rank_spearman = rank_rho,
        top10_overlap = length(intersect(top_a, top_b)),
        top10_jaccard = jaccard_similarity(top_a, top_b)
      )
      
      i <- i + 1
    }
  }
}

subtype_stability <- bind_rows(subtype_comparisons)

write_csv(
  subtype_stability,
  "results/robustness_subtype_ranking_agreement.csv"
)

# 5. LN_IC50 and AUC agreement within each drug

# Restrict responses to cell lines with calculated subtype scores.
scored_ids <- unique(scores$SANGER_MODEL_ID)

response_matched <- response %>%
  filter(SANGER_MODEL_ID %in% scored_ids) %>%
  mutate(
    LN_IC50 = suppressWarnings(as.numeric(LN_IC50)),
    AUC = suppressWarnings(as.numeric(AUC))
  )

# The QC established uniqueness of cell-line/drug pairs.
# Recheck the matched data rather than silently collapsing records.
duplicate_pairs <- response_matched %>%
  count(SANGER_MODEL_ID, DRUG_ID, name = "n_records") %>%
  filter(n_records > 1)

if (nrow(duplicate_pairs) > 0) {
  write_csv(
    duplicate_pairs,
    "results/robustness_duplicate_response_pairs.csv"
  )
  
  stop(
    "Repeated cell-line/drug response pairs found in the matched data. ",
    "See results/robustness_duplicate_response_pairs.csv."
  )
}

within_drug_agreement <- response_matched %>%
  group_by(DRUG_ID, DRUG_NAME) %>%
  summarise(
    n_models = n_distinct(
      SANGER_MODEL_ID[
        is.finite(LN_IC50) & is.finite(AUC)
      ]
    ),
    metric_spearman = {
      valid <- is.finite(LN_IC50) & is.finite(AUC)
      
      x <- LN_IC50[valid]
      y <- AUC[valid]
      
      if (
        length(x) >= 5 &&
        sd(x) > 0 &&
        sd(y) > 0
      ) {
        cor(x, y, method = "spearman")
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

write_csv(
  within_drug_agreement,
  "results/robustness_within_drug_metric_agreement.csv"
)

# 6. Print summaries

print(as.data.frame(threshold_stability), row.names = FALSE)
print(as.data.frame(metric_stability), row.names = FALSE)
print(as.data.frame(subtype_stability), row.names = FALSE)

cat("\n--- Within-drug LN_IC50 vs AUC agreement ---\n")

valid_within_drug <- within_drug_agreement %>%
  filter(is.finite(metric_spearman))

if (nrow(valid_within_drug) > 0) {
  
  print(
    valid_within_drug %>%
      summarise(
        drugs_assessed = n(),
        median_spearman = median(metric_spearman),
        minimum_spearman = min(metric_spearman),
        maximum_spearman = max(metric_spearman)
      )
  )
  
} else {
  cat("No drugs had sufficient non-constant paired response data.\n")
}

cat("\nScript 10 finished.\n")

