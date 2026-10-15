# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 13 - effect size robustness
# Author: Kazhal Ghiasi
# Notes: Compares drug-response association effect estimates from the
# full matched breast cancer cell-line panel and the Broad-only
# sensitivity analysis.
# Primary response metric: LN_IC50
# Primary drug-coverage threshold: >= 20 unique cell lines
# Outputs:
#   - broad_only_effect_size_robustness.csv
#   - broad_only_largest_effect_changes.csv
#   - broad_only_effect_size_robustness.png
#   - broad_only_effect_size_robustness_details.rds
# This is a sensitivity analysis of exploratory cell-line
# associations, not evidence of clinical efficacy.

# 0. Packages and input validation 

required_packages <- c("dplyr", "readr", "ggplot2", "tibble")

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
library(ggplot2)
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

# 1. Load saved association results 
# Previous analyses are not rerun in this script.

original <- read_csv(
  original_file,
  show_col_types = FALSE
)

broad <- read_csv(
  broad_file,
  show_col_types = FALSE
)

# 2. Standardise and validate input columns 

standardise_associations <- function(data, dataset_name) {
  
  if (!"metric" %in% names(data)) {
    if ("response_metric" %in% names(data)) {
      data <- data %>%
        rename(metric = response_metric)
    } else {
      stop(
        dataset_name,
        " results must contain 'metric' or 'response_metric'."
      )
    }
  }
  
  required_columns <- c(
    "DRUG_ID",
    "subtype_comparison",
    "metric",
    "n_models",
    "rho"
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
  
  if (!"DRUG_NAME" %in% names(data)) {
    data$DRUG_NAME <- NA_character_
  }
  
  data <- data %>%
    mutate(
      DRUG_ID = as.character(DRUG_ID),
      DRUG_NAME = as.character(DRUG_NAME),
      subtype_comparison = as.character(subtype_comparison),
      metric = as.character(metric),
      n_models = as.numeric(n_models),
      rho = as.numeric(rho)
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

# 3. Filter to the primary analysis 

primary_metric <- "LN_IC50"
primary_threshold <- 20

original <- original %>%
  filter(
    metric == primary_metric,
    n_models >= primary_threshold,
    is.finite(rho)
  )

broad <- broad %>%
  filter(
    metric == primary_metric,
    n_models >= primary_threshold,
    is.finite(rho)
  )

if (nrow(original) == 0) {
  stop(
    "No eligible full-panel LN_IC50 associations remain ",
    "at the >= 20-cell-line threshold."
  )
}

if (nrow(broad) == 0) {
  stop(
    "No eligible Broad-only LN_IC50 associations remain ",
    "at the >= 20-cell-line threshold."
  )
}

# 4. Check uniqueness of drug-association records 

check_duplicate_associations <- function(data, dataset_name) {
  
  duplicates <- data %>%
    count(
      subtype_comparison,
      DRUG_ID,
      name = "records"
    ) %>%
    filter(records > 1)
  
  if (nrow(duplicates) > 0) {
    audit_file <- file.path(
      results_dir,
      paste0(
        "script13_",
        dataset_name,
        "_duplicate_association_audit.csv"
      )
    )
    
    write_csv(duplicates, audit_file)
    
    stop(
      dataset_name,
      " results contain duplicate drug/subtype records. ",
      "Inspect: ",
      audit_file
    )
  }
}

check_duplicate_associations(
  original,
  "full_panel"
)

check_duplicate_associations(
  broad,
  "broad_only"
)

# 5. Match drugs within each shared subtype signature 

signatures <- intersect(
  unique(original$subtype_comparison),
  unique(broad$subtype_comparison)
)

if (length(signatures) == 0) {
  stop(
    "No shared subtype comparisons were found ",
    "between the full-panel and Broad-only analyses."
  )
}

comparison_list <- lapply(signatures, function(sig) {
  
  original_sig <- original %>%
    filter(subtype_comparison == sig) %>%
    select(
      DRUG_ID,
      DRUG_NAME,
      rho_original = rho
    )
  
  broad_sig <- broad %>%
    filter(subtype_comparison == sig) %>%
    select(
      DRUG_ID,
      rho_broad = rho
    )
  
  shared <- inner_join(
    original_sig,
    broad_sig,
    by = "DRUG_ID"
  ) %>%
    filter(
      is.finite(rho_original),
      is.finite(rho_broad)
    ) %>%
    mutate(
      signature = sig,
      effect_change = rho_broad - rho_original,
      absolute_change = abs(effect_change)
    )
  
  shared
})

effect_comparison <- bind_rows(comparison_list)

if (nrow(effect_comparison) == 0) {
  stop(
    "No shared drugs with finite effect estimates remain ",
    "after matching the two analyses."
  )
}

write_csv(
  effect_comparison,
  file.path(
    results_dir,
    "broad_only_effect_size_comparison_all_shared_drugs.csv"
  )
)

# 6. Safe correlation helper 

safe_correlation <- function(x, y, method = "spearman") {
  
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  
  if (
    length(x) < 3 ||
    dplyr::n_distinct(x) < 2 ||
    dplyr::n_distinct(y) < 2
  ) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(x, y, method = method)
  )
}

# 7. Summarise effect-estimate agreement 

effect_summary <- effect_comparison %>%
  group_by(signature) %>%
  summarise(
    shared_drugs = n(),
    
    spearman_correlation = safe_correlation(
      rho_original,
      rho_broad,
      method = "spearman"
    ),
    
    pearson_correlation = safe_correlation(
      rho_original,
      rho_broad,
      method = "pearson"
    ),
    
    median_absolute_change = median(
      absolute_change,
      na.rm = TRUE
    ),
    
    mean_absolute_change = mean(
      absolute_change,
      na.rm = TRUE
    ),
    
    proportion_change_ge_0_10_percent =
      mean(absolute_change >= 0.10) * 100,
    
    direction_comparable_drugs = sum(
      rho_original != 0 & rho_broad != 0
    ),
    
    proportion_same_direction_percent = {
      nonzero <- rho_original != 0 & rho_broad != 0
      
      if (sum(nonzero) > 0) {
        mean(
          sign(rho_original[nonzero]) ==
            sign(rho_broad[nonzero])
        ) * 100
      } else {
        NA_real_
      }
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

write_csv(
  effect_summary,
  file.path(
    results_dir,
    "broad_only_effect_size_robustness.csv"
  )
)

# 8. Identify the largest effect-estimate changes 

largest_changes <- effect_comparison %>%
  group_by(signature) %>%
  arrange(
    desc(absolute_change),
    DRUG_ID,
    .by_group = TRUE
  ) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  arrange(
    signature,
    desc(absolute_change)
  )

write_csv(
  largest_changes,
  file.path(
    results_dir,
    "broad_only_largest_effect_changes.csv"
  )
)

# 9. Plot full-panel versus Broad-only estimates 
# Each point represents one drug within one subtype signature.
# The dashed diagonal represents identical correlation estimates.
# The pink ribbon is the 95% confidence band around the fitted
# linear trend, not an individual-drug confidence interval.

p <- ggplot(
  effect_comparison,
  aes(
    x = rho_original,
    y = rho_broad
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#C2185B",
    fill = "#F8BBD0",
    alpha = 0.65,
    linewidth = 0.9
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.7
  ) +
  geom_point(
    color = "#AD1457",
    alpha = 0.55,
    size = 1.5
  ) +
  facet_wrap(
    ~ signature,
    ncol = 2
  ) +
  coord_equal(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  ) +
  labs(
    title = "Robustness of Drug–Response Associations",
    subtitle = paste(
      "Full matched panel versus Broad-only sensitivity analysis",
      "| GDSC2 breast cancer models"
    ),
    x = "Spearman rho in full-panel analysis",
    y = "Spearman rho in Broad-only analysis",
    caption = paste(
      "Each point represents a shared drug within a subtype signature.",
      "Dashed line: identical effect estimates.",
      "Pink band: 95% confidence band around the fitted linear trend."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      color = "#880E4F",
      size = 15
    ),
    plot.subtitle = element_text(
      color = "grey35"
    ),
    strip.text = element_text(
      face = "bold",
      color = "#880E4F"
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "#F8BBD0",
      fill = NA,
      linewidth = 0.7
    ),
    plot.caption = element_text(
      color = "grey40",
      size = 9
    )
  )

plot_file <- file.path(
  results_dir,
  "broad_only_effect_size_robustness.png"
)

ggsave(
  filename = plot_file,
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# 10. Save analysis settings and summary 

details <- list(
  script = "13_effect_size_robustness.R",
  author = "Kazhal Ghiasi",
  project = "From Cancer Transcriptomics to AI-Guided Drug Prioritization",
  full_panel_input = original_file,
  broad_only_input = broad_file,
  primary_metric = primary_metric,
  primary_minimum_cell_lines = primary_threshold,
  shared_signatures = signatures,
  number_of_shared_drug_signature_pairs = nrow(effect_comparison),
  effect_change_definition = "rho_broad - rho_original",
  absolute_change_definition = "abs(rho_broad - rho_original)",
  direction_analysis = "Exact zero correlations excluded",
  plot_confidence_band = "95% confidence band around fitted linear trend",
  effect_summary = effect_summary,
  largest_effect_changes = largest_changes,
  all_shared_drug_comparisons = effect_comparison
)

saveRDS(
  details,
  file.path(
    results_dir,
    "broad_only_effect_size_robustness_details.rds"
  )
)

# 11. Report results 

cat("\n=== EFFECT-SIZE ROBUSTNESS SUMMARY ===\n")
cat("Primary metric:", primary_metric, "\n")
cat("Minimum cell lines per drug:", primary_threshold, "\n")
cat(
  "Shared drug/subtype comparisons:",
  nrow(effect_comparison),
  "\n\n"
)

print(
  effect_summary,
  n = Inf,
  width = Inf
)

cat("\nSaved files:\n")
cat(
  "- results/broad_only_effect_size_robustness.csv\n",
  "- results/broad_only_largest_effect_changes.csv\n",
  "- results/broad_only_effect_size_robustness.png\n",
  "- results/broad_only_effect_size_robustness_details.rds\n",
  "- results/broad_only_effect_size_comparison_all_shared_drugs.csv\n"
)

print(p)

cat("\n=== SCRIPT 13 COMPLETE ===\n")

