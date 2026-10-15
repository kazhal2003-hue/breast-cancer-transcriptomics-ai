# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 04 - ml subtype benchmark cv safe
# Author: Kazhal Ghiasi
# Purpose: Re-evaluate PAM50 subtype classification with fold-specific
# feature selection to prevent cross-validation feature-selection leakage.
# Method: Multinomial LASSO, stratified 5-fold CV, held-out test evaluation.
# Note: All existing project inputs and benchmark outputs are preserved.

set.seed(42)

suppressPackageStartupMessages({
  library(glmnet)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. Load existing inputs without modifying them

vst <- readRDS("results/ml_tumour_vst_matrix_1095.rds")
meta <- read.csv(
  "results/ml_tumour_metadata_1095.csv",
  stringsAsFactors = FALSE
)
split_file <- "results/ml_subtype_train_test_split.csv"

stopifnot(
  is.matrix(vst),
  all(c("expression_id", "subtype") %in% names(meta)),
  !is.null(colnames(vst))
)

# Align metadata with the existing expression-matrix columns.
meta <- meta[match(colnames(vst), meta$expression_id), ]

if (anyNA(meta$expression_id) || anyNA(meta$subtype)) {
  stop("Could not align expression samples and subtype metadata.")
}

# Samples are rows and genes are columns for glmnet.
x <- t(vst)
storage.mode(x) <- "double"

y <- factor(meta$subtype)
if (anyNA(y)) stop("Subtype labels contain missing values.")

cat("Samples:", nrow(x), "\n")
cat("Genes:", ncol(x), "\n")
print(table(y))

# 2. Recover the original train/test split; do not generate a new one

if (!file.exists(split_file)) {
  stop("Original split file not found. No outputs were changed.")
}

sp <- read.csv(split_file, stringsAsFactors = FALSE)

split_col <- intersect(
  c("split", "set", "partition", "dataset", "subset"),
  names(sp)
)[1]

id_candidates <- intersect(
  c("expression_id", "sample_id", "full_barcode"),
  names(sp)
)

if (is.na(split_col) || length(id_candidates) == 0) {
  stop(
    "Could not identify ID and train/test columns in the saved split file. ",
    "No outputs were changed."
  )
}

# Find the ID column with the greatest overlap with the metadata.
meta_id_candidates <- intersect(
  c("expression_id", "sample_id", "full_barcode"),
  names(meta)
)

scores <- expand.grid(
  split_id = id_candidates,
  meta_id = meta_id_candidates,
  stringsAsFactors = FALSE
)

scores$overlap <- mapply(function(a, b) {
  length(intersect(as.character(sp[[a]]), as.character(meta[[b]])))
}, scores$split_id, scores$meta_id)

best <- scores[which.max(scores$overlap), ]

if (best$overlap < nrow(meta)) {
  stop(
    "The saved split does not match all expression samples. ",
    "No outputs were changed."
  )
}

split_map <- setNames(
  tolower(trimws(as.character(sp[[split_col]]))),
  as.character(sp[[best$split_id]])
)

split_label <- unname(split_map[as.character(meta[[best$meta_id]])])

train_idx <- which(grepl("^train", split_label))
test_idx  <- which(grepl("^test", split_label))

if (length(train_idx) != 877 || length(test_idx) != 218 ||
    length(intersect(train_idx, test_idx)) != 0) {
  stop(
    "The saved split did not resolve to the expected 877 training ",
    "and 218 test samples. No outputs were changed."
  )
}

x_train <- x[train_idx, , drop = FALSE]
x_test  <- x[test_idx, , drop = FALSE]
y_train <- droplevels(y[train_idx])
y_test  <- factor(y[test_idx], levels = levels(y_train))

cat("\nTraining samples:", length(train_idx), "\n")
cat("Held-out test samples:", length(test_idx), "\n")

# 3. Create stratified folds within the training cohort

set.seed(42)
k <- 5
foldid <- integer(length(y_train))

for (lev in levels(y_train)) {
  idx <- which(y_train == lev)
  idx <- sample(idx)
  foldid[idx] <- rep(seq_len(k), length.out = length(idx))
}

# Fixed, decreasing lambda grid; no held-out test data are used for tuning.
lambda_grid <- exp(seq(log(1), log(1e-5), length.out = 100))
n_select <- min(1000L, ncol(x_train))

# 4. Cross-validation with feature selection performed inside each fold

cv_error <- matrix(NA_real_, nrow = k, ncol = length(lambda_grid))

for (fold in seq_len(k)) {
  cat("Cross-validation fold", fold, "of", k, "\n")
  
  analysis_idx <- which(foldid != fold)
  validation_idx <- which(foldid == fold)
  
  xa <- x_train[analysis_idx, , drop = FALSE]
  xv <- x_train[validation_idx, , drop = FALSE]
  ya <- droplevels(y_train[analysis_idx])
  yv <- y_train[validation_idx]
  
  # Crucial correction: choose variable genes using this fold's
  # analysis samples only, never its validation samples.
  gene_var <- apply(xa, 2, var)
  gene_var[!is.finite(gene_var)] <- -Inf
  selected <- order(gene_var, decreasing = TRUE)[
    seq_len(min(n_select, sum(is.finite(gene_var))))
  ]
  
  if (length(selected) == 0) stop("No usable genes in a CV fold.")
  
  fit <- glmnet(
    x = xa[, selected, drop = FALSE],
    y = ya,
    family = "multinomial",
    alpha = 1,
    lambda = lambda_grid
  )
  
  for (j in seq_along(lambda_grid)) {
    pred <- predict(
      fit,
      newx = xv[, selected, drop = FALSE],
      s = lambda_grid[j],
      type = "class"
    )
    pred <- factor(as.character(pred), levels = levels(y_train))
    cv_error[fold, j] <- mean(pred != yv)
  }
}

mean_error <- colMeans(cv_error)
se_error <- apply(cv_error, 2, sd) / sqrt(k)

best_idx <- which.min(mean_error)
eligible <- which(mean_error <= mean_error[best_idx] + se_error[best_idx])

# Largest eligible lambda: the one-standard-error rule.
chosen_idx <- min(eligible)
chosen_lambda <- lambda_grid[chosen_idx]

cat("\nSelected lambda:", chosen_lambda, "\n")
cat("Mean CV classification error:", mean_error[chosen_idx], "\n")

# 5. Select genes using ALL training samples, then fit final model

gene_var <- apply(x_train, 2, var)
gene_var[!is.finite(gene_var)] <- -Inf
selected_genes <- order(gene_var, decreasing = TRUE)[
  seq_len(min(n_select, sum(is.finite(gene_var))))
]

final_fit <- glmnet(
  x = x_train[, selected_genes, drop = FALSE],
  y = y_train,
  family = "multinomial",
  alpha = 1,
  lambda = lambda_grid
)

pred <- predict(
  final_fit,
  newx = x_test[, selected_genes, drop = FALSE],
  s = chosen_lambda,
  type = "class"
)

predicted <- factor(as.character(pred), levels = levels(y_train))
actual <- factor(as.character(y_test), levels = levels(y_train))

# 6. Evaluate held-out test samples

confusion <- table(Actual = actual, Predicted = predicted)

class_metrics <- do.call(rbind, lapply(levels(y_train), function(cl) {
  tp <- sum(actual == cl & predicted == cl)
  fp <- sum(actual != cl & predicted == cl)
  fn <- sum(actual == cl & predicted != cl)
  
  precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  recall <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) &&
            precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  
  data.frame(
    subtype = cl,
    precision = precision,
    recall = recall,
    F1 = f1,
    support = sum(actual == cl)
  )
}))

accuracy <- mean(actual == predicted)
balanced_accuracy <- mean(class_metrics$recall, na.rm = TRUE)
macro_f1 <- mean(class_metrics$F1, na.rm = TRUE)

metrics <- data.frame(
  model = "Multinomial LASSO; fold-specific feature selection",
  samples_train = length(train_idx),
  samples_test = length(test_idx),
  genes_selected = length(selected_genes),
  CV_folds = k,
  selected_lambda = chosen_lambda,
  test_accuracy = accuracy,
  test_balanced_accuracy = balanced_accuracy,
  test_macro_F1 = macro_f1
)

cat("\n--- Corrected held-out test results ---\n")
print(metrics)
cat("\nConfusion matrix:\n")
print(confusion)
cat("\nPer-subtype metrics:\n")
print(class_metrics)

# 7. Save to NEW filenames; preserve every original output

write.csv(
  metrics,
  "results/ml_subtype_cv_safe_benchmark_metrics.csv",
  row.names = FALSE
)
write.csv(
  class_metrics,
  "results/ml_subtype_cv_safe_class_metrics.csv",
  row.names = FALSE
)
write.csv(
  as.data.frame.matrix(confusion),
  "results/ml_subtype_cv_safe_confusion_matrix.csv"
)
write.csv(
  data.frame(
    expression_id = meta$expression_id[test_idx],
    actual_subtype = actual,
    predicted_subtype = predicted
  ),
  "results/ml_subtype_cv_safe_predictions.csv",
  row.names = FALSE
)
write.csv(
  data.frame(gene = colnames(x_train)[selected_genes]),
  "results/ml_subtype_cv_safe_selected_genes.csv",
  row.names = FALSE
)
saveRDS(
  list(
    model = final_fit,
    selected_genes = colnames(x_train)[selected_genes],
    lambda_grid = lambda_grid,
    chosen_lambda = chosen_lambda,
    cv_error = cv_error,
    foldid = foldid,
    train_indices = train_idx,
    test_indices = test_idx,
    metrics = metrics,
    session_info = sessionInfo()
  ),
  "results/ml_subtype_cv_safe_model_details.rds"
)

cat("\nFinished. Existing files were not overwritten.\n")
