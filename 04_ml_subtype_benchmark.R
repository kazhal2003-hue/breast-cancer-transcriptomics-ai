# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script 04 — Machine-learning subtype benchmark
# Author: Kazhal Ghiasi
# Purpose:
# Benchmark multinomial LASSO classification of breast cancer subtypes
# using TCGA-BRCA tumour expression data.
# Design:
# - Use the saved VST expression matrix and matching subtype metadata.
# - Exclude missing or unknown subtype labels.
# - Create a reproducible, stratified training/test split.
# - Select variable genes using training data only.
# - Tune the multinomial LASSO using stratified cross-validation.
# - Evaluate predictions on the held-out test set.
# Important:
# This is a subtype-classification benchmark, not a drug-response predictor.

# 1. Check dependencies 

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop(
    "Package 'glmnet' is required. Install it with install.packages('glmnet'), ",
    "then rerun this script."
  )
}

library(glmnet)

# 2. Define inputs and outputs 

expr_file <- "results/ml_tumour_vst_matrix_1095.rds"
meta_file <- "results/ml_tumour_metadata_1095.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

if (!file.exists(expr_file)) {
  stop("Expression matrix not found: ", expr_file)
}

if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}

# 3. Load data 

expr <- readRDS(expr_file)
meta <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("expression_id", "sample_id", "subtype")

missing_columns <- setdiff(required_columns, colnames(meta))

if (length(missing_columns) > 0) {
  stop(
    "Metadata is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!is.matrix(expr) || !is.numeric(expr)) {
  stop("The expression object must be a numeric matrix.")
}

if (is.null(colnames(expr)) || is.null(rownames(expr))) {
  stop("The expression matrix must have gene row names and sample column names.")
}

if (anyDuplicated(colnames(expr))) {
  stop("Duplicate expression sample IDs were found.")
}

if (anyDuplicated(meta$expression_id)) {
  stop("Duplicate expression_id values were found in the metadata.")
}

if (ncol(expr) != nrow(meta)) {
  stop(
    "Expression/metadata dimensions do not match: ",
    ncol(expr), " expression samples versus ", nrow(meta), " metadata rows."
  )
}

# Match metadata to expression columns explicitly.
match_idx <- match(colnames(expr), meta$expression_id)

if (anyNA(match_idx)) {
  stop("Some expression columns have no matching metadata record.")
}

meta <- meta[match_idx, , drop = FALSE]

if (!identical(colnames(expr), meta$expression_id)) {
  stop("Expression and metadata sample order could not be aligned.")
}

if (any(!is.finite(expr))) {
  stop("The expression matrix contains missing or non-finite values.")
}

cat("\n=== INPUT DATA ===\n")
cat("Genes:", nrow(expr), "\n")
cat("Expression samples:", ncol(expr), "\n")
cat("Metadata rows:", nrow(meta), "\n")

# 4. Retain samples with known subtype labels 

subtype_text <- trimws(as.character(meta$subtype))

unknown_subtype <- is.na(subtype_text) |
  toupper(subtype_text) %in% c("", "NA", "N/A", "UNKNOWN", "NOT AVAILABLE")

keep <- !unknown_subtype

expr <- expr[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

meta$subtype <- droplevels(factor(trimws(as.character(meta$subtype))))

if (nlevels(meta$subtype) < 2) {
  stop("At least two known subtype classes are required.")
}

if (any(table(meta$subtype) < 2)) {
  stop(
    "At least one subtype has fewer than two samples; ",
    "a stratified train/test split cannot be performed safely."
  )
}

cat("\n=== SUBTYPE COUNTS ===\n")
print(table(meta$subtype))

# 5. Create a reproducible stratified train/test split 

set.seed(42)

test_idx <- unlist(
  lapply(
    split(seq_len(nrow(meta)), meta$subtype),
    function(idx) {
      n_test <- min(
        length(idx) - 1L,
        max(1L, round(0.20 * length(idx)))
      )
      
      sample(idx, size = n_test)
    }
  ),
  use.names = FALSE
)

test_idx <- sort(test_idx)
train_idx <- setdiff(seq_len(nrow(meta)), test_idx)

# glmnet expects samples in rows and genes in columns.
x <- t(expr)
y <- meta$subtype

y_train <- droplevels(y[train_idx])
y_test <- y[test_idx]

if (any(table(y_train) < 5)) {
  stop(
    "At least one subtype has fewer than five training samples. ",
    "Five-fold cross-validation would not be appropriate."
  )
}

cat("\n=== TRAIN/TEST SPLIT ===\n")
cat("Training samples:", length(train_idx), "\n")
cat("Test samples:", length(test_idx), "\n")

cat("\nTraining subtype counts:\n")
print(table(y_train))

cat("\nTest subtype counts:\n")
print(table(y_test))

# 6. Select variable genes using training data only 

gene_var <- apply(
  x[train_idx, , drop = FALSE],
  2,
  var
)

valid_genes <- is.finite(gene_var)

gene_var <- gene_var[valid_genes]

if (length(gene_var) == 0) {
  stop("No genes with valid training-set variance were found.")
}

n_top <- min(1000L, length(gene_var))

top_genes <- names(
  sort(gene_var, decreasing = TRUE)
)[seq_len(n_top)]

x_train <- x[train_idx, top_genes, drop = FALSE]
x_test <- x[test_idx, top_genes, drop = FALSE]

cat("\nGenes selected using training data:", length(top_genes), "\n")

# 7. Create stratified cross-validation folds 

set.seed(42)

nfolds <- 5L

foldid <- integer(length(train_idx))

for (class_name in levels(y_train)) {
  class_positions <- which(y_train == class_name)
  
  foldid[class_positions] <- sample(
    rep(seq_len(nfolds), length.out = length(class_positions))
  )
}

stopifnot(
  length(foldid) == length(train_idx),
  all(foldid %in% seq_len(nfolds))
)

# 8. Train multinomial LASSO with cross-validation 

set.seed(42)

cv_fit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "multinomial",
  alpha = 1,
  type.measure = "class",
  foldid = foldid
)

# 9. Predict the held-out test samples 

pred_raw <- predict(
  cv_fit,
  newx = x_test,
  s = "lambda.1se",
  type = "class"
)

pred <- factor(
  as.character(pred_raw[, 1]),
  levels = levels(y)
)

actual <- factor(
  as.character(y_test),
  levels = levels(y)
)

# 10. Calculate evaluation metrics 

confusion <- table(
  Actual = actual,
  Predicted = pred
)

class_support <- rowSums(confusion)
true_positive <- diag(confusion)

class_recall <- true_positive / class_support

predicted_support <- colSums(confusion)

class_precision <- ifelse(
  predicted_support > 0,
  true_positive / predicted_support,
  0
)

class_f1 <- ifelse(
  class_precision + class_recall > 0,
  2 * class_precision * class_recall /
    (class_precision + class_recall),
  0
)

class_metrics <- data.frame(
  subtype = levels(actual),
  support = as.integer(class_support),
  precision = as.numeric(class_precision),
  recall = as.numeric(class_recall),
  f1_score = as.numeric(class_f1)
)

overall_accuracy <- mean(pred == actual)
balanced_accuracy <- mean(class_recall)
macro_f1 <- mean(class_f1)

metrics <- data.frame(
  metric = c(
    "test_accuracy",
    "balanced_accuracy",
    "macro_f1",
    "selected_lambda_1se",
    "training_samples",
    "test_samples",
    "genes_selected",
    "cross_validation_folds"
  ),
  value = c(
    overall_accuracy,
    balanced_accuracy,
    macro_f1,
    cv_fit$lambda.1se,
    length(train_idx),
    length(test_idx),
    length(top_genes),
    nfolds
  )
)

cat("\n=== MODEL PERFORMANCE ===\n")
cat("Selected lambda (lambda.1se):", cv_fit$lambda.1se, "\n")
cat("Held-out test accuracy:", round(overall_accuracy, 4), "\n")
cat("Balanced accuracy:", round(balanced_accuracy, 4), "\n")
cat("Macro F1 score:", round(macro_f1, 4), "\n")

cat("\n=== CONFUSION MATRIX ===\n")
print(confusion)

cat("\n=== PER-SUBTYPE METRICS ===\n")
print(class_metrics, row.names = FALSE)

# 11. Save predictions, metrics, and model 

ml_results <- data.frame(
  sample_id = meta$sample_id[test_idx],
  expression_id = meta$expression_id[test_idx],
  actual_subtype = as.character(actual),
  predicted_subtype = as.character(pred),
  correct = pred == actual,
  stringsAsFactors = FALSE
)

split_results <- data.frame(
  sample_id = meta$sample_id,
  expression_id = meta$expression_id,
  subtype = as.character(y),
  partition = ifelse(
    seq_len(nrow(meta)) %in% test_idx,
    "test",
    "train"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  ml_results,
  "results/ml_subtype_predictions.csv",
  row.names = FALSE
)

write.csv(
  metrics,
  "results/ml_subtype_benchmark_metrics.csv",
  row.names = FALSE
)

write.csv(
  class_metrics,
  "results/ml_subtype_class_metrics.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(confusion),
  "results/ml_subtype_confusion_matrix.csv"
)

write.csv(
  split_results,
  "results/ml_subtype_train_test_split.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = top_genes),
  "results/ml_subtype_selected_genes.csv",
  row.names = FALSE
)

saveRDS(
  cv_fit,
  "results/ml_subtype_lasso_cv.rds"
)

saveRDS(
  list(
    seed = 42,
    nfolds = nfolds,
    selected_genes = top_genes,
    train_indices = train_idx,
    test_indices = test_idx,
    lambda_1se = cv_fit$lambda.1se,
    overall_accuracy = overall_accuracy,
    balanced_accuracy = balanced_accuracy,
    macro_f1 = macro_f1,
    confusion_matrix = confusion,
    class_metrics = class_metrics
  ),
  "results/ml_subtype_benchmark_details.rds"
)

cat("\n=== SCRIPT 04 COMPLETE ===\n")
cat("Predictions, metrics, split details, selected genes, and model saved.\n")
