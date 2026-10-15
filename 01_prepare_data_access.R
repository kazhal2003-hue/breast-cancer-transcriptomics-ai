# Project: From Cancer Transcriptomics to AI-Guided Drug Prioritization
# Script: 01_prepare_data_access.R
# Author: Kazhal Ghiasi
# Purpose: Prepare the R environment and dependencies for TCGA-BRCA analysis.
# Inputs:
#   - An R installation
#   - Internet access if required packages need to be installed
# Outputs:
#   - data/ and results/ directories, if they do not already exist
#   - R environment and package-version information in the console
#   - results/session_info.txt
# Note:
#   This script prepares the environment; it does not download TCGA-BRCA data.

# 1. Report the R environment 

cat("\n=== R ENVIRONMENT ===\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", normalizePath(getwd()), "\n")

# 2. Install BiocManager if necessary 

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop(
    "BiocManager could not be installed or loaded. ",
    "Check your internet connection and R installation."
  )
}

# 3. Install TCGAbiolinks if necessary 

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install(
    "TCGAbiolinks",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop(
    "TCGAbiolinks is unavailable. ",
    "Check the installation messages before continuing."
  )
}

# 4. Create the project directories 

project_directories <- c("data", "results")

for (directory in project_directories) {
  if (!dir.exists(directory)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  if (!dir.exists(directory)) {
    stop("Could not create or access directory: ", directory)
  }
}

cat("\n=== PROJECT DIRECTORIES ===\n")
print(data.frame(
  directory = project_directories,
  exists = dir.exists(project_directories)
))

# 5. Report package versions 

cat("\n=== PACKAGE VERSIONS ===\n")

package_versions <- c(
  BiocManager = as.character(packageVersion("BiocManager")),
  TCGAbiolinks = as.character(packageVersion("TCGAbiolinks"))
)

print(package_versions)

# 6. Record session information 

session_file <- file.path("results", "session_info.txt")

writeLines(
  capture.output(sessionInfo()),
  con = session_file
)

cat("\nSession information saved to:", session_file, "\n")

# 7. Confirm successful completion 

cat("\n=== ENVIRONMENT PREPARATION COMPLETE ===\n")
cat("TCGAbiolinks is available:", 
    requireNamespace("TCGAbiolinks", quietly = TRUE), "\n")
cat("Data directory:", normalizePath("data"), "\n")
cat("Results directory:", normalizePath("results"), "\n")
cat("Note: TCGA-BRCA data have not been downloaded by this script.\n")