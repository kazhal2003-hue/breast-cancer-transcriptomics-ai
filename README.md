# breast-cancer-transcriptomics-ai
A reproducible breast cancer transcriptomics and AI-guided drug prioritization project using TCGA-BRCA and GDSC2 data.
# From Cancer Transcriptomics to AI-Guided Drug Prioritization

An integrative breast cancer transcriptomics and machine-learning project investigating molecular subtypes, gene-expression signatures, and the robustness of computational drug prioritisation using public datasets.

## Overview

Breast cancer comprises molecularly distinct subtypes with different transcriptional profiles and biological characteristics. Understanding these differences can help identify subtype-associated genes, evaluate the predictive value of expression signatures, and explore potential therapeutic hypotheses.

This project combines differential gene-expression analysis, machine-learning classification, subtype-signature construction, and drug-response association analysis in a reproducible R-based workflow.

Rather than treating drug ranking as an isolated task, the project investigates how molecular subtype information and alternative expression signatures influence computational drug prioritisation, and whether the resulting rankings remain consistent under different analytical assumptions.

## Research objectives

1. Identify genes differentially expressed between matched breast tumour and normal tissue samples.
2. Evaluate machine-learning classification of breast cancer molecular subtypes using gene-expression data.
3. Construct subtype-associated expression signatures and assess their robustness.
4. Investigate associations between subtype signatures and drug-response measurements in cancer cell lines.
5. Evaluate the stability of drug rankings across response metrics, data-source restrictions, and sensitivity analyses.

## Data sources

| Dataset | Purpose |
|---|---|
| TCGA-BRCA | Breast cancer tumour and normal tissue gene-expression analysis |
| PAM50 molecular subtype annotations | Subtype labels for expression analysis and machine-learning classification |
| GDSC2 | Cancer cell-line drug-response measurements for exploratory drug prioritisation |
| Broad-source cell-line subset | Sensitivity analysis restricted to a defined cell-line data source |

The analyses use publicly available data and derived project-specific files. Original datasets are obtained from their respective repositories, with source information and processing details documented alongside the relevant scripts.

## Analytical workflow

### 1. Differential gene-expression analysis

Matched tumour–normal samples are analysed using DESeq2 with a paired design:

`~ patient + condition`

The analysis accounts for patient-specific differences when estimating tumour-associated expression changes. The resulting differential-expression results provide a biological starting point for downstream investigation.

### 2. Machine-learning subtype classification

A multinomial LASSO model implemented using `glmnet` classifies samples into PAM50-labelled molecular subtypes using variance-selected gene-expression features.

The benchmark uses:

- A stratified training/test split of 877 training samples and 218 held-out test samples.
- Five-fold cross-validation within the training cohort.
- Feature selection performed within each cross-validation fold.
- Final model evaluation on the held-out test set.
- Accuracy, balanced accuracy, macro F1, confusion matrices, and subtype-level performance metrics.

The corrected benchmark is designed to prevent validation-fold samples from influencing feature selection during cross-validation.

### 3. Subtype-associated expression signatures

Subtype signatures are constructed and evaluated to investigate molecular differences among breast cancer subtypes. Additional analyses assess the stability of these signatures and the extent to which analytical choices affect downstream findings.

### 4. Drug-response association and prioritisation

Expression signatures are evaluated against drug-response measurements from GDSC2 cell lines. Drug-level associations are used to explore candidate prioritisation and compare results across molecular subtypes.

These analyses assess associations between expression-derived scores and experimental cell-line response measurements. They do not establish clinical treatment efficacy or demonstrate that a model can predict an individual patient's response to therapy.

### 5. Robustness and sensitivity analyses

The project evaluates whether prioritisation results remain consistent when analytical assumptions change, including:

- Minimum cell-line coverage thresholds.
- Alternative drug-response metrics, including AUC and LN_IC50.
- Restriction to Broad-source cell lines.
- Alternative signature-scoring approaches and robustness analyses.

Ranking correlations, top-ranked candidate overlap, effect sizes, and direction concordance are used to assess the stability of the findings.

## Preliminary machine-learning results

The corrected multinomial LASSO benchmark produced the following held-out test results:

| Metric | Result |
|---|---:|
| Test samples | 218 |
| Selected genes | 1,000 |
| Accuracy | 90.83% |
| Balanced accuracy | 80.33% |
| Macro F1 | 0.830 |

Performance differed across subtypes. Basal and Luminal A samples were classified more consistently, whereas Normal-like classification remained challenging in the held-out test set.

These values describe the current benchmark on the specified split. They should be interpreted alongside the confusion matrix, class-level metrics, and limitations associated with the available subtype labels and class sizes.

## Reproducibility

The project is organised as a sequence of numbered R scripts. Each script documents its analytical purpose, inputs, outputs, and relevant methodological decisions.

The repository is intended to provide:

- Reproducible analysis scripts.
- Documentation of data provenance and preprocessing.
- Selected derived results and figures.
- Model evaluation and sensitivity-analysis summaries.
- Information about the R environment and package dependencies.

Large raw datasets and temporary intermediate files are not required to be stored in the Git repository. Users should obtain the original datasets from their official sources and follow the documented script order.

## Limitations

- Public transcriptomic datasets and cell-line drug-response data represent different experimental settings.
- Associations between expression signatures and drug response do not establish causality.
- Cell-line findings cannot be directly interpreted as clinical treatment recommendations.
- Subtype classification performance depends on the supplied labels, feature-selection procedure, and composition of the evaluation cohort.
- Results from drug prioritisation depend on data coverage, response metrics, and analytical assumptions.
- Computationally prioritised candidates require independent validation before biological or therapeutic conclusions can be established.

## Technologies

- **R and RStudio**
- **DESeq2** — differential gene-expression analysis
- **glmnet** — regularised multinomial classification
- **Bioconductor** — transcriptomic data analysis
- **TCGA-BRCA, PAM50 and GDSC2** — public molecular and drug-response data
- **Statistical evaluation and sensitivity analysis** — model performance, ranking consistency, and robustness

## Project status

The project includes paired differential-expression analysis, subtype-classification benchmarking, expression-signature analyses, drug-response association analyses, and sensitivity analyses. Repository organisation, selected figures, final documentation, and reproducibility checks are being consolidated.

## Scope

This is an exploratory computational biology project integrating cancer transcriptomics, machine learning, and expression-informed drug prioritisation. Its purpose is to investigate molecular patterns and the stability of computational findings, not to provide clinical decision support.

## Author

**Kazhal Ghiasi**

