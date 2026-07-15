# 02_limma.R
# Generalized microarray DEG workflow for breast, lung, and melanoma cohorts.
#
# Purpose:
#   Repository-facing limma workflow with a dataset configuration table that
#   preserves the finalized analysis class and criteria for each microarray
#   discovery cohort.
#
# Outputs:
#   If rerun, results are written under outputs/limma/.

dataset_config <- data.frame(
  dataset = c("GSE173839", "GSE33072", "GSE99898"),
  cancer = c("Breast cancer", "Lung cancer", "Melanoma"),
  analysis_class = c(
    "Primary_Strict_DEGs",
    "Exploratory_Candidates",
    "Exploratory_Candidates"
  ),
  finalized_criterion = c(
    "BH adjusted P < 0.05; |log2FC| >= 1.0",
    "Exploratory finalized lung criterion as documented in S11/S14 source outputs",
    "Exploratory finalized melanoma criterion as documented in S11/S14 source outputs"
  ),
  stringsAsFactors = FALSE
)

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")

message("Configured microarray datasets:")
print(dataset_config)

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning limma biological analyses.")
  quit(save = "no", status = 0)
}

required_packages <- c("limma", "readr", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(tibble)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "02_limma.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
out_dir <- file.path(repo_root, "outputs", "limma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_limma_dataset <- function(config_row) {
  dataset <- config_row[["dataset"]]
  expr_path <- file.path(repo_root, "data", "microarray", paste0(dataset, "_expression.tsv"))
  metadata_path <- file.path(repo_root, "data", "microarray", paste0(dataset, "_metadata.tsv"))

  if (!file.exists(expr_path) || !file.exists(metadata_path)) {
    stop("Missing expression or metadata input for ", dataset,
         ". Expected: ", expr_path, " and ", metadata_path)
  }

  expr <- readr::read_tsv(expr_path, show_col_types = FALSE)
  metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE)
  gene_col <- names(expr)[1]
  expr_matrix <- expr |>
    tibble::column_to_rownames(gene_col) |>
    as.matrix()

  if (!all(c("sample_id", "condition") %in% names(metadata))) {
    stop("Metadata for ", dataset, " must contain sample_id and condition columns.")
  }

  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sample_id
  expr_matrix <- expr_matrix[, rownames(metadata), drop = FALSE]

  design <- model.matrix(~ 0 + condition, data = metadata)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(expr_matrix, design)
  if (ncol(design) != 2) {
    stop("Expected two-condition design for ", dataset)
  }
  contrast <- makeContrasts(contrasts = paste(colnames(design), collapse = "-"), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  res <- topTable(fit2, number = Inf, adjust.method = "BH")
  res$Gene_Symbol <- rownames(res)
  readr::write_csv(res, file.path(out_dir, paste0(dataset, "_limma_results.csv")))
}

invisible(apply(dataset_config, 1, run_limma_dataset))
