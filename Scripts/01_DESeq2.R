# 01_DESeq2.R
# Generalized RNA-seq DEG workflow for ovarian, prostate, and SCLC cohorts.
#
# Purpose:
#   Repository-facing script for RNA-seq differential expression analysis.
#   Dataset-specific thresholds are documented in one configuration table so the
#   final manuscript does not imply that all RNA-seq cohorts used identical
#   criteria.
#
# Inputs:
#   Raw/read-count matrices and sample metadata are expected under data/ when
#   the analysis is rerun by the user.
#
# Outputs:
#   DEG result tables should be written under outputs/DESeq2/ if rerun.
#
# Reproducibility note:
#   This script is a generalized workflow scaffold. It preserves the finalized
#   criteria used downstream and does not rerun analyses unless RUN_ANALYSIS=TRUE.

dataset_config <- data.frame(
  dataset = c("GSE156699", "GSE279730", "GSE188705"),
  cancer = c("Ovarian cancer", "Prostate cancer", "Small-cell lung cancer"),
  analysis_class = c(
    "Primary_Strict_DEGs",
    "Primary_Strict_DEGs",
    "Exploratory_Candidates"
  ),
  finalized_criterion = c(
    "BH adjusted P < 0.05; |log2FC| >= 1.0",
    "BH adjusted P < 0.05; |log2FC| >= 1.0",
    "Exploratory finalized SCLC criterion as documented in S11/S14 source outputs"
  ),
  stringsAsFactors = FALSE
)

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")

message("Configured RNA-seq datasets:")
print(dataset_config)

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning DESeq2 biological analyses.")
  quit(save = "no", status = 0)
}

required_packages <- c("DESeq2", "readr", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "01_DESeq2.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
out_dir <- file.path(repo_root, "outputs", "DESeq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_deseq2_dataset <- function(config_row) {
  dataset <- config_row[["dataset"]]
  count_path <- file.path(repo_root, "data", "RNAseq", paste0(dataset, "_counts.tsv"))
  metadata_path <- file.path(repo_root, "data", "RNAseq", paste0(dataset, "_metadata.tsv"))

  if (!file.exists(count_path) || !file.exists(metadata_path)) {
    stop("Missing count or metadata input for ", dataset,
         ". Expected: ", count_path, " and ", metadata_path)
  }

  counts <- readr::read_tsv(count_path, show_col_types = FALSE)
  metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE)
  gene_col <- names(counts)[1]
  count_matrix <- counts |>
    tibble::column_to_rownames(gene_col) |>
    as.matrix()

  if (!"condition" %in% names(metadata)) {
    stop("Metadata for ", dataset, " must contain a condition column.")
  }
  if (!"sample_id" %in% names(metadata)) {
    stop("Metadata for ", dataset, " must contain a sample_id column.")
  }

  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sample_id
  count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = metadata,
    design = ~ condition
  )
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  res$Gene_Symbol <- rownames(res)
  readr::write_csv(res, file.path(out_dir, paste0(dataset, "_DESeq2_results.csv")))
}

invisible(apply(dataset_config, 1, run_deseq2_dataset))
