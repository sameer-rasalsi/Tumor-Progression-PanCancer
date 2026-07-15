# 03_GO_KEGG.R
# Final g:Profiler GO/KEGG enrichment workflow.
#
# Inputs:
#   Final DEG/hub gene tables from S11-compatible outputs.
# Outputs:
#   If rerun, enrichment tables are written to outputs/enrichment/.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning GO/KEGG enrichment.")
  quit(save = "no", status = 0)
}

required_packages <- c("gprofiler2", "readr", "dplyr", "purrr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(gprofiler2)
  library(readr)
  library(dplyr)
  library(purrr)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "03_GO_KEGG.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
hub_path <- file.path(repo_root, "work", "S11_Hub_Genes_MCC.csv")
out_dir <- file.path(repo_root, "outputs", "enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hubs <- readr::read_csv(hub_path, show_col_types = FALSE)

run_profile <- function(df) {
  genes <- unique(df$Gene_Symbol)
  if (length(genes) < 2) {
    return(NULL)
  }
  gost(
    query = genes,
    organism = "hsapiens",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    correction_method = "g_SCS"
  )$result
}

results <- hubs |>
  split(list(hubs$Cancer_Type, hubs$Regulation), drop = TRUE) |>
  purrr::imap(function(df, label) {
    res <- run_profile(df)
    if (!is.null(res)) {
      res$Analysis_Set <- label
    }
    res
  }) |>
  purrr::compact() |>
  dplyr::bind_rows()

readr::write_csv(results, file.path(out_dir, "S03_gProfiler_GO_KEGG_results.csv"))
