# 05_cytoHubba.R
# MCC hub prioritization plus pan-cancer recurrence and TF enrichment handoff.
#
# Notes:
#   Cytoscape/cytoHubba ranking is GUI/plugin-derived in the finalized workflow.
#   This repository-facing script standardizes downstream hub prioritization
#   inputs and records the finalized S11/S13 source files.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "05_cytoHubba.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
s11_path <- file.path(repo_root, "work", "S11_Hub_Genes_MCC.csv")
s13b_path <- file.path(repo_root, "work", "S13B_Conserved_Candidates_With_MCC_Hub_Support.csv")

if (!file.exists(s11_path)) {
  stop("Missing finalized S11 hub table: ", s11_path)
}
if (!file.exists(s13b_path)) {
  stop("Missing finalized S13B conserved candidate table: ", s13b_path)
}

message("Finalized MCC hubs: ", s11_path)
message("Finalized pan-cancer recurrence evidence: ", s13b_path)

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning cytoHubba/TF enrichment.")
  quit(save = "no", status = 0)
}

required_packages <- c("readr", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

out_dir <- file.path(repo_root, "outputs", "hub_prioritization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s11 <- readr::read_csv(s11_path, show_col_types = FALSE)
s13b <- readr::read_csv(s13b_path, show_col_types = FALSE)

summary <- s11 |>
  group_by(Cancer_Type, Dataset, Regulation) |>
  summarise(
    Finalized_MCC_Hub_Count = n_distinct(Gene_Symbol),
    Best_MCC_Rank = min(MCC_Rank, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(summary, file.path(out_dir, "S05_MCC_hub_summary.csv"))
readr::write_csv(s13b, file.path(out_dir, "S05_conserved_candidate_reference.csv"))
