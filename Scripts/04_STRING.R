# 04_STRING.R
# STRING v12 API network construction and statistics.
#
# Inputs:
#   work/S11_Hub_Genes_MCC.csv or configured DEG gene sets.
# Outputs:
#   STRING edge/node exports and network statistics under outputs/STRING/.
#
# This script is network-enabled only when RUN_ANALYSIS=TRUE.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning STRING network construction.")
  quit(save = "no", status = 0)
}

required_packages <- c("readr", "dplyr", "httr2", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(httr2)
  library(jsonlite)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "04_STRING.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
hub_path <- file.path(repo_root, "work", "S11_Hub_Genes_MCC.csv")
out_dir <- file.path(repo_root, "outputs", "STRING")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

string_required_score <- 700
string_version <- "12"
species <- 9606

hubs <- readr::read_csv(hub_path, show_col_types = FALSE)

query_string_network <- function(genes) {
  request("https://string-db.org/api/tsv/network") |>
    req_url_query(
      identifiers = paste(unique(genes), collapse = "%0d"),
      species = species,
      required_score = string_required_score
    ) |>
    req_perform() |>
    resp_body_string()
}

for (cancer in unique(hubs$Cancer_Type)) {
  genes <- hubs$Gene_Symbol[hubs$Cancer_Type == cancer]
  response <- query_string_network(genes)
  safe_name <- gsub("[^A-Za-z0-9]+", "_", cancer)
  writeLines(response, file.path(out_dir, paste0(safe_name, "_STRING_network.tsv")))
}

metadata <- data.frame(
  STRING_Version = string_version,
  Required_Score = string_required_score,
  Species = species,
  Query_Date = Sys.time()
)
readr::write_csv(metadata, file.path(out_dir, "STRING_query_metadata.csv"))
