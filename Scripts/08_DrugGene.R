# 08_DrugGene.R
# Existing evidence audit, DGIdb query, and final drug-gene integration.
#
# Finalized implementation sources:
#   ../17_Drug_Gene_Therapeutic_Relevance.py
#   ../17_Fresh_Drug_Gene_Query_Integration.py
#
# Drug-gene interaction evidence is not interpreted as proof of therapeutic
# efficacy, sensitivity, resistance, or cancer-specific clinical actionability.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "08_DrugGene.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
existing_audit_script <- file.path(repo_root, "17_Drug_Gene_Therapeutic_Relevance.py")
fresh_query_script <- file.path(repo_root, "17_Fresh_Drug_Gene_Query_Integration.py")

message("Existing drug-gene audit source: ", existing_audit_script)
message("Fresh DGIdb integration source: ", fresh_query_script)

if (!file.exists(existing_audit_script) || !file.exists(fresh_query_script)) {
  stop("Missing one or more finalized S17 Python workflows.")
}

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning S17 drug-gene workflows.")
  quit(save = "no", status = 0)
}

status1 <- system2("python", existing_audit_script)
if (!identical(status1, 0L)) {
  stop("S17 existing drug-gene audit failed with status ", status1)
}

status2 <- system2("python", fresh_query_script)
if (!identical(status2, 0L)) {
  stop("S17 fresh DGIdb query integration failed with status ", status2)
}
