# 07_scRNA_validation.R
# TISCH candidate prioritization and manually curated single-cell integration.
#
# Finalized implementation sources:
#   ../16_TISCH_Candidate_Prioritization.py
#   ../16_TISCH_Final_Integration.py
#
# TISCH values were manually curated after candidate selection. They are used
# only for localization/contextualization, not to select candidate genes.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "07_scRNA_validation.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
candidate_script <- file.path(repo_root, "16_TISCH_Candidate_Prioritization.py")
integration_script <- file.path(repo_root, "16_TISCH_Final_Integration.py")

message("scRNA candidate prioritization source: ", candidate_script)
message("Final TISCH integration source: ", integration_script)

if (!file.exists(candidate_script) || !file.exists(integration_script)) {
  stop("Missing one or more finalized S16 Python workflows.")
}

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning S16 single-cell workflows.")
  quit(save = "no", status = 0)
}

status1 <- system2("python", candidate_script)
if (!identical(status1, 0L)) {
  stop("S16 candidate prioritization failed with status ", status1)
}

status2 <- system2("python", integration_script)
if (!identical(status2, 0L)) {
  stop("S16 final TISCH integration failed with status ", status2)
}
