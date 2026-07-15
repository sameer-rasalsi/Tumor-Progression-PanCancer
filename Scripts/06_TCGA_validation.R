# 06_TCGA_validation.R
# TCGA tumor-normal, Kaplan-Meier, Cox, and stage validation.
#
# Finalized implementation source:
#   ../15_TCGA_Clinical_Validation.py
#
# This R-facing script preserves the supervisor's requested filename while
# invoking the validated Python workflow only when RUN_ANALYSIS=TRUE.

run_analysis <- identical(Sys.getenv("RUN_ANALYSIS"), "TRUE")
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)][1])
if (is.na(script_file)) {
  script_file <- file.path(getwd(), "Scripts", "06_TCGA_validation.R")
}
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
python_script <- file.path(repo_root, "15_TCGA_Clinical_Validation.py")
output_dir <- file.path(repo_root, "S15_TCGA_clinical_validation")

message("TCGA validation source script: ", python_script)
message("Finalized TCGA output directory: ", output_dir)

if (!file.exists(python_script)) {
  stop("Missing finalized TCGA Python workflow: ", python_script)
}

if (!run_analysis) {
  message("RUN_ANALYSIS is not TRUE; not rerunning TCGA validation.")
  quit(save = "no", status = 0)
}

status <- system2("python", python_script)
if (!identical(status, 0L)) {
  stop("TCGA validation Python workflow failed with status ", status)
}
