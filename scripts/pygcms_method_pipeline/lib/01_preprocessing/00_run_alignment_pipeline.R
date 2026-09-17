#!/usr/bin/env Rscript
# =============================================================================
# Purpose: Orchestrate sample grouping and within- and across-group alignment.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================


resolve_pipeline_dir <- function() {
  explicit_dir <- Sys.getenv("PIPELINE_DIR", unset = "")
  if (nzchar(trimws(explicit_dir))) {
    return(normalizePath(explicit_dir, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot resolve 00_run_alignment_pipeline.R location.")
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
}

pipeline_dir <- resolve_pipeline_dir()
module_paths <- file.path(pipeline_dir, c(
  "01_internal_standard_guided_sample_grouping.R",
  "06_alignid0_singleton_merge_core.R",
  "02_within_group_deconvolution_and_peak_alignment.R",
  "03_global_peak_alignment_across_sample_groups.R"
))
missing_modules <- module_paths[!file.exists(module_paths)]
if (length(missing_modules)) {
  stop("Missing alignment module(s): ", paste(missing_modules, collapse = ", "))
}
for (module_path in module_paths) sys.source(module_path, envir = globalenv())

config <- list(
  project_dir = Sys.getenv("PROJECT_DIR", unset = getwd()),
  cdf_dir = Sys.getenv("CDF_DIR", unset = ""),
  out_dir = Sys.getenv("OUT_DIR", unset = ""),
  preblock_dir = Sys.getenv("PREBLOCK_DIR", unset = ""),
  alignment_mode = Sys.getenv("ALIGNMENT_MODE", unset = "auto"),
  single_group_threshold = as.integer(Sys.getenv("ALIGNMENT_SINGLE_BLOCK_THRESHOLD", unset = "25")),
  single_group_time_dist_sec = as.numeric(Sys.getenv("ALIGNMENT_SINGLE_BLOCK_TIME_DIST_SEC", unset = "60"))
)

grouping <- run_internal_standard_guided_sample_grouping(config)

within_group <- run_within_group_deconvolution_and_peak_alignment(config, grouping)

global_alignment <- run_global_peak_alignment_across_sample_groups(config, within_group)
invisible(global_alignment)
