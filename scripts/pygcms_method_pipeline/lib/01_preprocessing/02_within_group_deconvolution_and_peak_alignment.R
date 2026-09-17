# =============================================================================
# Purpose: Deconvolve chromatograms and align cleaned peaks within sample groups.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

# Within-group deconvolution and peak alignment
# Self-contained stage module: no other script is loaded or executed here.

filter_peak_list_by_area_fraction <- function(peak_list, threshold, block_name) {
  if (!is.finite(threshold) || threshold < 0 || threshold >= 1) {
    stop("AREA_FRACTION_FILTER must be finite and in [0, 1), got: ", threshold)
  }
  sample_names <- names(peak_list)
  if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_along(peak_list))
  audits <- vector("list", length(peak_list))
  summaries <- vector("list", length(peak_list))
  filtered <- peak_list
  for (i in seq_along(peak_list)) {
    peaks <- as.data.frame(peak_list[[i]], stringsAsFactors = FALSE)
    if (!all(c("ID", "RT", "Area") %in% names(peaks))) {
      stop(
        "Missing ID, RT, or Area column for ",
        block_name, "/", sample_names[[i]]
      )
    }
    areas <- suppressWarnings(as.numeric(as.character(peaks$Area)))
    if (!nrow(peaks) || any(!is.finite(areas)) || any(areas < 0) || sum(areas) <= 0) {
      stop("Invalid post-cleanup Area values for ", block_name, "/", sample_names[[i]])
    }
    total_area <- sum(areas)
    fractions <- areas / total_area
    retained <- fractions >= threshold
    filtered[[i]] <- peak_list[[i]][retained, , drop = FALSE]
    audits[[i]] <- data.frame(
      rt_block = block_name,
      sample = sample_names[[i]],
      peak_id = peaks$ID,
      RT = peaks$RT,
      Area = areas,
      sample_total_area = total_area,
      area_fraction = fractions,
      threshold = threshold,
      retained = retained,
      removal_reason = ifelse(retained, "", "below_within_sample_area_fraction"),
      stringsAsFactors = FALSE
    )
    summaries[[i]] <- data.frame(
      rt_block = block_name,
      sample = sample_names[[i]],
      threshold = threshold,
      n_peaks_before = nrow(peaks),
      n_removed = sum(!retained),
      n_peaks_after = sum(retained),
      area_before = total_area,
      area_removed = sum(areas[!retained]),
      area_after = sum(areas[retained]),
      area_retained_pct = 100 * sum(areas[retained]) / total_area,
      stringsAsFactors = FALSE
    )
  }
  list(
    peak_list = filtered,
    audit = dplyr::bind_rows(audits),
    summary = dplyr::bind_rows(summaries)
  )
}

run_within_group_deconvolution_and_peak_alignment <- function(config = list(), grouping) {
  suppressPackageStartupMessages({
    library(erah)
    library(dplyr)
    library(readr)
  })
  
  # =============================================================================
  # Block-wise eRah pipeline for large GC-MS datasets
  #
  # This script implements the workflow discussed for reducing RT-mismatch-driven
  # false positives:
  #   1) Pre-analysis sample blocking using raw internal-standard EIC RTs
  #   2) eRah deconvolution within each RT block
  #   3) eRah alignment within each RT block using a reduced max.time.dist
  #   4) Unified post-alignment feature consolidation within each RT block
  #   5) Export block-level aligned tables for later strict cross-block merging
  #
  # Default optimized parameters:
  #   min.peak.width  = 2.5
  #   min.peak.height = 8000
  #   noise.threshold = 1000
  #   analysis.time   = 3-90 min
  #   min.spectra.cor = 0.90
  #   mz.range        = 46-650
  #
  # This module is loaded and called only by 00_run_alignment_pipeline.R.
  # =============================================================================
  
  # ---- Project root ------------------------------------------------------------
  # This makes the script safer to run from RStudio.
  # Run the repository launcher, or set PROJECT_DIR to your checkout.
  # The repository root and its scripts folder can also be detected from getwd().
  #
  # You can also override it manually before sourcing this script:
  #   Sys.setenv(PROJECT_DIR = "/path/to/repository")
  
  resolve_project_dir <- function() {
    env_project_dir <- Sys.getenv("PROJECT_DIR", unset = "")
    if (nzchar(trimws(env_project_dir))) {
      return(normalizePath(path.expand(env_project_dir), winslash = "/", mustWork = TRUE))
    }
    
    wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    
    # Try to infer the script location when running with R command-line runner.
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    script_candidates <- character(0)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg[[1]])
      script_path <- normalizePath(script_path, winslash = "/", mustWork = FALSE)
      script_candidates <- c(dirname(script_path), dirname(dirname(script_path)))
    }
    
    candidates <- unique(c(
      wd,
      dirname(wd),
      script_candidates
    ))
    
    has_pipeline_dir <- function(path) {
      file.exists(file.path(
        path,
        "scripts",
        "pygcms_method_pipeline",
        "lib",
        "01_preprocessing",
        "00_run_alignment_pipeline.R"
      ))
    }
    
    hits <- candidates[vapply(candidates, has_pipeline_dir, logical(1))]
    if (length(hits) > 0) {
      return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
    }
    
    stop(
      "Could not determine PROJECT_DIR automatically.\n",
      "Please set PROJECT_DIR before running the pipeline:\n",
      "  Sys.setenv(PROJECT_DIR = \"/path/to/repository\")\n",
      "Current working directory was: ", wd
    )
  }
  
  project_dir <- resolve_project_dir()
  setwd(project_dir)
  message("PROJECT_DIR: ", project_dir)
  
  
  
  cdf_dir <- Sys.getenv(
    "CDF_DIR",
    unset = file.path(project_dir, "data", "raw")
  )
  
  out_dir <- Sys.getenv(
    "OUT_DIR",
    unset = file.path(project_dir, "blockwise_erah_outputs")
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  preblock_dir <- Sys.getenv(
    "PREBLOCK_DIR",
    unset = file.path(project_dir, "cdf_qc", "pre_alignment_blocks")
  )
  dir.create(preblock_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Pre-blocking parameters ------------------------------------------------
  
  n_subset <- Sys.getenv("N_SUBSET", unset = "all")
  small_sample_n <- Sys.getenv("SMALL_SAMPLE_N", unset = "all")
  small_sample_mode <- Sys.getenv("SMALL_SAMPLE_MODE", unset = "first")
  small_sample_seed <- as.integer(Sys.getenv("SMALL_SAMPLE_SEED", unset = "1"))
  min_block_size <- Sys.getenv("MIN_BLOCK_SIZE", unset = "10")
  max_block_size <- Sys.getenv("MAX_BLOCK_SIZE", unset = "25")
  target_max_time_dist <- Sys.getenv("TARGET_MAX_TIME_DIST", unset = "30")
  requested_alignment_mode <- Sys.getenv("ALIGNMENT_MODE", unset = "auto")
  single_block_threshold <- as.integer(Sys.getenv(
    "ALIGNMENT_SINGLE_BLOCK_THRESHOLD", unset = "25"
  ))
  single_block_time_dist_sec <- as.numeric(Sys.getenv(
    "ALIGNMENT_SINGLE_BLOCK_TIME_DIST_SEC", unset = "60"
  ))
  
  # ---- eRah optimized processing parameters ----------------------------------
  
  min_peak_width <- as.numeric(Sys.getenv("MIN_PEAK_WIDTH", unset = "2.5"))
  min_peak_height <- as.numeric(Sys.getenv("MIN_PEAK_HEIGHT", unset = "8000"))
  area_fraction_filter <- as.numeric(Sys.getenv("AREA_FRACTION_FILTER", unset = "0"))
  if (!is.finite(area_fraction_filter) ||
      area_fraction_filter < 0 ||
      area_fraction_filter >= 1) {
    stop(
      "AREA_FRACTION_FILTER must be finite and in [0, 1), got: ",
      area_fraction_filter
    )
  }
  noise_threshold <- as.numeric(Sys.getenv("NOISE_THRESHOLD", unset = "1000"))
  analysis_start <- as.numeric(Sys.getenv("ANALYSIS_START_MIN", unset = "3"))
  analysis_end <- as.numeric(Sys.getenv("ANALYSIS_END_MIN", unset = "90"))
  
  parse_mz_ranges <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) {
      stop("AVOID_PROCESSING_MZ cannot be empty.")
    }
    parts <- trimws(unlist(strsplit(x, ",")))
    parts <- parts[nzchar(parts)]
    vals <- unlist(lapply(parts, function(part) {
      if (grepl(":", part, fixed = TRUE)) {
        bounds <- suppressWarnings(as.integer(trimws(strsplit(part, ":", fixed = TRUE)[[1]])))
        if (length(bounds) != 2 || any(is.na(bounds))) {
          stop("Invalid AVOID_PROCESSING_MZ range: ", part)
        }
        if (bounds[[2]] < bounds[[1]]) {
          stop("Invalid descending AVOID_PROCESSING_MZ range: ", part)
        }
        if ((bounds[[2]] - bounds[[1]] + 1L) > 20L) {
          stop(
            "Refusing very large AVOID_PROCESSING_MZ range: ", part,
            ". Expected small ranges such as 73:75 or 147:149."
          )
        }
        seq(bounds[[1]], bounds[[2]])
      } else {
        val <- suppressWarnings(as.integer(part))
        if (is.na(val)) stop("Invalid AVOID_PROCESSING_MZ value: ", part)
        val
      }
    }))
    vals <- sort(unique(as.integer(vals)))
    if (length(vals) > 30L) {
      stop(
        "Refusing ", length(vals), " avoid-processing m/z values: ",
        paste(vals, collapse = ","),
        ". This is probably a malformed AVOID_PROCESSING_MZ setting."
      )
    }
    vals
  }
  
  avoid_processing_mz_spec <- Sys.getenv("AVOID_PROCESSING_MZ", unset = "74:75,147:149,207,281")
  avoid_processing_mz <- parse_mz_ranges(avoid_processing_mz_spec)
  
  min_spectra_cor <- as.numeric(Sys.getenv("MIN_SPECTRA_COR", unset = "0.90"))
  mz_min <- as.integer(Sys.getenv("MZ_MIN", unset = "46"))
  mz_max <- as.integer(Sys.getenv("MZ_MAX", unset = "650"))
  
  # Pre-alignment over-deconvolution duplicate QC.
  # This step is intentionally conservative because it removes weaker peaks before
  # alignment. Profile correlation is calculated for audit only, but it is not used
  # as a deletion criterion. Final tested rule:
  #   within the same sample: RT diff <= 1 s, spectral cosine >= 0.95,
  #   BaseMz diff <= 1. The larger-area peak is retained.
  remove_split_peaks_pre_alignment <- tolower(Sys.getenv(
    "REMOVE_SPLIT_PEAKS_PRE_ALIGNMENT",
    unset = "true"
  )) %in% c("1", "true", "yes")
  pre_alignment_overdec_rt_sec <- as.numeric(Sys.getenv("PRE_ALIGNMENT_OVERDEC_RT_SEC", unset = "1.0"))
  pre_alignment_overdec_spectral_cos <- as.numeric(Sys.getenv("PRE_ALIGNMENT_OVERDEC_SPECTRAL_COS", unset = "0.95"))
  pre_alignment_overdec_base_mz_tol <- as.numeric(Sys.getenv("PRE_ALIGNMENT_OVERDEC_BASE_MZ_TOL", unset = "1"))
  
  # Unified post-alignment feature consolidation. This replaces the former two-step
  # cross-sample split merge + shoulder/residual consolidation in the main workflow.
  # Final tested rule within each RT block:
  #   RT diff <= 60 s, spectral cosine >= 0.90, BaseMz diff <= 1.
  # Within each final feature, sample-wise areas are summed, heights are taken as
  # maxima, and tmean is recalculated as an area-weighted mean.
  post_alignment_feature_consolidation <- tolower(Sys.getenv(
    "POST_ALIGNMENT_FEATURE_CONSOLIDATION",
    unset = "true"
  )) %in% c("1", "true", "yes")
  final_feature_merge_rt_sec <- as.numeric(Sys.getenv("FINAL_FEATURE_MERGE_RT_SEC", unset = "60"))
  final_feature_merge_spectral_cos <- as.numeric(Sys.getenv("FINAL_FEATURE_MERGE_SPECTRAL_COS", unset = "0.90"))
  final_feature_merge_base_mz_tol <- as.numeric(Sys.getenv("FINAL_FEATURE_MERGE_BASE_MZ_TOL", unset = "1"))
  final_feature_group_rt_span_sec <- as.numeric(Sys.getenv(
    "FINAL_FEATURE_GROUP_RT_SPAN_SEC",
    unset = as.character(final_feature_merge_rt_sec)
  ))
  final_feature_group_base_mz_span <- as.numeric(Sys.getenv(
    "FINAL_FEATURE_GROUP_BASE_MZ_SPAN",
    unset = as.character(final_feature_merge_base_mz_tol)
  ))
  final_feature_alignment_split_jaccard <- as.numeric(Sys.getenv(
    "FINAL_FEATURE_ALIGNMENT_SPLIT_JACCARD",
    unset = "0.15"
  ))

  # AlignID=0 is an alignment status, not a deletion criterion. Every such peak
  # is first offered to an existing block feature and is otherwise retained in
  # conservative same-block singleton clustering.
  alignid0_rescue_enabled <- tolower(Sys.getenv(
    "ALIGNID0_RESCUE_ENABLED", unset = "true"
  )) %in% c("1", "true", "yes")
  singleton_thresholds <- normalize_singleton_thresholds(list(
    max_rt_diff_sec = as.numeric(Sys.getenv(
      "ALIGNID0_WITHIN_BLOCK_RT_SEC", unset = as.character(final_feature_merge_rt_sec)
    )),
    min_spectral_cosine = as.numeric(Sys.getenv(
      "ALIGNID0_MIN_SPECTRAL_COSINE", unset = "0.90"
    )),
    base_mz_tol = as.numeric(Sys.getenv("ALIGNID0_BASE_MZ_TOL", unset = "1")),
    ambiguity_margin = as.numeric(Sys.getenv("ALIGNID0_AMBIGUITY_MARGIN", unset = "0.01"))
  ))
  
  # ALIGN_TIME_DIST can be:
  #   auto: use each block's recommended_max_time_dist_sec
  #   a number: use that fixed value for all blocks, e.g. 30
  align_time_dist_mode <- Sys.getenv("ALIGN_TIME_DIST", unset = "60")
  use_existing_block_folders <- tolower(Sys.getenv(
    "USE_EXISTING_BLOCK_FOLDERS",
    unset = "false"
  )) %in% c("1", "true", "yes")
  use_existing_preblock <- tolower(Sys.getenv(
    "USE_EXISTING_PREBLOCK",
    unset = "false"
  )) %in% c("1", "true", "yes")
  
  # Annotation is optional because the first goal is to test block-wise alignment.
  do_annotation <- tolower(Sys.getenv("DO_ANNOTATION", unset = "false")) %in% c("1", "true", "yes")
  msp_file <- Sys.getenv(
    "MSP_FILE",
    unset = file.path(project_dir, "external", "MassBank_NIST_20241126_rev.msp")
  )
  n_putative <- as.integer(Sys.getenv("N_PUTATIVE", unset = "3"))
  
  message("=== Block-wise eRah pipeline ===")
  message("CDF_DIR: ", cdf_dir)
  message("OUT_DIR: ", out_dir)
  message("PREBLOCK_DIR: ", preblock_dir)
  message("N_SUBSET: ", n_subset)
  message("SMALL_SAMPLE_N: ", small_sample_n)
  message("SMALL_SAMPLE_MODE: ", small_sample_mode)
  message("Block size: ", min_block_size, "-", max_block_size)
  message("Target max.time.dist: ", target_max_time_dist, " sec")
  message("Requested alignment mode: ", requested_alignment_mode)
  message("Single-block threshold: ", single_block_threshold)
  message("Single-block max.time.dist: ", single_block_time_dist_sec, " sec")
  message("ALIGN_TIME_DIST mode: ", align_time_dist_mode)
  message("USE_EXISTING_BLOCK_FOLDERS: ", use_existing_block_folders)
  message("USE_EXISTING_PREBLOCK: ", use_existing_preblock)
  message("REMOVE_SPLIT_PEAKS_PRE_ALIGNMENT: ", remove_split_peaks_pre_alignment)
  message("PRE_ALIGNMENT_OVERDEC_RT_SEC: ", pre_alignment_overdec_rt_sec)
  message("PRE_ALIGNMENT_OVERDEC_SPECTRAL_COS: ", pre_alignment_overdec_spectral_cos)
  message("PRE_ALIGNMENT_OVERDEC_BASE_MZ_TOL: ", pre_alignment_overdec_base_mz_tol)
  message("POST_ALIGNMENT_FEATURE_CONSOLIDATION: ", post_alignment_feature_consolidation)
  message("FINAL_FEATURE_MERGE_RT_SEC: ", final_feature_merge_rt_sec)
  message("FINAL_FEATURE_MERGE_SPECTRAL_COS: ", final_feature_merge_spectral_cos)
  message("FINAL_FEATURE_MERGE_BASE_MZ_TOL: ", final_feature_merge_base_mz_tol)
  message("AVOID_PROCESSING_MZ: ", paste(avoid_processing_mz, collapse = ","))
  
  # ---- Step 1: consume completed sample grouping ------------------------------
  selected_alignment_mode <- grouping$selected_alignment_mode
  if (identical(selected_alignment_mode, "single")) align_time_dist_mode <- "auto"
  message("[1/5] Using sample groups prepared by module 01.")
  block_summary_file <- file.path(preblock_dir, "auto_block_summary.csv")
  rt_by_sample_file <- file.path(preblock_dir, "auto_block_raw_internal_standard_rt_by_sample.csv")
  if (!file.exists(block_summary_file)) stop("Missing block summary: ", block_summary_file)
  if (!file.exists(rt_by_sample_file)) stop("Missing RT-by-sample file: ", rt_by_sample_file)
  
  block_summary <- read.csv(block_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  rt_by_sample <- read.csv(rt_by_sample_file, stringsAsFactors = FALSE, check.names = FALSE)
  blocks <- block_summary$rt_block
  only_block <- Sys.getenv("ONLY_BLOCK", unset = "")
  if (nzchar(trimws(only_block))) {
    requested_blocks <- trimws(unlist(strsplit(only_block, ",")))
    requested_blocks <- requested_blocks[nzchar(requested_blocks)]
    missing_blocks <- setdiff(requested_blocks, blocks)
    if (length(missing_blocks) > 0) {
      stop("ONLY_BLOCK contains unknown block(s): ", paste(missing_blocks, collapse = ", "))
    }
    blocks <- requested_blocks
  }
  
  write.csv(block_summary, file.path(out_dir, "00_block_summary.csv"), row.names = FALSE)
  write.csv(rt_by_sample, file.path(out_dir, "00_internal_standard_rt_by_sample.csv"), row.names = FALSE)
  
  message("Blocks generated: ", length(blocks))
  print(block_summary[, intersect(c(
    "rt_block", "n_samples", "max_internal_standard_range_sec",
    "recommended_max_time_dist_sec"
  ), names(block_summary))])
  
  # ---- Optional database import ----------------------------------------------
  
  id_database <- NULL
  if (do_annotation) {
    if (!file.exists(msp_file)) stop("MSP file not found: ", msp_file)
    message("\n[2/5] Importing MSP database...")
    id_database <- importMSP(
      file = msp_file,
      DB.name = "MassBank_NIST",
      DB.version = "2024.11",
      DB.info = "Full-library search; curation required"
    )
    save(id_database, file = file.path(out_dir, "00_id_database.rda"))
  } else {
    message("\n[2/5] Skipping annotation database import (DO_ANNOTATION=false).")
  }
  
  # ---- Helpers ----------------------------------------------------------------
  
  make_experiment <- function(cdf_files, block_name) {
    instrumental <- createInstrumentalTable(cdf_files)
    
    # Keep biological labels as metadata only. They should not define RT blocks.
    sample_names <- tools::file_path_sans_ext(basename(cdf_files))
    sample_group <- rep("Sample", length(cdf_files))
    phenotype <- createPhenoTable(cdf_files, sample_group)
    
    newExp(
      instrumental = instrumental,
      phenotype = phenotype,
      info = paste0("Block-wise eRah analysis: ", block_name)
    )
  }
  
  export_deconvolved_peaks <- function(exp_dec, block_name, file) {
    peak_list <- exp_dec@Data@FactorList
    sample_names <- names(peak_list)
    out <- bind_rows(lapply(seq_along(peak_list), function(i) {
      x <- as.data.frame(peak_list[[i]], stringsAsFactors = FALSE)
      if (nrow(x) == 0) return(data.frame())
      x$sample <- sample_names[i]
      x$rt_block <- block_name
      x
    }))
    write.csv(out, file, row.names = FALSE)
    out
  }
  
  parse_profile <- function(profile) {
    if (is.na(profile) || trimws(profile) == "") {
      return(data.frame(rt = numeric(0), y = numeric(0)))
    }
    pairs <- unlist(strsplit(trimws(profile), "\\s+"))
    spl <- strsplit(pairs, ",", fixed = TRUE)
    rt <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)
    )))
    y <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)
    )))
    keep <- !is.na(rt) & !is.na(y)
    data.frame(rt = rt[keep], y = y[keep])
  }
  
  profile_cor <- function(p1, p2) {
    if (nrow(p1) < 3 || nrow(p2) < 3) return(NA_real_)
    rmin <- max(min(p1$rt), min(p2$rt))
    rmax <- min(max(p1$rt), max(p2$rt))
    if (!is.finite(rmin) || !is.finite(rmax) || rmax <= rmin) return(NA_real_)
    grid <- seq(rmin, rmax, length.out = 20)
    y1 <- approx(p1$rt, p1$y, xout = grid, rule = 2)$y
    y2 <- approx(p2$rt, p2$y, xout = grid, rule = 2)$y
    suppressWarnings(cor(y1, y2))
  }
  
  parse_spectrum <- function(spectrum) {
    if (is.na(spectrum) || trimws(spectrum) == "") {
      return(data.frame(mz = numeric(0), intensity = numeric(0)))
    }
    pairs <- unlist(strsplit(trimws(spectrum), "\\s+"))
    spl <- strsplit(pairs, ",", fixed = TRUE)
    mz <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)
    )))
    intensity <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)
    )))
    keep <- !is.na(mz) & !is.na(intensity)
    data.frame(mz = mz[keep], intensity = intensity[keep])
  }
  
  base_mz_from_spectrum <- function(spectrum) {
    x <- parse_spectrum(spectrum)
    if (nrow(x) == 0) return(NA_real_)
    x$mz[which.max(x$intensity)]
  }
  
  spectral_cosine <- function(s1, s2, mz_min = 46, mz_max = 650) {
    x1 <- parse_spectrum(s1)
    x2 <- parse_spectrum(s2)
    if (nrow(x1) == 0 || nrow(x2) == 0) return(NA_real_)
    mz_grid <- mz_min:mz_max
    v1 <- numeric(length(mz_grid))
    v2 <- numeric(length(mz_grid))
    names(v1) <- mz_grid
    names(v2) <- mz_grid
    x1 <- x1[x1$mz %in% mz_grid, , drop = FALSE]
    x2 <- x2[x2$mz %in% mz_grid, , drop = FALSE]
    if (nrow(x1) > 0) v1[as.character(x1$mz)] <- x1$intensity
    if (nrow(x2) > 0) v2[as.character(x2$mz)] <- x2$intensity
    denom <- sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2))
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    sum(v1 * v2) / denom
  }
  
  detect_split_peaks_one_sample <- function(sample_df, sample_name, block_name) {
    if (nrow(sample_df) < 2) return(data.frame())
    x <- as.data.frame(sample_df, stringsAsFactors = FALSE)
    x$.row_index <- seq_len(nrow(x))
    x$RT <- suppressWarnings(as.numeric(as.character(x$RT)))
    x$`Peak Height` <- suppressWarnings(as.numeric(as.character(x$`Peak Height`)))
    x$Area <- suppressWarnings(as.numeric(as.character(x$Area)))
    x <- x[order(x$RT, x$.row_index), , drop = FALSE]
  
    profiles <- lapply(x$Profile, parse_profile)
    rt_tol_min <- pre_alignment_overdec_rt_sec / 60
    pairs <- list()
  
    for (i in seq_len(nrow(x) - 1)) {
      candidates <- which(seq_len(nrow(x)) > i & abs(x$RT - x$RT[i]) <= rt_tol_min)
      if (length(candidates) == 0) next
  
      for (j in candidates) {
        pc <- profile_cor(profiles[[i]], profiles[[j]])
        profile_points_i <- nrow(profiles[[i]])
        profile_points_j <- nrow(profiles[[j]])
        rt_diff_sec <- abs(x$RT[j] - x$RT[i]) * 60
        sp_cos <- spectral_cosine(x$Spectra[i], x$Spectra[j], mz_min = mz_min, mz_max = mz_max)
        base_mz_i <- base_mz_from_spectrum(x$Spectra[i])
        base_mz_j <- base_mz_from_spectrum(x$Spectra[j])
        base_mz_diff <- abs(base_mz_i - base_mz_j)
  
        pre_alignment_overdec_duplicate <- (
          rt_diff_sec <= pre_alignment_overdec_rt_sec &&
            is.finite(sp_cos) &&
            sp_cos >= pre_alignment_overdec_spectral_cos &&
            is.finite(base_mz_diff) &&
            base_mz_diff <= pre_alignment_overdec_base_mz_tol
        )
        if (!pre_alignment_overdec_duplicate) next
  
        area_i <- x$Area[i]
        area_j <- x$Area[j]
        if (is.finite(area_i) && is.finite(area_j) && area_i != area_j) {
          weaker <- ifelse(area_i >= area_j, j, i)
        } else {
          weaker <- ifelse(x$`Peak Height`[i] >= x$`Peak Height`[j], j, i)
        }
        stronger <- ifelse(weaker == i, j, i)
  
        pairs[[length(pairs) + 1]] <- data.frame(
          rt_block = block_name,
          sample = sample_name,
          peak_id_1 = x$ID[i],
          peak_id_2 = x$ID[j],
          row_index_1 = x$.row_index[i],
          row_index_2 = x$.row_index[j],
          rt_1 = x$RT[i],
          rt_2 = x$RT[j],
          rt_diff_sec = rt_diff_sec,
          area_1 = x$Area[i],
          area_2 = x$Area[j],
          height_1 = x$`Peak Height`[i],
          height_2 = x$`Peak Height`[j],
          profile_cor = pc,
          spectral_cosine = sp_cos,
          base_mz_1 = base_mz_i,
          base_mz_2 = base_mz_j,
          base_mz_diff = base_mz_diff,
          profile_points_1 = profile_points_i,
          profile_points_2 = profile_points_j,
          split_reason = "pre_alignment_overdeconvolution_rt1s_cos095_base_mz_supported",
          retained_peak_id = x$ID[stronger],
          removed_peak_id = x$ID[weaker],
          removed_row_index = x$.row_index[weaker],
          stringsAsFactors = FALSE
        )
      }
    }
  
    bind_rows(pairs)
  }
  
  apply_split_peak_qc <- function(exp_dec, block_name) {
    exp_clean <- exp_dec
    peak_list <- exp_dec@Data@FactorList
    sample_names <- names(peak_list)
    if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_along(peak_list))
    
    audit <- bind_rows(lapply(seq_along(peak_list), function(i) {
      detect_split_peaks_one_sample(peak_list[[i]], sample_names[i], block_name)
    }))
    
    summary <- data.frame(
      rt_block = block_name,
      sample = sample_names,
      n_peaks_before = vapply(peak_list, nrow, integer(1)),
      n_split_pairs = 0L,
      n_removed_pre_alignment = 0L,
      n_peaks_after = vapply(peak_list, nrow, integer(1)),
      stringsAsFactors = FALSE
    )
    
    if (nrow(audit) > 0) {
      for (i in seq_along(peak_list)) {
        sn <- sample_names[i]
        rm_idx <- unique(audit$removed_row_index[audit$sample == sn])
        rm_idx <- rm_idx[!is.na(rm_idx)]
        summary$n_split_pairs[i] <- sum(audit$sample == sn)
        summary$n_removed_pre_alignment[i] <- length(rm_idx)
        if (remove_split_peaks_pre_alignment && length(rm_idx) > 0) {
          exp_clean@Data@FactorList[[i]] <- peak_list[[i]][-rm_idx, , drop = FALSE]
        }
        summary$n_peaks_after[i] <- nrow(exp_clean@Data@FactorList[[i]])
      }
    }
    
    list(exp_clean = exp_clean, audit = audit, summary = summary)
  }
  
  get_block_time_dist <- function(block_row) {
    if (tolower(align_time_dist_mode) == "auto") {
      return(as.numeric(block_row$recommended_max_time_dist_sec))
    }
    val <- suppressWarnings(as.numeric(align_time_dist_mode))
    if (!is.finite(val) || val <= 0) {
      stop("ALIGN_TIME_DIST must be 'auto' or a positive number, got: ", align_time_dist_mode)
    }
    val
  }
  
  select_small_sample_files <- function(cdf_files, block_name) {
    if (tolower(trimws(small_sample_n)) %in% c("", "all", "none", "false")) {
      return(cdf_files)
    }
    n_keep <- suppressWarnings(as.integer(small_sample_n))
    if (!is.finite(n_keep) || n_keep < 2) {
      stop("SMALL_SAMPLE_N must be 'all' or an integer >= 2, got: ", small_sample_n)
    }
    if (length(cdf_files) <= n_keep) {
      return(cdf_files)
    }
    
    mode <- tolower(trimws(small_sample_mode))
    if (mode == "first") {
      selected <- cdf_files[seq_len(n_keep)]
    } else if (mode == "random") {
      if (!is.finite(small_sample_seed)) {
        stop("SMALL_SAMPLE_SEED must be an integer when SMALL_SAMPLE_MODE=random.")
      }
      set.seed(small_sample_seed + sum(utf8ToInt(block_name)))
      selected <- sample(cdf_files, n_keep)
    } else {
      stop("SMALL_SAMPLE_MODE must be 'first' or 'random', got: ", small_sample_mode)
    }
    selected
  }
  
  find_parent <- function(parent, x) {
    while (parent[[x]] != x) {
      parent[[x]] <- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  
  union_sets <- function(parent, a, b) {
    ra <- find_parent(parent, a)
    rb <- find_parent(parent, b)
    if (ra != rb) parent[[rb]] <- ra
    parent
  }
  
  consolidate_post_alignment_features <- function(aligned_height, aligned_area, block_name) {
    if (!post_alignment_feature_consolidation || nrow(aligned_height) < 2) {
      return(list(
        height = aligned_height,
        area = aligned_area,
        pair_audit = data.frame(),
        group_audit = data.frame()
      ))
    }
  
    height <- aligned_height
    area <- aligned_area
    height$tmean <- suppressWarnings(as.numeric(height$tmean))
    area$tmean <- suppressWarnings(as.numeric(area$tmean))
    if (!"block_feature_id" %in% names(height)) {
      height$block_feature_id <- paste(block_name, height$AlignID, sep = "__")
    }
    if (!"block_feature_id" %in% names(area)) {
      area$block_feature_id <- paste(block_name, area$AlignID, sep = "__")
    }
    height$BaseMz_tmp <- vapply(height$Spectra, base_mz_from_spectrum, numeric(1))
  
    sample_cols <- setdiff(
      names(area),
      c("AlignID", "Factor", "Spectra", "tmean", "FoundIn", "rt_block", "block_feature_id")
    )
    area_mat <- as.data.frame(lapply(area[, sample_cols, drop = FALSE], function(v) {
      suppressWarnings(as.numeric(as.character(v)))
    }), check.names = FALSE)
    area_mat[is.na(area_mat)] <- 0
  
    height_mat <- as.data.frame(lapply(height[, sample_cols, drop = FALSE], function(v) {
      suppressWarnings(as.numeric(as.character(v)))
    }), check.names = FALSE)
    height_mat[is.na(height_mat)] <- 0
  
    present <- as.matrix(area_mat > 0)
    present[is.na(present)] <- FALSE
    feature_ids <- as.character(height$block_feature_id)
  
    pairs <- list()
    ord <- order(height$tmean, height$AlignID)
    for (ii in seq_along(ord)) {
      if (ii >= length(ord)) next
      i <- ord[ii]
      later <- ord[(ii + 1):length(ord)]
      later <- later[!is.na(later)]
      if (length(later) == 0) next
      candidates <- later[(height$tmean[later] - height$tmean[i]) * 60 <= final_feature_merge_rt_sec]
      if (length(candidates) == 0) next
  
      for (j in candidates) {
        rt_diff_sec <- abs(height$tmean[j] - height$tmean[i]) * 60
        if (!is.finite(rt_diff_sec) || rt_diff_sec > final_feature_merge_rt_sec) next
        base_diff <- abs(height$BaseMz_tmp[j] - height$BaseMz_tmp[i])
        if (!is.finite(base_diff) || base_diff > final_feature_merge_base_mz_tol) next
        sp_cos <- spectral_cosine(height$Spectra[i], height$Spectra[j], mz_min = mz_min, mz_max = mz_max)
        if (!is.finite(sp_cos) || sp_cos < final_feature_merge_spectral_cos) next
  
        p1 <- present[i, , drop = TRUE]
        p2 <- present[j, , drop = TRUE]
        n_both <- sum(p1 & p2)
        n_either <- sum(p1 | p2)
        jaccard <- ifelse(n_either > 0, n_both / n_either, 0)
        merge_type <- ifelse(
          n_both == 0 || (is.finite(jaccard) && jaccard <= final_feature_alignment_split_jaccard),
          "alignment_split_like",
          "residual_or_shoulder_like"
        )
  
        pairs[[length(pairs) + 1L]] <- data.frame(
          rt_block = block_name,
          AlignID_1 = height$AlignID[i],
          AlignID_2 = height$AlignID[j],
          block_feature_id_1 = height$block_feature_id[i],
          block_feature_id_2 = height$block_feature_id[j],
          rt_1 = height$tmean[i],
          rt_2 = height$tmean[j],
          rt_diff_sec = rt_diff_sec,
          BaseMz_1 = height$BaseMz_tmp[i],
          BaseMz_2 = height$BaseMz_tmp[j],
          base_mz_diff = base_diff,
          spectral_cosine = sp_cos,
          FoundIn_1 = height$FoundIn[i],
          FoundIn_2 = height$FoundIn[j],
          n_samples_both_present = n_both,
          n_samples_either_present = n_either,
          sample_presence_jaccard = jaccard,
          merge_type = merge_type,
          auto_merge = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
  
    if (length(pairs) == 0) {
      height$BaseMz_tmp <- NULL
      return(list(
        height = height,
        area = area,
        pair_audit = data.frame(),
        group_audit = data.frame()
      ))
    }
  
    pair_audit <- bind_rows(pairs)
    merge_pair_audit <- pair_audit %>% filter(.data$auto_merge)
    if (nrow(merge_pair_audit) == 0) {
      height$BaseMz_tmp <- NULL
      return(list(
        height = height,
        area = area,
        pair_audit = pair_audit,
        group_audit = data.frame()
      ))
    }
  
    parent <- stats::setNames(as.list(feature_ids), feature_ids)
    merge_pair_audit <- merge_pair_audit[
      order(merge_pair_audit$rt_diff_sec, -merge_pair_audit$spectral_cosine),
      ,
      drop = FALSE
    ]
  
    for (k in seq_len(nrow(merge_pair_audit))) {
      fid1 <- merge_pair_audit$block_feature_id_1[k]
      fid2 <- merge_pair_audit$block_feature_id_2[k]
      root1 <- find_parent(parent, fid1)
      root2 <- find_parent(parent, fid2)
      if (identical(root1, root2)) next
  
      proposed_members <- feature_ids[vapply(feature_ids, function(fid) {
        find_parent(parent, fid) %in% c(root1, root2)
      }, logical(1))]
      proposed_idx <- match(proposed_members, feature_ids)
      rt_span_sec <- (max(height$tmean[proposed_idx], na.rm = TRUE) -
                        min(height$tmean[proposed_idx], na.rm = TRUE)) * 60
      base_mz_span <- max(height$BaseMz_tmp[proposed_idx], na.rm = TRUE) -
        min(height$BaseMz_tmp[proposed_idx], na.rm = TRUE)
  
      if (!is.finite(rt_span_sec) || rt_span_sec > final_feature_group_rt_span_sec) next
      if (!is.finite(base_mz_span) || base_mz_span > final_feature_group_base_mz_span) next
      parent <- union_sets(parent, fid1, fid2)
    }
  
    roots <- vapply(feature_ids, function(x) find_parent(parent, x), character(1))
    root_sizes <- table(roots)
    merge_roots <- names(root_sizes[root_sizes > 1])
    if (length(merge_roots) == 0) {
      height$BaseMz_tmp <- NULL
      return(list(
        height = height,
        area = area,
        pair_audit = pair_audit,
        group_audit = data.frame()
      ))
    }
  
    area_totals <- rowSums(area_mat[, sample_cols, drop = FALSE], na.rm = TRUE)
    keep_rows <- !(roots %in% merge_roots)
    new_height <- height[keep_rows, setdiff(names(height), "BaseMz_tmp"), drop = FALSE]
    new_area <- area[keep_rows, , drop = FALSE]
    group_audit <- list()
  
    for (root in merge_roots) {
      idx <- which(roots == root)
      group_present <- present[idx, , drop = FALSE]
      sample_counts <- colSums(group_present)
      n_overlap_samples <- sum(sample_counts > 1)
      n_union_samples <- sum(sample_counts > 0)
      group_jaccard <- ifelse(n_union_samples > 0, n_overlap_samples / n_union_samples, 0)
      merge_type <- ifelse(
        n_overlap_samples == 0 || group_jaccard <= final_feature_alignment_split_jaccard,
        "alignment_split_like",
        "residual_or_shoulder_like"
      )
  
      rep_idx <- idx[which.max(area_totals[idx] * pmax(suppressWarnings(as.numeric(height$FoundIn[idx])), 1))]
      rep_height <- height[rep_idx, setdiff(names(height), "BaseMz_tmp"), drop = FALSE]
      rep_area <- area[rep_idx, , drop = FALSE]
  
      group_area <- area_mat[idx, sample_cols, drop = FALSE]
      group_height <- height_mat[idx, sample_cols, drop = FALSE]
      if (nrow(group_area) == 1) {
        summed_area <- as.numeric(group_area[1, ])
        max_height <- as.numeric(group_height[1, ])
      } else {
        summed_area <- colSums(group_area, na.rm = TRUE)
        max_height <- apply(group_height, 2, max, na.rm = TRUE)
        max_height[!is.finite(max_height)] <- 0
      }
  
      rep_area[, sample_cols] <- as.list(summed_area)
      rep_height[, sample_cols] <- as.list(max_height)
      rep_area$FoundIn <- sum(summed_area > 0, na.rm = TRUE)
      rep_height$FoundIn <- rep_area$FoundIn
  
      weights <- area_totals[idx]
      if (sum(weights, na.rm = TRUE) > 0) {
        weighted_rt <- sum(height$tmean[idx] * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
        rep_area$tmean <- weighted_rt
        rep_height$tmean <- weighted_rt
      }
  
      final_id <- paste0(block_name, "__final_", length(group_audit) + 1L)
      rep_area$rt_block <- block_name
      rep_height$rt_block <- block_name
      rep_area$block_feature_id <- final_id
      rep_height$block_feature_id <- final_id
  
      new_area <- bind_rows(new_area, rep_area)
      new_height <- bind_rows(new_height, rep_height)
  
      group_audit[[length(group_audit) + 1L]] <- data.frame(
        rt_block = block_name,
        final_feature_id = final_id,
        representative_block_feature_id = final_id,
        representative_AlignID = rep_height$AlignID,
        n_erah_alignids = length(idx),
        original_AlignIDs = paste(height$AlignID[idx], collapse = ";"),
        original_block_feature_ids = paste(height$block_feature_id[idx], collapse = ";"),
        min_rt = min(height$tmean[idx], na.rm = TRUE),
        max_rt = max(height$tmean[idx], na.rm = TRUE),
        rt_span_sec = (max(height$tmean[idx], na.rm = TRUE) - min(height$tmean[idx], na.rm = TRUE)) * 60,
        BaseMz = paste(sort(unique(height$BaseMz_tmp[idx])), collapse = ";"),
        n_union_samples = n_union_samples,
        n_overlap_samples = n_overlap_samples,
        group_presence_jaccard = group_jaccard,
        merge_type = merge_type,
        total_FoundIn_before = sum(suppressWarnings(as.numeric(height$FoundIn[idx])), na.rm = TRUE),
        FoundIn_after = rep_height$FoundIn,
        stringsAsFactors = FALSE
      )
    }
  
    new_height <- new_height[order(suppressWarnings(as.numeric(new_height$tmean)), new_height$AlignID), , drop = FALSE]
    new_area <- new_area[order(suppressWarnings(as.numeric(new_area$tmean)), new_area$AlignID), , drop = FALSE]
  
    list(
      height = new_height,
      area = new_area,
      pair_audit = pair_audit,
      group_audit = bind_rows(group_audit)
    )
  }
  
  # ---- Step 3/4: deconvolve and align each block ------------------------------
  
  all_aligned_height <- list()
  all_aligned_area <- list()
  all_block_qc <- list()
  all_split_qc_summary <- list()
  all_split_qc_audit <- list()
  all_area_filter_summary <- list()
  all_area_filter_audit <- list()
  all_feature_consolidation_pair_audit <- list()
  all_feature_consolidation_group_audit <- list()
  all_alignid0_input <- list()
  all_alignid0_candidates <- list()
  all_alignid0_decisions <- list()
  all_alignid0_membership <- list()
  
  message("\n[3/5] Running block-wise deconvolution and alignment...")
  
  for (block_name in blocks) {
    block_row <- block_summary[block_summary$rt_block == block_name, , drop = FALSE]
    sample_list_file <- if (use_existing_block_folders) {
      file.path(preblock_dir, "block_folders", block_name, "sample_list.csv")
    } else {
      file.path(preblock_dir, paste0("auto_block_sample_list_", block_name, ".csv"))
    }
    if (!file.exists(sample_list_file)) {
      stop("Missing sample list for ", block_name, ": ", sample_list_file)
    }
    
    sample_list <- read.csv(sample_list_file, stringsAsFactors = FALSE, check.names = FALSE)
    cdf_files <- sample_list$cdf_file
    cdf_files <- cdf_files[file.exists(cdf_files)]
    cdf_files <- select_small_sample_files(cdf_files, block_name)
    if (length(cdf_files) < 2) {
      warning("Skipping ", block_name, ": fewer than 2 available CDF files.")
      next
    }
    
    block_dir <- file.path(out_dir, block_name)
    dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
    
    max_time_dist <- get_block_time_dist(block_row)
    
    message("\n--- ", block_name, " ---")
    message("Samples: ", length(cdf_files))
    message("max.time.dist: ", max_time_dist, " sec")
    write.csv(
      data.frame(cdf_file = cdf_files, stringsAsFactors = FALSE),
      file.path(block_dir, paste0("selected_sample_list_", block_name, ".csv")),
      row.names = FALSE
    )
    
    exp <- make_experiment(cdf_files, block_name)
    
    dec_params <- setDecPar(
      min.peak.width = min_peak_width,
      min.peak.height = min_peak_height,
      noise.threshold = noise_threshold,
      avoid.processing.mz = avoid_processing_mz,
      analysis.time = c(analysis_start, analysis_end)
    )
    
    exp_dec <- deconvolveComp(exp, decParameters = dec_params)
    save(exp_dec, file = file.path(block_dir, paste0("exp_dec_", block_name, ".rda")))
    dec_peaks <- export_deconvolved_peaks(
      exp_dec,
      block_name,
      file.path(block_dir, paste0("deconvolved_peaks_", block_name, ".csv"))
    )
    
    split_qc <- apply_split_peak_qc(exp_dec, block_name)
    area_filter_qc <- filter_peak_list_by_area_fraction(
      split_qc$exp_clean@Data@FactorList,
      area_fraction_filter,
      block_name
    )
    exp_dec_for_alignment <- split_qc$exp_clean
    exp_dec_for_alignment@Data@FactorList <- area_filter_qc$peak_list
    save(exp_dec_for_alignment, file = file.path(block_dir, paste0("exp_dec_for_alignment_", block_name, ".rda")))
    write.csv(split_qc$audit,
              file.path(block_dir, paste0("pre_alignment_split_peak_audit_", block_name, ".csv")),
              row.names = FALSE)
    write.csv(split_qc$summary,
              file.path(block_dir, paste0("pre_alignment_split_peak_summary_", block_name, ".csv")),
              row.names = FALSE)
    message("Pre-alignment split QC removed ",
            sum(split_qc$summary$n_removed_pre_alignment), " peaks.")
    write.csv(area_filter_qc$audit,
              file.path(block_dir, paste0("pre_alignment_area_filter_audit_", block_name, ".csv")),
              row.names = FALSE)
    write.csv(area_filter_qc$summary,
              file.path(block_dir, paste0("pre_alignment_area_filter_summary_", block_name, ".csv")),
              row.names = FALSE)
    message(
      "Pre-alignment within-sample area filter removed ",
      sum(area_filter_qc$summary$n_removed),
      " peaks at threshold ",
      area_fraction_filter,
      "."
    )
    
    al_params <- setAlPar(
      min.spectra.cor = min_spectra_cor,
      max.time.dist = max_time_dist,
      mz.range = mz_min:mz_max
    )
    
    # For these RT-drift blocks we expect moderate block sizes, so do not pass
    # blocks.size. eRah can fail when blocks.size exceeds the number of samples.
    exp_align_raw <- alignComp(exp_dec_for_alignment, alParameters = al_params)
    save(exp_align_raw,
         file = file.path(block_dir, paste0("exp_align_raw_before_feature_consolidation_", block_name, ".rda")))
  
    # Unified workflow: keep the raw eRah alignment object unchanged. Alignment-split-like
    # and residual/shoulder-like cases are handled together at the aligned table level.
    exp_align <- exp_align_raw
    save(exp_align, file = file.path(block_dir, paste0("exp_align_", block_name, ".rda")))
  
    aligned_height <- as.data.frame(exp_align@Results@Alignment, stringsAsFactors = FALSE)
    aligned_height$rt_block <- block_name
    aligned_height$block_feature_id <- paste(block_name, aligned_height$AlignID, sep = "__")
    
    aligned_area <- as.data.frame(alignList(exp_align, by.area = TRUE), stringsAsFactors = FALSE)
    aligned_area$rt_block <- block_name
    aligned_area$block_feature_id <- paste(block_name, aligned_area$AlignID, sep = "__")
    
    aligned_height_raw <- aligned_height
    aligned_area_raw <- aligned_area
    
    feature_consolidation_qc <- consolidate_post_alignment_features(aligned_height, aligned_area, block_name)
    aligned_height <- feature_consolidation_qc$height
    aligned_area <- feature_consolidation_qc$area
    n_removed_by_feature_consolidation <- nrow(aligned_height_raw) - nrow(aligned_height)
    write.csv(aligned_height_raw,
              file.path(block_dir, paste0("aligned_height_raw_before_feature_consolidation_", block_name, ".csv")),
              row.names = FALSE)
    write.csv(aligned_area_raw,
              file.path(block_dir, paste0("aligned_area_raw_before_feature_consolidation_", block_name, ".csv")),
              row.names = FALSE)
    write.csv(feature_consolidation_qc$pair_audit,
              file.path(block_dir, paste0("post_alignment_feature_consolidation_pair_audit_", block_name, ".csv")),
              row.names = FALSE)
    write.csv(feature_consolidation_qc$group_audit,
              file.path(block_dir, paste0("post_alignment_feature_consolidation_group_audit_", block_name, ".csv")),
              row.names = FALSE)
    message("Unified post-alignment feature consolidation removed ",
            n_removed_by_feature_consolidation, " aligned features.")

    singleton_input <- extract_alignid0_from_experiment(exp_align_raw, block_name)
    if (alignid0_rescue_enabled) {
      singleton_repair <- merge_singletons_within_block(
        aligned_height = aligned_height,
        aligned_area = aligned_area,
        singletons = singleton_input,
        block_name = block_name,
        thresholds = singleton_thresholds
      )
      aligned_height <- singleton_repair$height
      aligned_area <- singleton_repair$area
    } else {
      stop(
        "ALIGNID0_RESCUE_ENABLED=false would silently drop ", nrow(singleton_input),
        " quality-filtered AlignID=0 peaks in ", block_name,
        ". This complete-peak pipeline requires singleton rescue."
      )
    }
    write.csv(
      singleton_input,
      file.path(block_dir, paste0("alignid0_singleton_input_", block_name, ".csv")),
      row.names = FALSE
    )
    write.csv(
      singleton_repair$candidates,
      file.path(block_dir, paste0("alignid0_singleton_candidate_audit_", block_name, ".csv")),
      row.names = FALSE
    )
    write.csv(
      singleton_repair$decisions,
      file.path(block_dir, paste0("alignid0_singleton_decision_audit_", block_name, ".csv")),
      row.names = FALSE
    )
    write.csv(
      singleton_repair$membership,
      file.path(block_dir, paste0("alignid0_singleton_cluster_membership_", block_name, ".csv")),
      row.names = FALSE
    )
    message(
      "AlignID=0 complete-peak repair represented all ", nrow(singleton_input),
      " peaks; ", sum(singleton_repair$decisions$status == "merged_existing"),
      " attached to existing block features and ",
      nrow(singleton_repair$groups), " retained singleton-derived features were created."
    )
    
    write.csv(aligned_height, file.path(block_dir, paste0("aligned_height_", block_name, ".csv")), row.names = FALSE)
    write.csv(aligned_area, file.path(block_dir, paste0("aligned_area_", block_name, ".csv")), row.names = FALSE)
    
    if (do_annotation) {
      message("Annotating ", block_name, " with n.putative=", n_putative)
      exp_annot <- identifyComp(
        exp_align,
        id.database = id_database,
        mz.range = mz_min:mz_max,
        n.putative = n_putative
      )
      save(exp_annot, file = file.path(block_dir, paste0("exp_annot_", block_name, ".rda")))
      annot <- idList(exp_annot, id.database = id_database)
      annot$rt_block <- block_name
      annot$block_feature_id <- paste(block_name, annot$AlignID, sep = "__")
      write.csv(annot, file.path(block_dir, paste0("annotation_", block_name, ".csv")), row.names = FALSE)
    }
    
    sample_cols <- setdiff(
      names(aligned_height),
      c("AlignID", "Factor", "Spectra", "tmean", "FoundIn", "rt_block",
        "block_feature_id", "feature_origin", "BaseMz")
    )
    mat <- as.data.frame(lapply(aligned_height[, sample_cols, drop = FALSE], function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    }))
    found_in <- rowSums(mat > 0, na.rm = TRUE)
    
    qc <- data.frame(
      rt_block = block_name,
      n_samples = length(cdf_files),
      max_time_dist_sec = max_time_dist,
      n_deconvolved_peaks = nrow(dec_peaks),
      n_removed_split_peaks = sum(split_qc$summary$n_removed_pre_alignment),
      area_fraction_filter = area_fraction_filter,
      n_removed_area_fraction_peaks = sum(area_filter_qc$summary$n_removed),
      area_retained_pct_after_fraction_filter =
        100 * sum(area_filter_qc$summary$area_after) /
        sum(area_filter_qc$summary$area_before),
      n_post_alignment_feature_consolidation_pairs = nrow(feature_consolidation_qc$pair_audit),
      n_post_alignment_feature_consolidation_groups = nrow(feature_consolidation_qc$group_audit),
      n_removed_post_alignment_feature_consolidation_features = n_removed_by_feature_consolidation,
      n_alignment_split_like_groups = if ("merge_type" %in% names(feature_consolidation_qc$group_audit)) sum(feature_consolidation_qc$group_audit$merge_type == "alignment_split_like", na.rm = TRUE) else 0L,
      n_residual_or_shoulder_like_groups = if ("merge_type" %in% names(feature_consolidation_qc$group_audit)) sum(feature_consolidation_qc$group_audit$merge_type == "residual_or_shoulder_like", na.rm = TRUE) else 0L,
      mean_deconvolved_peaks_per_sample = nrow(dec_peaks) / length(cdf_files),
      n_aligned_features = nrow(aligned_height),
      n_alignid0_input_peaks = nrow(singleton_input),
      n_alignid0_attached_existing = sum(singleton_repair$decisions$status == "merged_existing"),
      n_singleton_derived_block_features = nrow(singleton_repair$groups),
      mean_found_in = mean(found_in),
      median_found_in = median(found_in),
      found_ge_5 = sum(found_in >= 5),
      found_ge_10 = sum(found_in >= 10),
      stringsAsFactors = FALSE
    )
    write.csv(qc, file.path(block_dir, paste0("qc_", block_name, ".csv")), row.names = FALSE)
    
    all_aligned_height[[block_name]] <- aligned_height
    all_aligned_area[[block_name]] <- aligned_area
    all_block_qc[[block_name]] <- qc
    all_split_qc_summary[[block_name]] <- split_qc$summary
    all_split_qc_audit[[block_name]] <- split_qc$audit
    all_area_filter_summary[[block_name]] <- area_filter_qc$summary
    all_area_filter_audit[[block_name]] <- area_filter_qc$audit
    all_feature_consolidation_pair_audit[[block_name]] <- feature_consolidation_qc$pair_audit
    all_feature_consolidation_group_audit[[block_name]] <- feature_consolidation_qc$group_audit
    all_alignid0_input[[block_name]] <- singleton_input
    all_alignid0_candidates[[block_name]] <- singleton_repair$candidates
    all_alignid0_decisions[[block_name]] <- singleton_repair$decisions
    all_alignid0_membership[[block_name]] <- singleton_repair$membership
  }
  
  # ---- Step 5: combined exports for strict cross-block consolidation ----------
  
  message("\n[4/5] Combining block-level outputs...")
  combined_height <- bind_rows(all_aligned_height)
  combined_area <- bind_rows(all_aligned_area)
  combined_qc <- bind_rows(all_block_qc)
  combined_split_qc_summary <- bind_rows(all_split_qc_summary)
  combined_split_qc_audit <- bind_rows(all_split_qc_audit)
  combined_area_filter_summary <- bind_rows(all_area_filter_summary)
  combined_area_filter_audit <- bind_rows(all_area_filter_audit)
  combined_feature_consolidation_pair_audit <- bind_rows(all_feature_consolidation_pair_audit)
  combined_feature_consolidation_group_audit <- bind_rows(all_feature_consolidation_group_audit)
  combined_alignid0_input <- bind_rows(all_alignid0_input)
  combined_alignid0_candidates <- bind_rows(all_alignid0_candidates)
  combined_alignid0_decisions <- bind_rows(all_alignid0_decisions)
  combined_alignid0_membership <- bind_rows(all_alignid0_membership)
  
  write.csv(combined_height, file.path(out_dir, "01_all_blocks_aligned_height_unmerged.csv"), row.names = FALSE)
  write.csv(combined_area, file.path(out_dir, "01_all_blocks_aligned_area_unmerged.csv"), row.names = FALSE)
  write.csv(combined_height, file.path(out_dir, "01_all_blocks_repaired_height_unmerged.csv"), row.names = FALSE)
  write.csv(combined_area, file.path(out_dir, "01_all_blocks_repaired_area_unmerged.csv"), row.names = FALSE)
  write.csv(combined_qc, file.path(out_dir, "01_blockwise_qc_summary.csv"), row.names = FALSE)
  write.csv(combined_split_qc_summary, file.path(out_dir, "01_pre_alignment_split_peak_summary.csv"), row.names = FALSE)
  write.csv(combined_split_qc_audit, file.path(out_dir, "01_pre_alignment_split_peak_audit.csv"), row.names = FALSE)
  write.csv(combined_area_filter_summary,
            file.path(out_dir, "01_pre_alignment_area_filter_summary.csv"),
            row.names = FALSE)
  write.csv(combined_area_filter_audit,
            file.path(out_dir, "01_pre_alignment_area_filter_audit.csv"),
            row.names = FALSE)
  write.csv(combined_feature_consolidation_pair_audit,
            file.path(out_dir, "01_post_alignment_feature_consolidation_pair_audit.csv"),
            row.names = FALSE)
  write.csv(combined_feature_consolidation_group_audit,
            file.path(out_dir, "01_post_alignment_feature_consolidation_group_audit.csv"),
            row.names = FALSE)
  write.csv(combined_alignid0_input,
            file.path(out_dir, "01_alignid0_singleton_input.csv"), row.names = FALSE)
  write.csv(combined_alignid0_candidates,
            file.path(out_dir, "01_alignid0_singleton_candidate_audit.csv"), row.names = FALSE)
  write.csv(combined_alignid0_decisions,
            file.path(out_dir, "01_alignid0_singleton_decision_audit.csv"), row.names = FALSE)
  write.csv(combined_alignid0_membership,
            file.path(out_dir, "01_alignid0_singleton_cluster_membership.csv"), row.names = FALSE)
  
  run_parameters <- data.frame(
    parameter = c(
      "cdf_dir", "requested_alignment_mode", "selected_alignment_mode",
      "single_block_threshold", "single_block_time_dist_sec",
      "n_subset", "min_block_size", "max_block_size", "target_max_time_dist",
      "small_sample_n", "small_sample_mode", "small_sample_seed",
      "min_peak_width", "min_peak_height", "area_fraction_filter",
      "noise_threshold", "analysis_time",
      "avoid_processing_mz", "min_spectra_cor", "align_time_dist_mode", "mz_range",
      "preblock_dir",
      "remove_split_peaks_pre_alignment",
      "pre_alignment_overdec_rt_sec",
      "pre_alignment_overdec_spectral_cos",
      "pre_alignment_overdec_base_mz_tol",
      "post_alignment_feature_consolidation",
      "final_feature_merge_rt_sec",
      "final_feature_merge_spectral_cos",
      "final_feature_merge_base_mz_tol",
      "final_feature_group_rt_span_sec",
      "final_feature_group_base_mz_span",
      "final_feature_alignment_split_jaccard",
      "alignid0_rescue_enabled",
      "alignid0_within_block_rt_sec",
      "alignid0_min_spectral_cosine",
      "alignid0_base_mz_tol",
      "alignid0_ambiguity_margin",
      "do_annotation", "n_putative"
    ),
    value = c(
      cdf_dir, requested_alignment_mode, selected_alignment_mode,
      single_block_threshold, single_block_time_dist_sec,
      n_subset, min_block_size, max_block_size, target_max_time_dist,
      small_sample_n, small_sample_mode, small_sample_seed,
      min_peak_width, min_peak_height, area_fraction_filter, noise_threshold,
      paste(analysis_start, analysis_end, sep = "-"),
      paste(avoid_processing_mz, collapse = ","),
      min_spectra_cor, align_time_dist_mode, paste(mz_min, mz_max, sep = "-"),
      preblock_dir,
      remove_split_peaks_pre_alignment,
      pre_alignment_overdec_rt_sec,
      pre_alignment_overdec_spectral_cos,
      pre_alignment_overdec_base_mz_tol,
      post_alignment_feature_consolidation,
      final_feature_merge_rt_sec,
      final_feature_merge_spectral_cos,
      final_feature_merge_base_mz_tol,
      final_feature_group_rt_span_sec,
      final_feature_group_base_mz_span,
      final_feature_alignment_split_jaccard,
      alignid0_rescue_enabled,
      singleton_thresholds$max_rt_diff_sec,
      singleton_thresholds$min_spectral_cosine,
      singleton_thresholds$base_mz_tol,
      singleton_thresholds$ambiguity_margin,
      do_annotation, n_putative
    )
  )
  write.csv(run_parameters, file.path(out_dir, "00_run_parameters.csv"), row.names = FALSE)
  
  message("\n[5/5] Done.")
  message("Important outputs:")
  message("  ", file.path(out_dir, "00_block_summary.csv"))
  message("  ", file.path(out_dir, "01_blockwise_qc_summary.csv"))
  message("  ", file.path(out_dir, "01_all_blocks_aligned_height_unmerged.csv"))
  message("  ", file.path(out_dir, "01_all_blocks_aligned_area_unmerged.csv"))
  message("\nNote: The combined tables are intentionally unmerged across blocks.")
  message("Use strict RT-corrected + spectral-similarity rules for cross-block consolidation.")
  
  invisible(list(
    combined_aligned_height = file.path(out_dir, "01_all_blocks_aligned_height_unmerged.csv"),
    combined_aligned_area = file.path(out_dir, "01_all_blocks_aligned_area_unmerged.csv"),
    combined_repaired_height = file.path(out_dir, "01_all_blocks_repaired_height_unmerged.csv"),
    combined_repaired_area = file.path(out_dir, "01_all_blocks_repaired_area_unmerged.csv"),
    combined_feature_metadata = file.path(out_dir, "01_all_blocks_aligned_height_unmerged.csv"),
    out_dir = out_dir
  ))
}
