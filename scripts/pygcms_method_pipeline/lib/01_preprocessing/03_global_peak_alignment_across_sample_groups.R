# =============================================================================
# Purpose: Correct group retention times and consolidate features across groups.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

# Global peak alignment across sample groups
# Self-contained stage module: no other script is loaded or executed here.

choose_adaptive_landmark_correction <- function(
    pairs, min_pairs = 8, max_identity_median_sec = 5,
    max_identity_max_sec = 15, max_linear_residual_sec = 20,
    min_slope = 0.5, max_slope = 1.5) {
  if (nrow(pairs) < min_pairs) {
    return(list(
      intercept = 0, slope = 1, n_landmarks = nrow(pairs), n_landmarks_used = 0,
      r_squared = NA_real_, median_abs_residual_sec = NA_real_,
      max_abs_residual_sec = NA_real_,
      identity_median_abs_residual_sec = NA_real_, identity_max_abs_residual_sec = NA_real_,
      linear_correction_required = FALSE,
      status = "identity_fallback_insufficient_landmarks", pairs = pairs
    ))
  }
  identity_residuals <- abs(pairs$reference_rt - pairs$block_rt) * 60
  identity_median <- median(identity_residuals, na.rm = TRUE)
  identity_max <- max(identity_residuals, na.rm = TRUE)
  if (identity_median <= max_identity_median_sec && identity_max <= max_identity_max_sec) {
    pairs$residual_sec <- identity_residuals
    return(list(
      intercept = 0, slope = 1, n_landmarks = nrow(pairs), n_landmarks_used = nrow(pairs),
      r_squared = NA_real_, median_abs_residual_sec = identity_median,
      max_abs_residual_sec = identity_max,
      identity_median_abs_residual_sec = identity_median,
      identity_max_abs_residual_sec = identity_max,
      linear_correction_required = FALSE,
      status = "adaptive_identity_residuals_acceptable", pairs = pairs
    ))
  }
  fit1 <- stats::lm(reference_rt ~ block_rt, data = pairs)
  pairs$residual_sec <- abs(stats::residuals(fit1)) * 60
  keep <- pairs$residual_sec <= max_linear_residual_sec
  pairs_used <- if (sum(keep) >= min_pairs) pairs[keep, , drop = FALSE] else pairs
  fit <- stats::lm(reference_rt ~ block_rt, data = pairs_used)
  pairs_used$residual_sec <- abs(stats::residuals(fit)) * 60
  co <- stats::coef(fit)
  intercept <- unname(co[["(Intercept)"]])
  slope <- unname(co[["block_rt"]])
  if (!is.finite(intercept) || !is.finite(slope) || slope < min_slope || slope > max_slope) {
    return(list(
      intercept = 0, slope = 1, n_landmarks = nrow(pairs), n_landmarks_used = nrow(pairs_used),
      r_squared = NA_real_, median_abs_residual_sec = identity_median,
      max_abs_residual_sec = identity_max,
      identity_median_abs_residual_sec = identity_median,
      identity_max_abs_residual_sec = identity_max,
      linear_correction_required = FALSE,
      status = "identity_fallback_invalid_linear_fit", pairs = pairs_used
    ))
  }
  list(
    intercept = intercept, slope = slope, n_landmarks = nrow(pairs),
    n_landmarks_used = nrow(pairs_used), r_squared = summary(fit)$r.squared,
    median_abs_residual_sec = median(pairs_used$residual_sec, na.rm = TRUE),
    max_abs_residual_sec = max(pairs_used$residual_sec, na.rm = TRUE),
    identity_median_abs_residual_sec = identity_median,
    identity_max_abs_residual_sec = identity_max,
    linear_correction_required = TRUE,
    status = "adaptive_landmark_linear", pairs = pairs_used
  )
}

late_anchor_weight <- function(rt_min, start_min = 30, full_min = 45) {
  if (!is.finite(start_min) || !is.finite(full_min) || full_min <= start_min) {
    stop("Late-anchor full RT must be greater than its start RT.")
  }
  x <- pmin(1, pmax(0, (rt_min - start_min) / (full_min - start_min)))
  x * x * (3 - 2 * x)
}

combine_late_anchor_shifts <- function(c24_shift_min, chrysene_shift_min,
                                       warn_sec = 15, block_sec = 30) {
  disagreement_sec <- abs(c24_shift_min - chrysene_shift_min) * 60
  shift_min <- stats::median(c(c24_shift_min, chrysene_shift_min), na.rm = TRUE)
  if (!is.finite(disagreement_sec) || !is.finite(shift_min)) {
    return(list(
      shift_min = 0, disagreement_sec = NA_real_, apply = FALSE,
      status = "blocked_missing_anchor"
    ))
  }
  if (disagreement_sec > block_sec) {
    return(list(
      shift_min = shift_min, disagreement_sec = disagreement_sec, apply = FALSE,
      status = paste0("blocked_anchor_disagreement_gt", format(block_sec, trim = TRUE), "sec")
    ))
  }
  status <- if (disagreement_sec > warn_sec) {
    paste0("warning_anchor_disagreement_", format(warn_sec, trim = TRUE), "to",
           format(block_sec, trim = TRUE), "sec")
  } else {
    paste0("accepted_anchor_agreement_le", format(warn_sec, trim = TRUE), "sec")
  }
  list(
    shift_min = shift_min, disagreement_sec = disagreement_sec,
    apply = TRUE, status = status
  )
}

parse_spectrum <- function(spectrum) {
  if (is.na(spectrum) || trimws(spectrum) == "") {
    return(data.frame(mz = numeric(0), intensity = numeric(0)))
  }
  parts <- unlist(strsplit(trimws(spectrum), "\\s+"))
  spl <- strsplit(parts, "[:,]", perl = TRUE)
  mz <- suppressWarnings(as.numeric(vapply(
    spl, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)
  )))
  intensity <- suppressWarnings(as.numeric(vapply(
    spl, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)
  )))
  keep <- !is.na(mz) & !is.na(intensity)
  data.frame(mz = round(mz[keep]), intensity = intensity[keep])
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

base_mz_from_spectrum <- function(spectrum) {
  x <- parse_spectrum(spectrum)
  if (nrow(x) == 0) return(NA_real_)
  x$mz[which.max(x$intensity)]
}

select_global_sample_columns <- function(table) {
  metadata_columns <- c(
    "AlignID", "Factor", "Spectra", "tmean", "FoundIn", "rt_block",
    "block_feature_id", "original_tmean", "feature_origin", "BaseMz",
    "corrected_tmean", "rt_shift_sec", "legacy_FoundIn",
    "legacy_median_nonzero_area"
  )
  setdiff(names(table), metadata_columns)
}

# Preserve legacy GF numbering by building ordinary eRah-derived clusters first,
# then attaching or retaining singleton-derived block features in a second pass.
stable_two_pass_global_clusters <- function(metadata, thresholds) {
  thresholds <- normalize_singleton_thresholds(thresholds)
  x <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "block_feature_id", "rt_block", "corrected_tmean", "tmean_raw",
    "BaseMz", "Spectra", "FoundIn", "median_nonzero_area"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Global clustering metadata is missing: ", paste(missing, collapse = ", "))
  if (!"TopIons" %in% names(x)) x$TopIons <- NA_character_
  if (!"feature_origin" %in% names(x)) x$feature_origin <- "erah_aligned"
  x$feature_origin[is.na(x$feature_origin) | !nzchar(x$feature_origin)] <- "erah_aligned"
  if (anyDuplicated(x$block_feature_id)) stop("Duplicate block_feature_id values in global clustering.")

  legacy_origin <- x$feature_origin %in% c("erah_aligned", "alignid0_attached")
  stable_order <- function(idx) {
    idx[order(
      as.numeric(x$corrected_tmean[idx]), as.numeric(x$BaseMz[idx]),
      as.character(x$block_feature_id[idx])
    )]
  }
  legacy_idx <- stable_order(which(legacy_origin))
  singleton_idx <- stable_order(which(!legacy_origin))

  clusters <- list()
  assignments <- vector("list", nrow(x))
  singleton_candidates <- list()
  singleton_decisions <- list()
  next_id <- 1L

  cluster_span <- function(ids, column, multiplier = 1) {
    values <- suppressWarnings(as.numeric(x[[column]][match(ids, x$block_feature_id)]))
    values <- values[is.finite(values)]
    if (length(values) < 2L) return(0)
    diff(range(values)) * multiplier
  }
  row_score <- function(row) {
    use_legacy <- as.character(row$feature_origin) %in% c("erah_aligned", "alignid0_attached") &&
      "legacy_median_nonzero_area" %in% names(row) &&
      is.finite(suppressWarnings(as.numeric(row$legacy_median_nonzero_area)))
    area_score <- if (use_legacy) {
      suppressWarnings(as.numeric(row$legacy_median_nonzero_area))
    } else {
      suppressWarnings(as.numeric(row$median_nonzero_area))
    }
    found_score <- if (use_legacy && "legacy_FoundIn" %in% names(row)) {
      suppressWarnings(as.numeric(row$legacy_FoundIn))
    } else {
      suppressWarnings(as.numeric(row$FoundIn))
    }
    if (!is.finite(area_score)) area_score <- 0
    if (!is.finite(found_score)) found_score <- 1
    area_score * max(found_score, 1)
  }
  create_cluster <- function(row) {
    gid <- sprintf("GF%05d", next_id)
    next_id <<- next_id + 1L
    clusters[[length(clusters) + 1L]] <<- list(
      GlobalFeatureID = gid,
      corrected_tmean = as.numeric(row$corrected_tmean),
      original_tmean = as.numeric(row$tmean_raw),
      BaseMz = as.numeric(row$BaseMz), Spectra = as.character(row$Spectra),
      TopIons = as.character(row$TopIons), rt_blocks = as.character(row$rt_block),
      block_feature_ids = as.character(row$block_feature_id), n_block_features = 1L,
      representative_block_feature_id = as.character(row$block_feature_id),
      representative_score = row_score(row),
      contains_legacy = as.character(row$feature_origin) %in% c("erah_aligned", "alignid0_attached")
    )
    length(clusters)
  }
  candidate_rows <- function(row) {
    out <- list()
    if (!length(clusters)) return(data.frame())
    for (cid in seq_along(clusters)) {
      cl <- clusters[[cid]]
      # Same-block ambiguities/collisions were already deliberately resolved or
      # retained during block repair and must not be undone globally.
      if (as.character(row$rt_block) %in% cl$rt_blocks) next
      rt_diff <- abs(as.numeric(row$corrected_tmean) - cl$corrected_tmean) * 60
      base_diff <- abs(as.numeric(row$BaseMz) - cl$BaseMz)
      if (!is.finite(rt_diff) || rt_diff > thresholds$max_rt_diff_sec ||
          !is.finite(base_diff) || base_diff > thresholds$base_mz_tol) next
      cosine <- spectral_cosine(as.character(row$Spectra), cl$Spectra)
      if (!is.finite(cosine) || cosine < thresholds$min_spectral_cosine) next
      proposed <- c(cl$block_feature_ids, as.character(row$block_feature_id))
      rt_span <- cluster_span(proposed, "corrected_tmean", 60)
      base_span <- cluster_span(proposed, "BaseMz", 1)
      if (!is.finite(rt_span) || rt_span > thresholds$max_rt_diff_sec ||
          !is.finite(base_span) || base_span > thresholds$base_mz_tol) next
      out[[length(out) + 1L]] <- data.frame(
        block_feature_id = as.character(row$block_feature_id),
        candidate_GlobalFeatureID = cl$GlobalFeatureID, cluster_index = cid,
        spectral_cosine = cosine, corrected_rt_diff_sec = rt_diff,
        base_mz_diff = base_diff, proposed_rt_span_sec = rt_span,
        proposed_base_mz_span = base_span,
        score = cosine - 0.05 * (rt_diff / thresholds$max_rt_diff_sec),
        stringsAsFactors = FALSE
      )
    }
    if (!length(out)) return(data.frame())
    ranked <- do.call(rbind, out)
    ranked <- ranked[order(-ranked$score, -ranked$spectral_cosine,
                           ranked$corrected_rt_diff_sec,
                           ranked$candidate_GlobalFeatureID), , drop = FALSE]
    ranked$candidate_rank <- seq_len(nrow(ranked))
    ranked
  }
  attach_to_cluster <- function(row, candidate) {
    cid <- as.integer(candidate$cluster_index)
    cl <- clusters[[cid]]
    ids <- c(cl$block_feature_ids, as.character(row$block_feature_id))
    idx <- match(ids, x$block_feature_id)
    corrected <- as.numeric(x$corrected_tmean[idx])
    original <- as.numeric(x$tmean_raw[idx])
    score <- row_score(row)
    row_is_legacy <- as.character(row$feature_origin) %in% c("erah_aligned", "alignid0_attached")
    preserve_legacy_anchor <- isTRUE(cl$contains_legacy) && !row_is_legacy
    if (!preserve_legacy_anchor && is.finite(score) && score > cl$representative_score) {
      cl$Spectra <- as.character(row$Spectra)
      cl$TopIons <- as.character(row$TopIons)
      cl$BaseMz <- as.numeric(row$BaseMz)
      cl$representative_block_feature_id <- as.character(row$block_feature_id)
      cl$representative_score <- score
    }
    if (!preserve_legacy_anchor) {
      cl$corrected_tmean <- mean(corrected, na.rm = TRUE)
      cl$original_tmean <- mean(original, na.rm = TRUE)
    }
    cl$rt_blocks <- unique(c(cl$rt_blocks, as.character(row$rt_block)))
    cl$block_feature_ids <- ids
    cl$n_block_features <- length(ids)
    cl$contains_legacy <- isTRUE(cl$contains_legacy) || row_is_legacy
    clusters[[cid]] <<- cl
    cid
  }
  assignment_row <- function(row, gid, candidate = NULL) {
    if (is.null(candidate)) {
      data.frame(
        block_feature_id = as.character(row$block_feature_id), GlobalFeatureID = gid,
        merge_spectral_cosine = NA_real_, merge_rt_diff_sec = NA_real_,
        merge_base_mz_diff = NA_real_, proposed_group_corrected_rt_span_sec = NA_real_,
        proposed_group_base_mz_span = NA_real_, stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        block_feature_id = as.character(row$block_feature_id), GlobalFeatureID = gid,
        merge_spectral_cosine = candidate$spectral_cosine,
        merge_rt_diff_sec = candidate$corrected_rt_diff_sec,
        merge_base_mz_diff = candidate$base_mz_diff,
        proposed_group_corrected_rt_span_sec = candidate$proposed_rt_span_sec,
        proposed_group_base_mz_span = candidate$proposed_base_mz_span,
        stringsAsFactors = FALSE
      )
    }
  }

  # Pass 1 intentionally reproduces the historical greedy ordering and winner rule.
  for (i in legacy_idx) {
    row <- x[i, , drop = FALSE]
    candidates <- candidate_rows(row)
    if (!nrow(candidates)) {
      cid <- create_cluster(row)
      assignments[[i]] <- assignment_row(row, clusters[[cid]]$GlobalFeatureID)
    } else {
      winner <- candidates[1, , drop = FALSE]
      cid <- attach_to_cluster(row, winner)
      assignments[[i]] <- assignment_row(row, clusters[[cid]]$GlobalFeatureID, winner)
    }
  }

  # Pass 2 permits attachment only when the best destination is clearly superior.
  for (i in singleton_idx) {
    row <- x[i, , drop = FALSE]
    candidates <- candidate_rows(row)
    if (nrow(candidates)) {
      singleton_candidates[[length(singleton_candidates) + 1L]] <- candidates
    }
    runner <- if (nrow(candidates) >= 2L) candidates$score[2] else NA_real_
    margin <- if (nrow(candidates) >= 2L) candidates$score[1] - runner else Inf
    ambiguous <- nrow(candidates) >= 2L && margin < thresholds$ambiguity_margin
    if (!nrow(candidates) || ambiguous) {
      cid <- create_cluster(row)
      gid <- clusters[[cid]]$GlobalFeatureID
      assignments[[i]] <- assignment_row(row, gid)
      status <- if (ambiguous) {
        "ambiguous_retained_global_singleton"
      } else {
        "unmatched_retained_global_singleton"
      }
      singleton_decisions[[length(singleton_decisions) + 1L]] <- data.frame(
        block_feature_id = as.character(row$block_feature_id), GlobalFeatureID = gid,
        status = status, best_score = if (nrow(candidates)) candidates$score[1] else NA_real_,
        runner_up_score = runner, score_margin = if (is.finite(margin)) margin else NA_real_,
        n_candidates = nrow(candidates), stringsAsFactors = FALSE
      )
    } else {
      winner <- candidates[1, , drop = FALSE]
      cid <- attach_to_cluster(row, winner)
      gid <- clusters[[cid]]$GlobalFeatureID
      assignments[[i]] <- assignment_row(row, gid, winner)
      singleton_decisions[[length(singleton_decisions) + 1L]] <- data.frame(
        block_feature_id = as.character(row$block_feature_id), GlobalFeatureID = gid,
        status = "merged_global_feature", best_score = winner$score,
        runner_up_score = runner, score_margin = if (is.finite(margin)) margin else NA_real_,
        n_candidates = nrow(candidates), stringsAsFactors = FALSE
      )
    }
  }
  list(
    clusters = clusters,
    assignments = if (length(assignments)) do.call(rbind, assignments) else data.frame(),
    singleton_candidates = if (length(singleton_candidates)) do.call(rbind, singleton_candidates) else data.frame(),
    singleton_decisions = if (length(singleton_decisions)) do.call(rbind, singleton_decisions) else data.frame()
  )
}



run_global_peak_alignment_across_sample_groups <- function(config = list(), within_group) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })
  
  # =============================================================================
  # Project paths
  # =============================================================================
  # Default project root. You can override before load_module():
  #   Sys.setenv(PROJECT_DIR = "/path/to/repository")
  #   Sys.setenv(BLOCKWISE_OUT_DIR = "/path/to/repository/blockwise_erah_outputs")
  
  project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
  project_dir <- normalizePath(path.expand(project_dir), winslash = "/", mustWork = FALSE)
  
  root_out_dir <- file.path(project_dir, "blockwise_erah_outputs")
  scripts_out_dir <- file.path(project_dir, "scripts", "blockwise_erah_outputs")
  
  has_blockwise_outputs <- function(path) {
    file.exists(file.path(path, "01_all_blocks_aligned_area_unmerged.csv")) &&
      file.exists(file.path(path, "01_all_blocks_aligned_height_unmerged.csv")) &&
      file.exists(file.path(path, "00_block_summary.csv"))
  }
  
  env_out_dir <- Sys.getenv("BLOCKWISE_OUT_DIR", unset = "")
  if (nzchar(trimws(env_out_dir))) {
    out_dir <- normalizePath(path.expand(env_out_dir), winslash = "/", mustWork = FALSE)
  } else if (has_blockwise_outputs(root_out_dir)) {
    out_dir <- root_out_dir
  } else if (has_blockwise_outputs(scripts_out_dir)) {
    out_dir <- scripts_out_dir
  } else {
    out_dir <- root_out_dir
  }
  
  message("PROJECT_DIR: ", project_dir)
  message("BLOCKWISE_OUT_DIR: ", out_dir)
  
  area_file <- file.path(out_dir, "01_all_blocks_aligned_area_unmerged.csv")
  height_file <- file.path(out_dir, "01_all_blocks_aligned_height_unmerged.csv")
  block_summary_file <- file.path(out_dir, "00_block_summary.csv")
  
  global_metadata_file <- file.path(out_dir, "02_global_feature_metadata.csv")
  global_area_file <- file.path(out_dir, "02_global_feature_area_matrix.csv")
  global_map_file <- file.path(out_dir, "02_block_feature_to_global_feature_map.csv")
  global_summary_file <- file.path(out_dir, "02_global_merge_summary.csv")
  cross_block_pair_audit_file <- file.path(out_dir, "02_cross_block_pair_audit.csv")
  cross_block_group_audit_file <- file.path(out_dir, "02_cross_block_group_audit.csv")
  singleton_global_candidate_audit_file <- file.path(out_dir, "02_singleton_global_candidate_audit.csv")
  singleton_global_decision_audit_file <- file.path(out_dir, "02_singleton_global_decision_audit.csv")
  
  max_rt_diff_sec <- as.numeric(Sys.getenv("GLOBAL_MERGE_RT_SEC", unset = "60"))
  min_spectral_cosine <- as.numeric(Sys.getenv("GLOBAL_MERGE_SPECTRAL_COS", unset = "0.90"))
  base_mz_tol <- as.numeric(Sys.getenv("GLOBAL_MERGE_BASE_MZ_TOL", unset = "1"))
  
  # RT correction is based on landmark feature pairs, not internal standards.
  # For each block, features with highly similar EI spectra to the reference block
  # are used to fit: reference_RT = intercept + slope * observed_RT.
  rt_correction_method <- Sys.getenv("RT_CORRECTION_METHOD", unset = "landmark_linear")
  configured_reference_block <- config$reference_block
  reference_block <- if (!is.null(configured_reference_block) &&
                         length(configured_reference_block) == 1L &&
                         !is.na(configured_reference_block) &&
                         nzchar(trimws(as.character(configured_reference_block)))) {
    as.character(configured_reference_block)
  } else {
    Sys.getenv("RT_REFERENCE_BLOCK", unset = "")
  }
  landmark_min_spectral_cosine <- as.numeric(Sys.getenv("LANDMARK_MIN_SPECTRAL_COS", unset = "0.90"))
  landmark_rt_window_sec <- as.numeric(Sys.getenv("LANDMARK_RT_WINDOW_SEC", unset = "240"))
  landmark_min_pairs <- as.integer(Sys.getenv("LANDMARK_MIN_PAIRS", unset = "8"))
  landmark_max_residual_sec <- as.numeric(Sys.getenv("LANDMARK_MAX_RESIDUAL_SEC", unset = "20"))
  landmark_min_slope <- as.numeric(Sys.getenv("LANDMARK_MIN_SLOPE", unset = "0.50"))
  landmark_max_slope <- as.numeric(Sys.getenv("LANDMARK_MAX_SLOPE", unset = "1.50"))
  landmark_identity_median_sec <- as.numeric(Sys.getenv("LANDMARK_IDENTITY_MEDIAN_SEC", unset = "5"))
  landmark_identity_max_sec <- as.numeric(Sys.getenv("LANDMARK_IDENTITY_MAX_SEC", unset = "15"))
  late_anchor_start_min <- as.numeric(Sys.getenv("ALIGNMENT_LATE_ANCHOR_START_MIN", unset = "30"))
  late_anchor_full_min <- as.numeric(Sys.getenv("ALIGNMENT_LATE_ANCHOR_FULL_MIN", unset = "45"))
  late_anchor_warn_sec <- as.numeric(Sys.getenv("ALIGNMENT_LATE_ANCHOR_WARN_SEC", unset = "15"))
  late_anchor_block_sec <- as.numeric(Sys.getenv("ALIGNMENT_LATE_ANCHOR_BLOCK_SEC", unset = "30"))
  corrected_rt_min <- as.numeric(Sys.getenv("CORRECTED_RT_MIN", unset = "3"))
  protect_corrected_rt_lower_bound <- tolower(Sys.getenv(
    "PROTECT_CORRECTED_RT_LOWER_BOUND",
    unset = "false"
  )) %in% c("1", "true", "yes")
  corrected_rt_lower_bound_mode <- Sys.getenv("CORRECTED_RT_LOWER_BOUND_MODE", unset = "block_shift")
  
  rt_correction_file <- file.path(out_dir, "02_rt_correction_by_block.csv")
  rt_landmark_file <- file.path(out_dir, "02_rt_correction_landmark_pairs.csv")
  rt_landmark_plot_dir <- file.path(out_dir, "02_rt_correction_landmark_fit_plots")
  late_anchor_file <- file.path(out_dir, "02_late_anchor_correction_by_block.csv")
  
  parse_spectrum <- function(spectrum) {
    if (is.na(spectrum) || trimws(spectrum) == "") {
      return(data.frame(mz = numeric(0), intensity = numeric(0)))
    }
    parts <- unlist(strsplit(trimws(spectrum), "\\s+"))
    spl <- strsplit(parts, "[:,]", perl = TRUE)
    mz <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)
    )))
    intensity <- suppressWarnings(as.numeric(vapply(
      spl, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)
    )))
    keep <- !is.na(mz) & !is.na(intensity)
    data.frame(mz = round(mz[keep]), intensity = intensity[keep])
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
  
  base_mz_from_spectrum <- function(spectrum) {
    x <- parse_spectrum(spectrum)
    if (nrow(x) == 0) return(NA_real_)
    x$mz[which.max(x$intensity)]
  }
  
  top_ions_from_spectrum <- function(spectrum, n = 8L) {
    x <- parse_spectrum(spectrum)
    if (nrow(x) == 0) return(NA_character_)
    x <- x[order(-x$intensity, x$mz), , drop = FALSE]
    x <- head(x, n)
    paste(paste0(x$mz, ":", x$intensity), collapse = " ")
  }
  
  median_nonzero <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    v <- v[is.finite(v) & v > 0]
    if (!length(v)) NA_real_ else median(v)
  }
  
  find_landmark_pairs <- function(block_df, ref_df) {
    pairs <- list()
    rt_window_min <- landmark_rt_window_sec / 60
    for (i in seq_len(nrow(block_df))) {
      row <- block_df[i, , drop = FALSE]
      candidates <- ref_df[
        abs(ref_df$tmean - row$tmean) <= rt_window_min &
          is.finite(ref_df$BaseMz) & is.finite(row$BaseMz) &
          abs(ref_df$BaseMz - row$BaseMz) <= base_mz_tol,
        ,
        drop = FALSE
      ]
      if (nrow(candidates) == 0) next
      cosines <- vapply(candidates$Spectra, function(s) {
        spectral_cosine(row$Spectra, s)
      }, numeric(1))
      keep <- is.finite(cosines) & cosines >= landmark_min_spectral_cosine
      if (!any(keep)) next
      candidates <- candidates[keep, , drop = FALSE]
      cosines <- cosines[keep]
      best <- which.max(cosines)
      pairs[[length(pairs) + 1L]] <- data.frame(
        rt_block = row$rt_block,
        block_feature_id = row$block_feature_id,
        block_rt = row$tmean,
        reference_block = candidates$rt_block[best],
        reference_block_feature_id = candidates$block_feature_id[best],
        reference_rt = candidates$tmean[best],
        BaseMz = row$BaseMz,
        spectral_cosine = cosines[best],
        raw_rt_diff_sec = abs(candidates$tmean[best] - row$tmean) * 60,
        stringsAsFactors = FALSE
      )
    }
    out <- bind_rows(pairs)
    if (nrow(out) == 0) return(out)
    out %>%
      arrange(.data$reference_block_feature_id, desc(.data$spectral_cosine), .data$raw_rt_diff_sec) %>%
      group_by(.data$reference_block_feature_id) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(.data$block_feature_id, desc(.data$spectral_cosine), .data$raw_rt_diff_sec) %>%
      group_by(.data$block_feature_id) %>%
      slice(1) %>%
      ungroup()
  }
  
  fit_rt_correction <- function(pairs) {
    if (rt_correction_method %in% c("adaptive_landmark", "adaptive_landmark_late_anchor")) {
      return(choose_adaptive_landmark_correction(
        pairs,
        min_pairs = landmark_min_pairs,
        max_identity_median_sec = landmark_identity_median_sec,
        max_identity_max_sec = landmark_identity_max_sec,
        max_linear_residual_sec = landmark_max_residual_sec,
        min_slope = landmark_min_slope,
        max_slope = landmark_max_slope
      ))
    }
    if (nrow(pairs) < landmark_min_pairs) {
      return(list(
        intercept = 0,
        slope = 1,
        n_landmarks = nrow(pairs),
        n_landmarks_used = 0,
        r_squared = NA_real_,
        median_abs_residual_sec = NA_real_,
        max_abs_residual_sec = NA_real_,
        status = "identity_fallback_insufficient_landmarks",
        pairs = pairs
      ))
    }
    
    fit1 <- stats::lm(reference_rt ~ block_rt, data = pairs)
    pairs$residual_sec <- abs(stats::residuals(fit1)) * 60
    keep <- pairs$residual_sec <= landmark_max_residual_sec
    if (sum(keep) >= landmark_min_pairs) {
      pairs_used <- pairs[keep, , drop = FALSE]
      fit <- stats::lm(reference_rt ~ block_rt, data = pairs_used)
    } else {
      pairs_used <- pairs
      fit <- fit1
    }
    pairs_used$residual_sec <- abs(stats::residuals(fit)) * 60
    
    co <- stats::coef(fit)
    intercept <- unname(co[["(Intercept)"]])
    slope <- unname(co[["block_rt"]])
    if (!is.finite(intercept) || !is.finite(slope) ||
        slope < landmark_min_slope || slope > landmark_max_slope) {
      pairs_used$residual_sec <- NA_real_
      return(list(
        intercept = 0,
        slope = 1,
        n_landmarks = nrow(pairs),
        n_landmarks_used = nrow(pairs_used),
        r_squared = NA_real_,
        median_abs_residual_sec = NA_real_,
        max_abs_residual_sec = NA_real_,
        status = "identity_fallback_invalid_linear_fit",
        pairs = pairs_used
      ))
    }
    
    list(
      intercept = intercept,
      slope = slope,
      n_landmarks = nrow(pairs),
      n_landmarks_used = nrow(pairs_used),
      r_squared = summary(fit)$r.squared,
      median_abs_residual_sec = median(pairs_used$residual_sec, na.rm = TRUE),
      max_abs_residual_sec = max(pairs_used$residual_sec, na.rm = TRUE),
      status = ifelse(sum(keep) >= landmark_min_pairs, "landmark_linear", "landmark_linear_no_outlier_filter"),
      pairs = pairs_used
    )
  }
  
  plot_rt_correction_landmarks <- function(landmark_table, rt_correction_table, plot_dir) {
    if (nrow(landmark_table) == 0) {
      warning("No RT correction landmark pairs available for plotting.")
      return(invisible(NULL))
    }
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    plot_df <- landmark_table %>%
      left_join(rt_correction_table, by = c("rt_block", "reference_block")) %>%
      mutate(
        fitted_reference_rt = .data$intercept + .data$slope * .data$block_rt,
        label = paste0(
          .data$rt_block,
          " vs ", .data$reference_block,
          "\nmethod = ", .data$rt_correction_method,
          "\nR2 = ", ifelse(is.na(.data$r_squared), "NA", sprintf("%.4f", .data$r_squared)),
          "\nn = ", .data$n_landmarks_used, "/", .data$n_landmarks,
          "\nmedian abs residual = ",
          ifelse(is.na(.data$median_abs_residual_sec), "NA", sprintf("%.2f s", .data$median_abs_residual_sec)),
          "\nmax abs residual = ",
          ifelse(is.na(.data$max_abs_residual_sec), "NA", sprintf("%.2f s", .data$max_abs_residual_sec))
        )
      )
    
    make_plot <- function(df) {
      ggplot(df, aes(x = .data$block_rt, y = .data$reference_rt)) +
        geom_point(aes(color = .data$spectral_cosine), size = 2.1, alpha = 0.85) +
        geom_abline(aes(intercept = .data$intercept, slope = .data$slope),
                    linewidth = 0.7, color = "#1F2937") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                    linewidth = 0.45, color = "#9CA3AF") +
        annotate("text", x = -Inf, y = Inf, label = df$label[1],
                 hjust = -0.02, vjust = 1.05, size = 3.2) +
        scale_color_viridis_c(option = "C", limits = c(landmark_min_spectral_cosine, 1),
                              name = "Spectral cosine") +
        coord_equal() +
        labs(
          x = "Observed RT in block (min)",
          y = "Reference-block RT (min)",
          title = paste0("RT correction landmarks: ", df$rt_block[1])
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
    }
    
    plot_list <- lapply(split(plot_df, plot_df$rt_block), make_plot)
    for (block_name in names(plot_list)) {
      ggsave(
        filename = file.path(plot_dir, paste0("rt_landmark_fit_", block_name, ".png")),
        plot = plot_list[[block_name]],
        width = 5.6,
        height = 5.0,
        dpi = 240
      )
    }
    
    pdf_file <- file.path(plot_dir, "rt_landmark_fit_all_blocks.pdf")
    grDevices::pdf(pdf_file, width = 5.6, height = 5.0)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (block_name in names(plot_list)) {
      print(plot_list[[block_name]])
    }
    invisible(NULL)
  }
  
  required <- c(area_file, height_file, block_summary_file)
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop("Missing required blockwise output(s):\n", paste(missing, collapse = "\n"))
  }
  
  area <- read.csv(area_file, check.names = FALSE)
  height <- read.csv(height_file, check.names = FALSE)
  block_summary <- read.csv(block_summary_file, check.names = FALSE)
  if (!"feature_origin" %in% names(area)) area$feature_origin <- "erah_aligned"
  if (!"feature_origin" %in% names(height)) height$feature_origin <- "erah_aligned"

  # Preserve the within-group RT as an audit value before landmark correction.
  if (!"original_tmean" %in% names(area)) area$original_tmean <- area$tmean
  
  if (!"block_feature_id" %in% names(area) || !"block_feature_id" %in% names(height)) {
    stop("Both area and height tables must contain block_feature_id.")
  }
  
  sample_cols <- select_global_sample_columns(area)
  
  metadata <- area %>%
    select(any_of(c(
      "block_feature_id", "rt_block", "AlignID", "Factor", "tmean", "FoundIn",
      "original_tmean", "feature_origin", "legacy_FoundIn",
      "legacy_median_nonzero_area"
    )), all_of(sample_cols)) %>%
    left_join(height %>% select(.data$block_feature_id, .data$Spectra), by = "block_feature_id") %>%
    mutate(
      tmean = suppressWarnings(as.numeric(.data$tmean)),
      original_tmean = suppressWarnings(as.numeric(.data$original_tmean)),
      tmean_raw = .data$tmean,
      FoundIn = suppressWarnings(as.numeric(.data$FoundIn)),
      BaseMz = vapply(.data$Spectra, base_mz_from_spectrum, numeric(1)),
      TopIons = vapply(.data$Spectra, top_ions_from_spectrum, character(1)),
      median_nonzero_area = apply(select(., all_of(sample_cols)), 1, median_nonzero)
    )
  
  blocks <- unique(metadata$rt_block)
  landmark_metadata <- metadata %>%
    filter(.data$feature_origin %in% c("erah_aligned", "alignid0_attached"))
  if (!nrow(landmark_metadata)) {
    stop("No ordinary eRah-aligned features are available for stable RT correction.")
  }
  if (!nzchar(reference_block)) {
    reference_block <- landmark_metadata %>%
      count(.data$rt_block, name = "n_features") %>%
      arrange(desc(.data$n_features), .data$rt_block) %>%
      slice(1) %>%
      pull(.data$rt_block)
  }
  if (!reference_block %in% blocks) {
    stop("RT_REFERENCE_BLOCK not found in data: ", reference_block)
  }
  
  ref_df <- landmark_metadata %>% filter(.data$rt_block == reference_block)
  rt_corrections <- list()
  all_landmarks <- list()
  
  for (block_name in blocks) {
    block_df <- landmark_metadata %>% filter(.data$rt_block == block_name)
    if (block_name == reference_block || rt_correction_method == "none") {
      rt_corrections[[block_name]] <- data.frame(
        rt_block = block_name,
        reference_block = reference_block,
        rt_correction_method = ifelse(block_name == reference_block, "reference_identity", "none_identity"),
        intercept = 0,
        slope = 1,
        n_landmarks = NA_integer_,
        n_landmarks_used = NA_integer_,
        r_squared = NA_real_,
        median_abs_residual_sec = NA_real_,
        max_abs_residual_sec = NA_real_,
        identity_median_abs_residual_sec = NA_real_,
        identity_max_abs_residual_sec = NA_real_,
        linear_correction_required = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    
    landmark_pairs <- find_landmark_pairs(block_df, ref_df)
    fit <- fit_rt_correction(landmark_pairs)
    rt_corrections[[block_name]] <- data.frame(
      rt_block = block_name,
      reference_block = reference_block,
      rt_correction_method = fit$status,
      intercept = fit$intercept,
      slope = fit$slope,
      n_landmarks = fit$n_landmarks,
      n_landmarks_used = fit$n_landmarks_used,
      r_squared = fit$r_squared,
      median_abs_residual_sec = fit$median_abs_residual_sec,
      max_abs_residual_sec = fit$max_abs_residual_sec,
      identity_median_abs_residual_sec = if (is.null(fit$identity_median_abs_residual_sec)) NA_real_ else fit$identity_median_abs_residual_sec,
      identity_max_abs_residual_sec = if (is.null(fit$identity_max_abs_residual_sec)) NA_real_ else fit$identity_max_abs_residual_sec,
      linear_correction_required = if (is.null(fit$linear_correction_required)) TRUE else fit$linear_correction_required,
      stringsAsFactors = FALSE
    )
    if (nrow(fit$pairs) > 0) {
      fit$pairs$used_for_fit <- fit$status != "identity_fallback_insufficient_landmarks"
      all_landmarks[[block_name]] <- fit$pairs
    }
  }
  
  rt_correction_table <- bind_rows(rt_corrections)
  landmark_table <- bind_rows(all_landmarks)

  late_anchor_table <- NULL
  if (identical(rt_correction_method, "adaptive_landmark_late_anchor")) {
    internal_standard_file <- file.path(out_dir, "00_internal_standard_rt_by_sample.csv")
    if (!file.exists(internal_standard_file)) {
      stop("Late-anchor correction requires: ", internal_standard_file)
    }
    internal_standards <- read.csv(internal_standard_file, check.names = FALSE)
    required_is_columns <- c("rt_block", "sample", "c24_rt", "chrysene_rt")
    missing_is_columns <- setdiff(required_is_columns, names(internal_standards))
    if (length(missing_is_columns)) {
      stop("Internal-standard table is missing: ", paste(missing_is_columns, collapse = ", "))
    }
    is_medians <- internal_standards %>%
      group_by(.data$rt_block) %>%
      summarise(
        n_samples = dplyr::n(),
        raw_c24_median_min = median(as.numeric(.data$c24_rt), na.rm = TRUE),
        raw_chrysene_median_min = median(as.numeric(.data$chrysene_rt), na.rm = TRUE),
        .groups = "drop"
      )
    reference_is <- is_medians %>% filter(.data$rt_block == reference_block)
    if (nrow(reference_is) != 1L) {
      stop("Reference block is absent or duplicated in internal-standard table: ", reference_block)
    }
    late_anchor_rows <- lapply(seq_len(nrow(rt_correction_table)), function(i) {
      correction <- rt_correction_table[i, , drop = FALSE]
      observed <- is_medians %>% filter(.data$rt_block == correction$rt_block)
      if (nrow(observed) != 1L) {
        combined <- combine_late_anchor_shifts(NA_real_, NA_real_, late_anchor_warn_sec, late_anchor_block_sec)
        mapped_c24 <- mapped_chrysene <- c24_shift <- chrysene_shift <- NA_real_
        n_samples <- 0L
      } else {
        mapped_c24 <- correction$intercept + correction$slope * observed$raw_c24_median_min
        mapped_chrysene <- correction$intercept + correction$slope * observed$raw_chrysene_median_min
        c24_shift <- reference_is$raw_c24_median_min - mapped_c24
        chrysene_shift <- reference_is$raw_chrysene_median_min - mapped_chrysene
        combined <- combine_late_anchor_shifts(
          c24_shift, chrysene_shift, late_anchor_warn_sec, late_anchor_block_sec
        )
        n_samples <- observed$n_samples
      }
      if (identical(as.character(correction$rt_block), reference_block)) {
        combined <- list(
          shift_min = 0, disagreement_sec = 0, apply = TRUE,
          status = "reference_block_identity"
        )
      }
      data.frame(
        rt_block = correction$rt_block,
        reference_block = reference_block,
        n_samples = n_samples,
        landmark_intercept = correction$intercept,
        landmark_slope = correction$slope,
        raw_c24_median_min = if (nrow(observed)) observed$raw_c24_median_min else NA_real_,
        raw_chrysene_median_min = if (nrow(observed)) observed$raw_chrysene_median_min else NA_real_,
        mapped_c24_before_min = mapped_c24,
        mapped_chrysene_before_min = mapped_chrysene,
        c24_residual_shift_sec = c24_shift * 60,
        chrysene_residual_shift_sec = chrysene_shift * 60,
        late_anchor_shift_sec = combined$shift_min * 60,
        anchor_disagreement_sec = combined$disagreement_sec,
        late_anchor_status = combined$status,
        late_anchor_applied = combined$apply,
        stringsAsFactors = FALSE
      )
    })
    late_anchor_table <- bind_rows(late_anchor_rows)
    write.csv(late_anchor_table, late_anchor_file, row.names = FALSE)
  }
  
  metadata <- metadata %>%
    left_join(rt_correction_table %>% select(.data$rt_block, .data$intercept, .data$slope,
                                             .data$rt_correction_method),
              by = "rt_block") %>%
    mutate(
      corrected_tmean_raw = .data$intercept + .data$slope * .data$tmean_raw
    )

  if (!is.null(late_anchor_table)) {
    metadata <- metadata %>%
      left_join(
        late_anchor_table %>% select(
          .data$rt_block, .data$late_anchor_shift_sec,
          .data$late_anchor_status, .data$late_anchor_applied
        ),
        by = "rt_block"
      ) %>%
      mutate(
        late_anchor_weight = late_anchor_weight(
          .data$corrected_tmean_raw, late_anchor_start_min, late_anchor_full_min
        ),
        late_anchor_applied_shift_sec = ifelse(
          .data$late_anchor_applied, .data$late_anchor_shift_sec, 0
        ),
        corrected_tmean_raw = .data$corrected_tmean_raw +
          .data$late_anchor_weight * .data$late_anchor_applied_shift_sec / 60
      )
  } else {
    metadata <- metadata %>% mutate(
      late_anchor_shift_sec = 0,
      late_anchor_status = "not_requested",
      late_anchor_applied = FALSE,
      late_anchor_weight = 0,
      late_anchor_applied_shift_sec = 0
    )
  }
  
  if (protect_corrected_rt_lower_bound) {
    if (tolower(corrected_rt_lower_bound_mode) == "block_shift") {
      rt_lower_offsets <- metadata %>%
        group_by(.data$rt_block) %>%
        summarise(
          corrected_rt_lower_offset_min = max(0, corrected_rt_min - min(.data$corrected_tmean_raw, na.rm = TRUE)),
          .groups = "drop"
        )
      metadata <- metadata %>%
        left_join(rt_lower_offsets, by = "rt_block") %>%
        mutate(corrected_tmean = .data$corrected_tmean_raw + .data$corrected_rt_lower_offset_min)
      rt_correction_table <- rt_correction_table %>%
        left_join(rt_lower_offsets, by = "rt_block") %>%
        mutate(corrected_rt_lower_offset_min = ifelse(
          is.na(.data$corrected_rt_lower_offset_min), 0, .data$corrected_rt_lower_offset_min
        ))
    } else if (tolower(corrected_rt_lower_bound_mode) == "floor") {
      metadata <- metadata %>%
        mutate(
          corrected_tmean = pmax(.data$corrected_tmean_raw, corrected_rt_min),
          corrected_rt_lower_offset_min = .data$corrected_tmean - .data$corrected_tmean_raw
        )
      rt_lower_offsets <- metadata %>%
        group_by(.data$rt_block) %>%
        summarise(
          corrected_rt_lower_offset_min = max(.data$corrected_rt_lower_offset_min, na.rm = TRUE),
          .groups = "drop"
        )
      rt_correction_table <- rt_correction_table %>%
        left_join(rt_lower_offsets, by = "rt_block") %>%
        mutate(corrected_rt_lower_offset_min = ifelse(
          is.na(.data$corrected_rt_lower_offset_min), 0, .data$corrected_rt_lower_offset_min
        ))
    } else {
      stop("CORRECTED_RT_LOWER_BOUND_MODE must be 'floor' or 'block_shift'.")
    }
  } else {
    metadata <- metadata %>%
      mutate(
        corrected_tmean = .data$corrected_tmean_raw,
        corrected_rt_lower_offset_min = 0
      )
    rt_correction_table$corrected_rt_lower_offset_min <- 0
  }
  
  metadata <- metadata %>%
    mutate(
      rt_correction_delta_sec = (.data$corrected_tmean - .data$tmean_raw) * 60
    )
  
  write.csv(rt_correction_table, rt_correction_file, row.names = FALSE)
  write.csv(landmark_table, rt_landmark_file, row.names = FALSE)
  plot_rt_correction_landmarks(landmark_table, rt_correction_table, rt_landmark_plot_dir)
  
  metadata <- metadata %>%
    arrange(.data$corrected_tmean, .data$BaseMz, .data$block_feature_id)
  
  cluster_rt_span_sec <- function(block_feature_ids, meta) {
    vals <- meta$corrected_tmean[meta$block_feature_id %in% block_feature_ids]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) return(0)
    (max(vals) - min(vals)) * 60
  }
  
  cluster_original_rt_span_sec <- function(block_feature_ids, meta) {
    vals <- meta$tmean_raw[meta$block_feature_id %in% block_feature_ids]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) return(0)
    (max(vals) - min(vals)) * 60
  }
  
  cluster_base_mz_span <- function(block_feature_ids, meta) {
    vals <- meta$BaseMz[meta$block_feature_id %in% block_feature_ids]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) return(0)
    max(vals) - min(vals)
  }
  
  cluster_min_pairwise_cosine <- function(block_feature_ids, meta) {
    idx <- which(meta$block_feature_id %in% block_feature_ids)
    if (length(idx) < 2) return(NA_real_)
    cmb <- utils::combn(idx, 2)
    vals <- apply(cmb, 2, function(z) {
      spectral_cosine(meta$Spectra[z[1]], meta$Spectra[z[2]])
    })
    vals <- vals[is.finite(vals)]
    if (!length(vals)) NA_real_ else min(vals)
  }
  
  global_cluster_result <- stable_two_pass_global_clusters(
    metadata,
    list(
      max_rt_diff_sec = max_rt_diff_sec,
      min_spectral_cosine = min_spectral_cosine,
      base_mz_tol = base_mz_tol,
      ambiguity_margin = as.numeric(Sys.getenv("ALIGNID0_AMBIGUITY_MARGIN", unset = "0.01"))
    )
  )
  clusters <- global_cluster_result$clusters
  assignments <- split(
    global_cluster_result$assignments,
    seq_len(nrow(global_cluster_result$assignments))
  )
  
  make_cross_block_pair_audit <- function(feature_map, meta) {
    out <- list()
    groups <- split(feature_map$block_feature_id, feature_map$GlobalFeatureID)
    for (gid in names(groups)) {
      ids <- groups[[gid]]
      if (length(ids) < 2) next
      idx <- match(ids, meta$block_feature_id)
      idx <- idx[!is.na(idx)]
      if (length(idx) < 2) next
      cmb <- utils::combn(idx, 2)
      for (k in seq_len(ncol(cmb))) {
        i1 <- cmb[1, k]
        i2 <- cmb[2, k]
        out[[length(out) + 1L]] <- data.frame(
          GlobalFeatureID = gid,
          block_feature_id_1 = meta$block_feature_id[i1],
          block_feature_id_2 = meta$block_feature_id[i2],
          rt_block_1 = meta$rt_block[i1],
          rt_block_2 = meta$rt_block[i2],
          raw_rt_1 = meta$tmean_raw[i1],
          raw_rt_2 = meta$tmean_raw[i2],
          corrected_rt_1 = meta$corrected_tmean[i1],
          corrected_rt_2 = meta$corrected_tmean[i2],
          raw_rt_diff_sec = abs(meta$tmean_raw[i1] - meta$tmean_raw[i2]) * 60,
          corrected_rt_diff_sec = abs(meta$corrected_tmean[i1] - meta$corrected_tmean[i2]) * 60,
          BaseMz_1 = meta$BaseMz[i1],
          BaseMz_2 = meta$BaseMz[i2],
          BaseMz_diff = abs(meta$BaseMz[i1] - meta$BaseMz[i2]),
          spectral_cosine = spectral_cosine(meta$Spectra[i1], meta$Spectra[i2]),
          auto_merge = TRUE,
          stringsAsFactors = FALSE
        )
      }
    }
    bind_rows(out)
  }
  
  feature_map <- bind_rows(assignments) %>%
    left_join(metadata %>% select(.data$block_feature_id, .data$rt_block, .data$AlignID,
                                  .data$original_tmean,
                                  .data$tmean_raw, .data$corrected_tmean_raw,
                                  .data$corrected_tmean, .data$rt_correction_delta_sec,
                                  .data$BaseMz, .data$FoundIn, .data$median_nonzero_area,
                                  .data$feature_origin),
              by = "block_feature_id") %>%
    arrange(.data$GlobalFeatureID, .data$rt_block, .data$AlignID)
  
  global_metadata <- bind_rows(lapply(clusters, function(cl) {
    member_meta <- metadata[metadata$block_feature_id %in% cl$block_feature_ids, , drop = FALSE]
    member_is_singleton <- !member_meta$feature_origin %in% c("erah_aligned", "alignid0_attached")
    model_residuals <- rt_correction_table$median_abs_residual_sec[
      match(unique(member_meta$rt_block), rt_correction_table$rt_block)
    ]
    model_residuals <- model_residuals[is.finite(model_residuals)]
    data.frame(
      GlobalFeatureID = cl$GlobalFeatureID,
      corrected_tmean = cl$corrected_tmean,
      mean_original_tmean = cl$original_tmean,
      min_original_tmean = min(member_meta$tmean_raw, na.rm = TRUE),
      max_original_tmean = max(member_meta$tmean_raw, na.rm = TRUE),
      original_rt_span_sec = cluster_original_rt_span_sec(cl$block_feature_ids, metadata),
      min_corrected_tmean = min(member_meta$corrected_tmean, na.rm = TRUE),
      max_corrected_tmean = max(member_meta$corrected_tmean, na.rm = TRUE),
      corrected_rt_span_sec = cluster_rt_span_sec(cl$block_feature_ids, metadata),
      BaseMz = cl$BaseMz,
      BaseMz_span = cluster_base_mz_span(cl$block_feature_ids, metadata),
      min_pairwise_spectral_cosine = cluster_min_pairwise_cosine(cl$block_feature_ids, metadata),
      TopIons = cl$TopIons,
      Spectra = cl$Spectra,
      n_blocks = length(unique(cl$rt_blocks)),
      n_block_features = cl$n_block_features,
      rt_blocks = paste(unique(cl$rt_blocks), collapse = ";"),
      block_feature_ids = paste(cl$block_feature_ids, collapse = ";"),
      representative_block_feature_id = cl$representative_block_feature_id,
      feature_origins = paste(sort(unique(member_meta$feature_origin)), collapse = ";"),
      contains_singleton_derived = any(member_is_singleton),
      n_singleton_derived_block_features = sum(member_is_singleton),
      max_block_rt_model_median_abs_residual_sec = if (length(model_residuals)) max(model_residuals) else NA_real_,
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(.data$corrected_tmean, .data$BaseMz, .data$GlobalFeatureID)
  
  cross_block_pair_audit <- make_cross_block_pair_audit(feature_map, metadata)
  cross_block_group_audit <- global_metadata %>%
    filter(.data$n_block_features > 1)
  
  area_numeric <- area %>%
    select(.data$block_feature_id, all_of(sample_cols)) %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    left_join(feature_map %>% select(.data$block_feature_id, .data$GlobalFeatureID),
              by = "block_feature_id")
  
  global_area <- area_numeric %>%
    select(.data$GlobalFeatureID, all_of(sample_cols)) %>%
    group_by(.data$GlobalFeatureID) %>%
    summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    right_join(global_metadata %>% select(.data$GlobalFeatureID), by = "GlobalFeatureID") %>%
    arrange(match(.data$GlobalFeatureID, global_metadata$GlobalFeatureID))
  
  write.csv(global_metadata, global_metadata_file, row.names = FALSE)
  write.csv(global_area, global_area_file, row.names = FALSE)
  write.csv(feature_map, global_map_file, row.names = FALSE)
  write.csv(cross_block_pair_audit, cross_block_pair_audit_file, row.names = FALSE)
  write.csv(cross_block_group_audit, cross_block_group_audit_file, row.names = FALSE)
  write.csv(global_cluster_result$singleton_candidates,
            singleton_global_candidate_audit_file, row.names = FALSE)
  write.csv(global_cluster_result$singleton_decisions,
            singleton_global_decision_audit_file, row.names = FALSE)
  
  summary_table <- data.frame(
    n_block_features = nrow(metadata),
    n_global_features = nrow(global_metadata),
    n_cross_block_merged_features = nrow(metadata) - nrow(global_metadata),
    n_cross_block_merge_groups = nrow(cross_block_group_audit),
    n_singleton_derived_block_features = sum(!metadata$feature_origin %in% c("erah_aligned", "alignid0_attached")),
    n_singleton_global_attachments = sum(global_cluster_result$singleton_decisions$status == "merged_global_feature"),
    n_singleton_global_ambiguous_retained = sum(global_cluster_result$singleton_decisions$status == "ambiguous_retained_global_singleton"),
    n_singleton_global_unmatched_retained = sum(global_cluster_result$singleton_decisions$status == "unmatched_retained_global_singleton"),
    n_samples = length(sample_cols),
    max_rt_diff_sec = max_rt_diff_sec,
    min_spectral_cosine = min_spectral_cosine,
    base_mz_tol = base_mz_tol,
    protect_corrected_rt_lower_bound = protect_corrected_rt_lower_bound,
    corrected_rt_lower_bound_mode = corrected_rt_lower_bound_mode,
    rt_correction_method = rt_correction_method,
    reference_block = reference_block,
    landmark_min_spectral_cosine = landmark_min_spectral_cosine,
    landmark_rt_window_sec = landmark_rt_window_sec,
    landmark_min_pairs = landmark_min_pairs,
    landmark_max_residual_sec = landmark_max_residual_sec,
    landmark_min_slope = landmark_min_slope,
    landmark_max_slope = landmark_max_slope,
    stringsAsFactors = FALSE
  )
  write.csv(summary_table, global_summary_file, row.names = FALSE)
  
  cat("Wrote:", global_metadata_file, "\n")
  cat("Wrote:", global_area_file, "\n")
  cat("Wrote:", global_map_file, "\n")
  cat("Wrote:", cross_block_pair_audit_file, "\n")
  cat("Wrote:", cross_block_group_audit_file, "\n")
  cat("Wrote:", singleton_global_candidate_audit_file, "\n")
  cat("Wrote:", singleton_global_decision_audit_file, "\n")
  cat("Wrote:", rt_correction_file, "\n")
  cat("Wrote:", rt_landmark_file, "\n")
  cat("Wrote:", rt_landmark_plot_dir, "\n")
  print(summary_table)
  
  invisible(list(
    global_feature_metadata = global_metadata_file,
    global_feature_area = global_area_file,
    global_feature_map = global_map_file,
    out_dir = out_dir
  ))
}
