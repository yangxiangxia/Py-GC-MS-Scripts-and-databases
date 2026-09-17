# =============================================================================
# Purpose: Provide conservative matching and merging rules for AlignID=0 peaks.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

# Shared conservative merge rules for eRah FactorList peaks with AlignID=0.
# This module deliberately has no package dependency so it can be unit tested alone.

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_singleton_thresholds <- function(thresholds) {
  defaults <- list(
    max_rt_diff_sec = 60,
    min_spectral_cosine = 0.90,
    base_mz_tol = 1,
    ambiguity_margin = 0.01
  )
  out <- defaults
  for (name in names(thresholds)) out[[name]] <- thresholds[[name]]
  out <- lapply(out, function(x) suppressWarnings(as.numeric(x)[1]))
  if (!is.finite(out$max_rt_diff_sec) || out$max_rt_diff_sec <= 0 ||
      !is.finite(out$min_spectral_cosine) || out$min_spectral_cosine < 0 ||
      out$min_spectral_cosine > 1 || !is.finite(out$base_mz_tol) ||
      out$base_mz_tol < 0 || !is.finite(out$ambiguity_margin) ||
      out$ambiguity_margin < 0) {
    stop("Invalid AlignID=0 singleton merge thresholds.", call. = FALSE)
  }
  out
}

parse_singleton_spectrum <- function(spectrum) {
  if (length(spectrum) != 1L || is.na(spectrum) || !nzchar(trimws(spectrum))) {
    return(stats::setNames(numeric(), character()))
  }
  tokens <- strsplit(trimws(spectrum), "\\s+")[[1]]
  pairs <- strsplit(tokens, "[:,]", perl = TRUE)
  mz <- suppressWarnings(as.integer(round(as.numeric(vapply(
    pairs, function(x) if (length(x)) x[[1]] else NA_character_, character(1)
  )))))
  intensity <- suppressWarnings(as.numeric(vapply(
    pairs, function(x) if (length(x) >= 2L) x[[2]] else NA_character_, character(1)
  )))
  keep <- is.finite(mz) & is.finite(intensity) & intensity >= 0
  if (!any(keep)) return(stats::setNames(numeric(), character()))
  summed <- tapply(intensity[keep], mz[keep], sum)
  stats::setNames(as.numeric(summed), names(summed))
}

singleton_base_mz <- function(spectrum) {
  parsed <- parse_singleton_spectrum(spectrum)
  if (!length(parsed)) return(NA_real_)
  as.numeric(names(parsed)[which.max(parsed)])
}

singleton_spectral_cosine <- function(spectrum_1, spectrum_2, mz_min = 46, mz_max = 650) {
  x <- parse_singleton_spectrum(spectrum_1)
  y <- parse_singleton_spectrum(spectrum_2)
  x <- x[as.numeric(names(x)) >= mz_min & as.numeric(names(x)) <= mz_max]
  y <- y[as.numeric(names(y)) >= mz_min & as.numeric(names(y)) <= mz_max]
  if (!length(x) || !length(y)) return(NA_real_)
  common <- intersect(names(x), names(y))
  numerator <- if (length(common)) sum(x[common] * y[common]) else 0
  denominator <- sqrt(sum(x * x)) * sqrt(sum(y * y))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  as.numeric(numerator / denominator)
}

empty_singleton_table <- function() {
  data.frame(
    singleton_id = character(), rt_block = character(), sample = character(),
    peak_id = character(), raw_rt = numeric(), corrected_rt = numeric(),
    area = numeric(), peak_height = numeric(), Spectra = character(),
    BaseMz = numeric(), stringsAsFactors = FALSE
  )
}

extract_alignid0_from_experiment <- function(exp_align_raw, block_name) {
  if (!isS4(exp_align_raw) || !"Data" %in% methods::slotNames(exp_align_raw)) {
    stop("exp_align_raw must be an eRah MetaboSet-like S4 object.", call. = FALSE)
  }
  factor_list <- exp_align_raw@Data@FactorList
  if (is.null(names(factor_list)) || any(!nzchar(names(factor_list)))) {
    stop("FactorList samples must be named.", call. = FALSE)
  }
  required <- c("ID", "RT", "Area", "Peak Height", "Spectra", "AlignID")
  rows <- list()
  for (sample_name in names(factor_list)) {
    peaks <- as.data.frame(factor_list[[sample_name]], stringsAsFactors = FALSE)
    missing <- setdiff(required, names(peaks))
    if (length(missing)) {
      stop("FactorList table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    align_id <- suppressWarnings(as.numeric(as.character(peaks$AlignID)))
    selected <- peaks[is.finite(align_id) & align_id == 0, , drop = FALSE]
    if (!nrow(selected)) next
    spectra <- as.character(selected$Spectra)
    peak_ids <- as.character(selected$ID)
    rows[[length(rows) + 1L]] <- data.frame(
      singleton_id = paste(block_name, sample_name, peak_ids, sep = "__"),
      rt_block = as.character(block_name), sample = as.character(sample_name),
      peak_id = peak_ids,
      raw_rt = suppressWarnings(as.numeric(as.character(selected$RT))),
      corrected_rt = suppressWarnings(as.numeric(as.character(selected$RT))),
      area = suppressWarnings(as.numeric(as.character(selected$Area))),
      peak_height = suppressWarnings(as.numeric(as.character(selected[["Peak Height"]]))),
      Spectra = spectra,
      BaseMz = vapply(spectra, singleton_base_mz, numeric(1)),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  if (!length(rows)) return(empty_singleton_table())
  out <- do.call(rbind, rows)
  out <- out[order(out$sample, out$raw_rt, out$peak_id), , drop = FALSE]
  rownames(out) <- NULL
  if (anyDuplicated(out$singleton_id)) stop("Duplicate singleton IDs extracted.", call. = FALSE)
  numeric_fields <- c("raw_rt", "corrected_rt", "area", "peak_height", "BaseMz")
  if (any(!vapply(out[numeric_fields], function(x) all(is.finite(x)), logical(1)))) {
    stop("AlignID=0 peaks contain invalid RT, area, height, or spectrum values.", call. = FALSE)
  }
  out
}

validate_peak_accounting <- function(input_peak_ids, accepted_decisions, retained_membership) {
  input <- as.character(input_peak_ids)
  if (anyNA(input) || any(!nzchar(input)) || anyDuplicated(input)) {
    stop("Input singleton IDs must be non-empty and unique.", call. = FALSE)
  }
  accepted <- if (nrow(accepted_decisions)) as.character(accepted_decisions$singleton_id) else character()
  retained <- if (nrow(retained_membership)) as.character(retained_membership$singleton_id) else character()
  represented <- c(accepted, retained)
  if (anyNA(represented) || any(!nzchar(represented)) || anyDuplicated(represented) ||
      !identical(sort(represented), sort(input))) {
    missing <- setdiff(input, represented)
    extra <- setdiff(represented, input)
    stop(
      "Singleton peak accounting failed. missing=", paste(missing, collapse = ";"),
      " extra=", paste(extra, collapse = ";"), call. = FALSE
    )
  }
  TRUE
}

rank_singleton_candidates <- function(singletons, destinations, thresholds) {
  thresholds <- normalize_singleton_thresholds(thresholds)
  singleton_required <- c("singleton_id", "Spectra", "BaseMz")
  destination_required <- c("block_feature_id", "tmean", "Spectra", "BaseMz")
  if (length(setdiff(singleton_required, names(singletons)))) {
    stop("Singleton table lacks matching fields.", call. = FALSE)
  }
  if (length(setdiff(destination_required, names(destinations)))) {
    stop("Destination table lacks matching fields.", call. = FALSE)
  }
  rt_col <- if ("corrected_rt" %in% names(singletons)) "corrected_rt" else "raw_rt"
  candidates <- list()
  decisions <- list()
  for (i in seq_len(nrow(singletons))) {
    singleton <- singletons[i, , drop = FALSE]
    rt_diff <- abs(as.numeric(destinations$tmean) - as.numeric(singleton[[rt_col]])) * 60
    base_diff <- abs(as.numeric(destinations$BaseMz) - as.numeric(singleton$BaseMz))
    eligible <- which(is.finite(rt_diff) & rt_diff <= thresholds$max_rt_diff_sec &
                        is.finite(base_diff) & base_diff <= thresholds$base_mz_tol)
    rows <- list()
    for (j in eligible) {
      cosine <- singleton_spectral_cosine(singleton$Spectra, destinations$Spectra[j])
      if (!is.finite(cosine) || cosine < thresholds$min_spectral_cosine) next
      rows[[length(rows) + 1L]] <- data.frame(
        singleton_id = as.character(singleton$singleton_id),
        block_feature_id = as.character(destinations$block_feature_id[j]),
        spectral_cosine = cosine, rt_diff_sec = rt_diff[j],
        base_mz_diff = base_diff[j],
        score = cosine - 0.05 * (rt_diff[j] / thresholds$max_rt_diff_sec),
        stringsAsFactors = FALSE
      )
    }
    ranked <- if (length(rows)) do.call(rbind, rows) else data.frame()
    if (!nrow(ranked)) {
      decisions[[length(decisions) + 1L]] <- data.frame(
        singleton_id = as.character(singleton$singleton_id),
        block_feature_id = NA_character_, status = "unmatched_retained_singleton",
        best_score = NA_real_, runner_up_score = NA_real_, score_margin = NA_real_,
        n_candidates = 0L, stringsAsFactors = FALSE
      )
      next
    }
    ranked <- ranked[order(-ranked$score, -ranked$spectral_cosine,
                           ranked$rt_diff_sec, ranked$block_feature_id), , drop = FALSE]
    ranked$candidate_rank <- seq_len(nrow(ranked))
    candidates[[length(candidates) + 1L]] <- ranked
    runner <- if (nrow(ranked) >= 2L) ranked$score[2] else NA_real_
    margin <- if (is.finite(runner)) ranked$score[1] - runner else Inf
    status <- if (nrow(ranked) >= 2L && margin < thresholds$ambiguity_margin) {
      "ambiguous_retained_singleton"
    } else {
      "matched_candidate"
    }
    decisions[[length(decisions) + 1L]] <- data.frame(
      singleton_id = as.character(singleton$singleton_id),
      block_feature_id = if (status == "matched_candidate") ranked$block_feature_id[1] else NA_character_,
      status = status, best_score = ranked$score[1], runner_up_score = runner,
      score_margin = if (is.finite(margin)) margin else NA_real_,
      n_candidates = nrow(ranked), stringsAsFactors = FALSE
    )
  }
  list(
    candidates = if (length(candidates)) do.call(rbind, candidates) else data.frame(),
    decisions = if (length(decisions)) do.call(rbind, decisions) else data.frame()
  )
}

cluster_singleton_features <- function(singletons, thresholds, id_prefix = "S") {
  thresholds <- normalize_singleton_thresholds(thresholds)
  if (!nrow(singletons)) {
    return(list(metadata = data.frame(), membership = data.frame(), area = data.frame(), height = data.frame()))
  }
  required <- c("singleton_id", "sample", "area", "peak_height", "Spectra", "BaseMz")
  if (length(setdiff(required, names(singletons)))) stop("Singleton cluster input lacks fields.", call. = FALSE)
  rt_col <- if ("corrected_rt" %in% names(singletons)) "corrected_rt" else "raw_rt"
  ord <- order(as.numeric(singletons[[rt_col]]), as.numeric(singletons$BaseMz),
               as.character(singletons$sample), as.character(singletons$singleton_id))
  x <- singletons[ord, , drop = FALSE]
  clusters <- list()
  for (i in seq_len(nrow(x))) {
    eligible <- list()
    for (k in seq_along(clusters)) {
      idx <- clusters[[k]]
      if (x$sample[i] %in% x$sample[idx]) next
      proposed <- c(idx, i)
      rt_span <- diff(range(as.numeric(x[[rt_col]][proposed]))) * 60
      base_span <- diff(range(as.numeric(x$BaseMz[proposed])))
      if (!is.finite(rt_span) || rt_span > thresholds$max_rt_diff_sec ||
          !is.finite(base_span) || base_span > thresholds$base_mz_tol) next
      cosines <- vapply(idx, function(j) singleton_spectral_cosine(x$Spectra[i], x$Spectra[j]), numeric(1))
      if (any(!is.finite(cosines)) || min(cosines) < thresholds$min_spectral_cosine) next
      rep_idx <- idx[which.max(as.numeric(x$area[idx]))]
      rt_diff <- abs(as.numeric(x[[rt_col]][i]) - as.numeric(x[[rt_col]][rep_idx])) * 60
      score <- singleton_spectral_cosine(x$Spectra[i], x$Spectra[rep_idx]) -
        0.05 * (rt_diff / thresholds$max_rt_diff_sec)
      eligible[[length(eligible) + 1L]] <- c(cluster = k, score = score)
    }
    if (!length(eligible)) {
      clusters[[length(clusters) + 1L]] <- i
      next
    }
    eligible_mat <- do.call(rbind, eligible)
    eligible_mat <- eligible_mat[order(-eligible_mat[, "score"], eligible_mat[, "cluster"]), , drop = FALSE]
    ambiguous <- nrow(eligible_mat) >= 2L &&
      (eligible_mat[1, "score"] - eligible_mat[2, "score"]) < thresholds$ambiguity_margin
    if (ambiguous) {
      clusters[[length(clusters) + 1L]] <- i
    } else {
      k <- as.integer(eligible_mat[1, "cluster"])
      clusters[[k]] <- c(clusters[[k]], i)
    }
  }
  metadata <- list()
  membership <- list()
  area_rows <- list()
  height_rows <- list()
  samples <- sort(unique(as.character(x$sample)))
  for (k in seq_along(clusters)) {
    idx <- clusters[[k]]
    feature_id <- paste0(id_prefix, sprintf("%05d", k))
    weights <- pmax(as.numeric(x$area[idx]), 0)
    tmean <- if (sum(weights) > 0) {
      sum(as.numeric(x[[rt_col]][idx]) * weights) / sum(weights)
    } else {
      mean(as.numeric(x[[rt_col]][idx]))
    }
    rep_idx <- idx[which.max(as.numeric(x$area[idx]))]
    origin <- if (length(idx) == 1L) "singleton_unmatched" else "singleton_cluster"
    metadata[[k]] <- data.frame(
      block_feature_id = feature_id, tmean = tmean, FoundIn = length(unique(x$sample[idx])),
      Spectra = as.character(x$Spectra[rep_idx]), BaseMz = as.numeric(x$BaseMz[rep_idx]),
      feature_origin = origin, stringsAsFactors = FALSE
    )
    membership[[k]] <- data.frame(
      block_feature_id = feature_id, singleton_id = as.character(x$singleton_id[idx]),
      sample = as.character(x$sample[idx]), rt_block = as.character(x$rt_block[idx] %||% NA_character_),
      stringsAsFactors = FALSE
    )
    area_values <- stats::setNames(rep(0, length(samples)), samples)
    height_values <- stats::setNames(rep(0, length(samples)), samples)
    for (j in idx) {
      area_values[[as.character(x$sample[j])]] <- as.numeric(x$area[j])
      height_values[[as.character(x$sample[j])]] <- as.numeric(x$peak_height[j])
    }
    area_rows[[k]] <- data.frame(block_feature_id = feature_id, as.list(area_values), check.names = FALSE)
    height_rows[[k]] <- data.frame(block_feature_id = feature_id, as.list(height_values), check.names = FALSE)
  }
  list(
    metadata = do.call(rbind, metadata), membership = do.call(rbind, membership),
    area = do.call(rbind, area_rows), height = do.call(rbind, height_rows)
  )
}

singleton_sample_columns <- function(table) {
  metadata <- c(
    "AlignID", "Factor", "Spectra", "tmean", "FoundIn", "rt_block",
    "block_feature_id", "feature_origin", "BaseMz", "original_tmean",
    "corrected_tmean", "rt_shift_sec", "legacy_FoundIn",
    "legacy_median_nonzero_area"
  )
  setdiff(names(table), metadata)
}

make_singleton_feature_rows <- function(template, clustered, block_name, value_table, kind) {
  if (!nrow(clustered$metadata)) return(template[0, , drop = FALSE])
  rows <- vector("list", nrow(clustered$metadata))
  sample_cols <- singleton_sample_columns(template)
  for (i in seq_len(nrow(clustered$metadata))) {
    row <- if (nrow(template)) {
      template[1L, , drop = FALSE]
    } else {
      as.data.frame(
        stats::setNames(lapply(names(template), function(name) NA), names(template)),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
    for (name in names(row)) row[[name]] <- NA
    row[1, intersect(sample_cols, names(row))] <- 0
    meta <- clustered$metadata[i, , drop = FALSE]
    if ("AlignID" %in% names(row)) row$AlignID <- -i
    if ("Factor" %in% names(row)) row$Factor <- paste(
      clustered$membership$singleton_id[clustered$membership$block_feature_id == meta$block_feature_id],
      collapse = ";"
    )
    for (name in intersect(c("Spectra", "tmean", "FoundIn", "BaseMz", "feature_origin"), names(row))) {
      row[[name]] <- meta[[name]]
    }
    if ("rt_block" %in% names(row)) row$rt_block <- block_name
    if ("block_feature_id" %in% names(row)) row$block_feature_id <- meta$block_feature_id
    values <- value_table[value_table$block_feature_id == meta$block_feature_id, , drop = FALSE]
    for (sample_name in intersect(sample_cols, names(values))) row[[sample_name]] <- values[[sample_name]][1]
    rows[[i]] <- row
  }
  do.call(rbind, rows)
}

merge_singletons_within_block <- function(aligned_height, aligned_area, singletons,
                                          block_name, thresholds) {
  thresholds <- normalize_singleton_thresholds(thresholds)
  height <- as.data.frame(aligned_height, stringsAsFactors = FALSE, check.names = FALSE)
  area <- as.data.frame(aligned_area, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(height) != nrow(area) || !identical(as.character(height$block_feature_id),
                                                as.character(area$block_feature_id))) {
    stop("Aligned height and area tables are not row-aligned.", call. = FALSE)
  }
  if (!"feature_origin" %in% names(height)) height$feature_origin <- rep("erah_aligned", nrow(height))
  if (!"feature_origin" %in% names(area)) area$feature_origin <- rep("erah_aligned", nrow(area))
  if (!"BaseMz" %in% names(height)) height$BaseMz <- vapply(height$Spectra, singleton_base_mz, numeric(1))
  if (!"BaseMz" %in% names(area)) {
    area$BaseMz <- height$BaseMz[match(area$block_feature_id, height$block_feature_id)]
  }
  pre_repair_sample_cols <- singleton_sample_columns(area)
  pre_repair_area <- as.data.frame(lapply(area[, pre_repair_sample_cols, drop = FALSE], function(values) {
    values <- suppressWarnings(as.numeric(as.character(values)))
    values[!is.finite(values)] <- 0
    values
  }), check.names = FALSE)
  if (!"legacy_FoundIn" %in% names(area)) {
    area$legacy_FoundIn <- rowSums(pre_repair_area > 0)
    height$legacy_FoundIn <- area$legacy_FoundIn
  }
  if (!"legacy_median_nonzero_area" %in% names(area)) {
    area$legacy_median_nonzero_area <- if (nrow(pre_repair_area)) {
      apply(pre_repair_area, 1, function(values) {
        values <- values[values > 0]
        if (length(values)) median(values) else NA_real_
      })
    } else numeric()
    height$legacy_median_nonzero_area <- area$legacy_median_nonzero_area
  }
  if (!nrow(singletons)) {
    return(list(
      height = height, area = area, remaining_singletons = singletons,
      candidates = data.frame(), decisions = data.frame(), groups = data.frame(),
      membership = data.frame(), accepted_membership = data.frame(),
      retained_membership = data.frame()
    ))
  }
  ranked <- rank_singleton_candidates(singletons, height, thresholds)
  decisions <- ranked$decisions
  decisions$destination_block_feature_id <- decisions$block_feature_id
  sample_cols <- union(singleton_sample_columns(area), as.character(singletons$sample))
  for (sample_name in setdiff(sample_cols, names(area))) {
    area[[sample_name]] <- 0
    height[[sample_name]] <- 0
  }
  for (sample_name in sample_cols) {
    area[[sample_name]] <- suppressWarnings(as.numeric(as.character(area[[sample_name]])))
    height[[sample_name]] <- suppressWarnings(as.numeric(as.character(height[[sample_name]])))
    area[[sample_name]][!is.finite(area[[sample_name]])] <- 0
    height[[sample_name]][!is.finite(height[[sample_name]])] <- 0
  }
  for (i in seq_len(nrow(singletons))) {
    decision_idx <- match(singletons$singleton_id[i], decisions$singleton_id)
    if (decisions$status[decision_idx] != "matched_candidate") next
    feature_idx <- match(decisions$block_feature_id[decision_idx], area$block_feature_id)
    sample_name <- as.character(singletons$sample[i])
    if (area[[sample_name]][feature_idx] > 0 || height[[sample_name]][feature_idx] > 0) {
      decisions$status[decision_idx] <- "collision_retained_singleton"
      decisions$block_feature_id[decision_idx] <- NA_character_
      next
    }
    area[[sample_name]][feature_idx] <- as.numeric(singletons$area[i])
    height[[sample_name]][feature_idx] <- as.numeric(singletons$peak_height[i])
    decisions$status[decision_idx] <- "merged_existing"
    height$feature_origin[feature_idx] <- "alignid0_attached"
    area$feature_origin[feature_idx] <- "alignid0_attached"
  }
  retained_ids <- decisions$singleton_id[decisions$status != "merged_existing"]
  retained <- singletons[match(retained_ids, singletons$singleton_id), , drop = FALSE]
  clustered <- cluster_singleton_features(retained, thresholds, paste0(block_name, "__S"))
  if (nrow(clustered$metadata)) {
    height_rows <- make_singleton_feature_rows(height, clustered, block_name, clustered$height, "height")
    area_rows <- make_singleton_feature_rows(area, clustered, block_name, clustered$area, "area")
    height <- rbind(height, height_rows)
    area <- rbind(area, area_rows)
  }
  area_samples <- singleton_sample_columns(area)
  for (table_name in c("height", "area")) {
    table <- get(table_name)
    numeric_matrix <- as.data.frame(lapply(table[, area_samples, drop = FALSE], function(x) {
      values <- suppressWarnings(as.numeric(as.character(x)))
      values[!is.finite(values)] <- 0
      values
    }), check.names = FALSE)
    table$FoundIn <- rowSums(numeric_matrix > 0)
    assign(table_name, table)
  }
  accepted <- decisions[decisions$status == "merged_existing", , drop = FALSE]
  retained_membership <- clustered$membership
  validate_peak_accounting(singletons$singleton_id, accepted, retained_membership)
  list(
    height = height, area = area, remaining_singletons = retained,
    candidates = ranked$candidates, decisions = decisions,
    groups = clustered$metadata, membership = clustered$membership,
    accepted_membership = accepted[, c("singleton_id", "destination_block_feature_id"), drop = FALSE],
    retained_membership = retained_membership
  )
}
