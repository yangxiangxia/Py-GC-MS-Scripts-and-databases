# =============================================================================
# Purpose: Group CDF samples using internal-standard retention-time variation.
# Associated manuscript:
# "A reproducible workflow for high-throughput pyrolysis gas
# chromatography–mass spectrometry analysis of soil organic matter".
# Usage, inputs, outputs, and parameters: see README.md.
# =============================================================================

# Internal-standard-guided sample grouping
# Self-contained module: it does not source or execute another script.

discover_alignment_cdf_files <- function(cdf_dir, minimum = 2L) {
  if (!dir.exists(cdf_dir)) {
    stop("CDF directory not found: ", cdf_dir, call. = FALSE)
  }
  entries <- list.files(cdf_dir, full.names = TRUE, all.files = FALSE)
  entries <- entries[file.info(entries)$isdir %in% FALSE]
  cdf_files <- entries[grepl("\\.(cdf|netcdf)$", entries, ignore.case = TRUE)]
  cdf_files <- sort(normalizePath(cdf_files, winslash = "/", mustWork = TRUE))
  minimum <- as.integer(minimum)
  if (length(cdf_files) < minimum) {
    stop(
      "Alignment route requires at least ", minimum,
      " CDF/NetCDF files; found ", length(cdf_files), " in ", cdf_dir,
      call. = FALSE
    )
  }
  cdf_files
}

select_alignment_mode <- function(requested_mode, n_samples, threshold) {
  mode <- tolower(trimws(as.character(requested_mode)))
  if (!mode %in% c("auto", "single", "blockwise")) {
    stop("Unsupported ALIGNMENT_MODE: ", requested_mode, call. = FALSE)
  }
  if (mode != "auto") return(mode)

  threshold <- as.integer(threshold)
  if (is.na(threshold) || threshold < 2L) {
    stop("ALIGNMENT_SINGLE_BLOCK_THRESHOLD must be at least 2", call. = FALSE)
  }
  if (as.integer(n_samples) <= threshold) "single" else "blockwise"
}

write_single_block_preblock <- function(preblock_dir, cdf_files, align_time_dist_sec) {
  dir.create(preblock_dir, recursive = TRUE, showWarnings = FALSE)
  cdf_files <- sort(normalizePath(cdf_files, winslash = "/", mustWork = TRUE))
  sample_names <- tools::file_path_sans_ext(basename(cdf_files))
  block_name <- "block1"

  summary_path <- file.path(preblock_dir, "auto_block_summary.csv")
  rt_path <- file.path(preblock_dir, "auto_block_raw_internal_standard_rt_by_sample.csv")
  sample_path <- file.path(preblock_dir, paste0("auto_block_sample_list_", block_name, ".csv"))

  write.csv(data.frame(
    rt_block = block_name,
    n_samples = length(cdf_files),
    max_internal_standard_range_sec = 0,
    recommended_max_time_dist_sec = as.numeric(align_time_dist_sec),
    stringsAsFactors = FALSE
  ), summary_path, row.names = FALSE)

  write.csv(data.frame(
    sample = sample_names,
    cdf_file = cdf_files,
    rt_block = block_name,
    stringsAsFactors = FALSE
  ), rt_path, row.names = FALSE)

  write.csv(data.frame(
    sample = sample_names,
    cdf_file = cdf_files,
    stringsAsFactors = FALSE
  ), sample_path, row.names = FALSE)

  invisible(list(summary = summary_path, retention_times = rt_path, samples = sample_path))
}


build_internal_standard_sample_groups <- function() {
  suppressPackageStartupMessages({
    library(ncdf4)
    library(dplyr)
    library(ggplot2)
  })
  
  # =============================================================================
  # Pre-analysis RT-drift block assignment for GC-MS CDF samples
  #
  # Purpose:
  #   Before formal eRah processing, extract internal-standard EIC apex RTs from
  #   raw CDF files, then split samples into blocks that minimize within-block
  #   internal-standard RT distance. Each block can then be processed separately
  #   with a smaller max.time.dist to reduce alignment false positives.
  #
  # Internal standards used here:
  #   N-Tetracosane-d50: m/z 66, expected RT ~50.7 min
  #   Chrysene-d12:     m/z 240, expected RT ~52.7 min
  #
  # Main environment variables:
  #   N_SUBSET=50 or all
  #   MIN_BLOCK_SIZE=8
  #   MAX_BLOCK_SIZE=20
  #   TARGET_MAX_TIME_DIST=30
  #   CDF_DIR=/path/to/gcms_cdf_data
  #   COPY_BLOCK_FILES=true   # true copies CDF files; false creates symlinks
  # =============================================================================
  
  cdf_dir <- Sys.getenv(
    "CDF_DIR",
    unset = file.path(Sys.getenv("PROJECT_DIR", unset = getwd()), "data", "raw")
  )
  out_dir <- Sys.getenv(
    "PREBLOCK_DIR",
    unset = file.path(getwd(), "cdf_qc", "pre_alignment_blocks")
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  block_folder_root <- file.path(out_dir, "block_folders")
  dir.create(block_folder_root, recursive = TRUE, showWarnings = FALSE)
  
  all_files <- sort(list.files(cdf_dir, pattern = "\\.cdf$", full.names = TRUE))
  if (length(all_files) == 0) stop("No CDF files found in: ", cdf_dir)
  
  n_subset_env <- Sys.getenv("N_SUBSET", unset = "all")
  n_subset <- if (tolower(n_subset_env) == "all") length(all_files) else as.integer(n_subset_env)
  if (is.na(n_subset) || n_subset < 1) n_subset <- length(all_files)
  n_subset <- min(n_subset, length(all_files))
  subset_idx <- unique(round(seq(1, length(all_files), length.out = n_subset)))
  cdf_files <- all_files[subset_idx]
  
  min_block_size <- as.integer(Sys.getenv("MIN_BLOCK_SIZE", unset = "8"))
  max_block_size <- as.integer(Sys.getenv("MAX_BLOCK_SIZE", unset = "20"))
  target_max_time_dist <- as.numeric(Sys.getenv("TARGET_MAX_TIME_DIST", unset = "30"))
  copy_block_files <- tolower(Sys.getenv("COPY_BLOCK_FILES", unset = "true")) %in% c("1", "true", "yes")
  if (is.na(min_block_size) || min_block_size < 2) min_block_size <- 8
  if (is.na(max_block_size) || max_block_size < min_block_size) max_block_size <- max(min_block_size, 20)
  if (is.na(target_max_time_dist) || target_max_time_dist <= 0) target_max_time_dist <- 30
  
  moving_average <- function(x, n = 5) {
    y <- stats::filter(x, rep(1 / n, n), sides = 2)
    y[is.na(y)] <- x[is.na(y)]
    as.numeric(y)
  }
  
  read_eic <- function(file, target_mz, mz_tol = 0.5) {
    nc <- nc_open(file)
    on.exit(nc_close(nc), add = TRUE)
  
    rt <- ncvar_get(nc, "scan_acquisition_time") / 60
    scan_index <- ncvar_get(nc, "scan_index")
    point_count <- ncvar_get(nc, "point_count")
    mz_all <- ncvar_get(nc, "mass_values")
    int_all <- ncvar_get(nc, "intensity_values")
  
    eic <- numeric(length(rt))
    for (i in seq_along(rt)) {
      n <- point_count[i]
      if (is.na(n) || n <= 0) next
      start <- scan_index[i] + 1
      idx <- start:(start + n - 1)
      keep <- abs(mz_all[idx] - target_mz) <= mz_tol
      if (any(keep)) eic[i] <- sum(int_all[idx][keep], na.rm = TRUE)
    }
    data.frame(rt = rt, intensity = eic)
  }
  
  find_eic_apex <- function(file, target_mz, expected_rt, rt_window = 5,
                            mz_tol = 0.5, smooth_n = 5, min_snr = 5) {
    eic <- read_eic(file, target_mz = target_mz, mz_tol = mz_tol)
    x <- eic %>%
      filter(.data$rt >= expected_rt - rt_window,
             .data$rt <= expected_rt + rt_window) %>%
      arrange(.data$rt)
  
    if (nrow(x) < 10) {
      return(data.frame(found = FALSE, rt = NA_real_, height = NA_real_,
                        baseline = NA_real_, snr = NA_real_))
    }
  
    x$y <- moving_average(x$intensity, smooth_n)
    baseline <- as.numeric(stats::quantile(x$y, 0.10, na.rm = TRUE))
    noise <- stats::mad(diff(x$y), na.rm = TRUE)
    if (!is.finite(noise) || noise <= 0) noise <- stats::median(abs(diff(x$y)), na.rm = TRUE)
    if (!is.finite(noise) || noise <= 0) noise <- 1
  
    apex_idx <- which.max(x$y)
    height <- x$y[apex_idx] - baseline
    snr <- height / noise
    found <- is.finite(height) && height > 0 && snr >= min_snr
  
    data.frame(
      found = found,
      rt = if (found) x$rt[apex_idx] else NA_real_,
      height = if (found) height else NA_real_,
      baseline = baseline,
      snr = snr
    )
  }
  
  extract_internal_standard_rts <- function(files) {
    message("Extracting raw internal-standard EIC RTs for ", length(files), " samples...")
    bind_rows(lapply(files, function(file) {
      c24 <- find_eic_apex(file, target_mz = 66, expected_rt = 50.7, rt_window = 5)
      chrys <- find_eic_apex(file, target_mz = 240, expected_rt = 52.7, rt_window = 5)
      data.frame(
        sample = tools::file_path_sans_ext(basename(file)),
        cdf_file = file,
        c24_found = c24$found,
        c24_rt = c24$rt,
        c24_height = c24$height,
        c24_snr = c24$snr,
        chrysene_found = chrys$found,
        chrysene_rt = chrys$rt,
        chrysene_height = chrys$height,
        chrysene_snr = chrys$snr,
        stringsAsFactors = FALSE
      )
    }))
  }
  
  block_cost <- function(df, i, j) {
    x <- df[i:j, , drop = FALSE]
    max(
      diff(range(x$c24_drift_sec, na.rm = TRUE)),
      diff(range(x$chrysene_drift_sec, na.rm = TRUE))
    )
  }
  
  partition_minimax <- function(df, min_size, max_size) {
    n <- nrow(df)
    if (n == 0) return(list(blocks = list(), objective = NA_real_, total_cost = NA_real_))
    if (n < min_size) {
      return(list(blocks = list(seq_len(n)), objective = block_cost(df, 1, n), total_cost = block_cost(df, 1, n)))
    }
  
    dp_max <- rep(Inf, n + 1)
    dp_sum <- rep(Inf, n + 1)
    prev <- rep(NA_integer_, n + 1)
    dp_max[1] <- 0
    dp_sum[1] <- 0
  
    for (end in seq_len(n)) {
      for (size in min_size:max_size) {
        start <- end - size + 1
        if (start < 1) next
        prev_idx <- start
        if (!is.finite(dp_max[prev_idx])) next
  
        cst <- block_cost(df, start, end)
        cand_max <- max(dp_max[prev_idx], cst)
        cand_sum <- dp_sum[prev_idx] + cst
  
        # Lexicographic optimization:
        #   1) minimize the worst within-block RT distance
        #   2) if tied, minimize the total within-block RT distance
        if (cand_max < dp_max[end + 1] ||
            (isTRUE(all.equal(cand_max, dp_max[end + 1])) && cand_sum < dp_sum[end + 1])) {
          dp_max[end + 1] <- cand_max
          dp_sum[end + 1] <- cand_sum
          prev[end + 1] <- start
        }
      }
    }
  
    if (!is.finite(dp_max[n + 1])) {
      warning("No feasible partition with the requested block sizes. Falling back to one block.")
      return(list(blocks = list(seq_len(n)), objective = block_cost(df, 1, n), total_cost = block_cost(df, 1, n)))
    }
  
    blocks <- list()
    pos <- n + 1
    while (pos > 1) {
      start <- prev[pos]
      end <- pos - 1
      blocks[[length(blocks) + 1]] <- start:end
      pos <- start
    }
    blocks <- rev(blocks)
    list(blocks = blocks, objective = dp_max[n + 1], total_cost = dp_sum[n + 1])
  }
  
  rt_df <- extract_internal_standard_rts(cdf_files) %>%
    mutate(
      c24_drift_sec = (.data$c24_rt - median(.data$c24_rt, na.rm = TRUE)) * 60,
      chrysene_drift_sec = (.data$chrysene_rt - median(.data$chrysene_rt, na.rm = TRUE)) * 60,
      mean_drift_sec = rowMeans(cbind(.data$c24_drift_sec, .data$chrysene_drift_sec), na.rm = TRUE)
    )
  
  complete <- rt_df %>%
    filter(is.finite(.data$c24_drift_sec), is.finite(.data$chrysene_drift_sec))
  
  if (nrow(complete) > 1) {
    pca <- prcomp(complete[, c("c24_drift_sec", "chrysene_drift_sec")], center = TRUE, scale. = TRUE)
    complete$rt_order_score <- pca$x[, 1]
  } else {
    complete$rt_order_score <- complete$mean_drift_sec
  }
  
  complete <- complete %>% arrange(.data$rt_order_score)
  partition <- partition_minimax(complete, min_block_size, max_block_size)
  
  complete$rt_block <- NA_character_
  for (b in seq_along(partition$blocks)) {
    complete$rt_block[partition$blocks[[b]]] <- paste0("block", b)
  }
  
  rt_df <- rt_df %>%
    left_join(complete[, c("sample", "rt_order_score", "rt_block")], by = "sample")
  rt_df$rt_block[is.na(rt_df$rt_block)] <- "needs_manual_check"
  rt_df$rt_block_num <- suppressWarnings(as.integer(sub("^block", "", rt_df$rt_block)))
  
  block_summary <- rt_df %>%
    filter(.data$rt_block != "needs_manual_check") %>%
    group_by(.data$rt_block, .data$rt_block_num) %>%
    summarise(
      n_samples = n(),
      c24_rt_min = min(.data$c24_rt, na.rm = TRUE),
      c24_rt_max = max(.data$c24_rt, na.rm = TRUE),
      c24_rt_median = median(.data$c24_rt, na.rm = TRUE),
      c24_rt_range_sec = diff(range(.data$c24_rt, na.rm = TRUE)) * 60,
      chrysene_rt_min = min(.data$chrysene_rt, na.rm = TRUE),
      chrysene_rt_max = max(.data$chrysene_rt, na.rm = TRUE),
      chrysene_rt_median = median(.data$chrysene_rt, na.rm = TRUE),
      chrysene_rt_range_sec = diff(range(.data$chrysene_rt, na.rm = TRUE)) * 60,
      max_internal_standard_range_sec = pmax(.data$c24_rt_range_sec, .data$chrysene_rt_range_sec),
      recommended_max_time_dist_sec = case_when(
        .data$max_internal_standard_range_sec <= target_max_time_dist ~ target_max_time_dist,
        .data$max_internal_standard_range_sec <= 45 ~ 45,
        .data$max_internal_standard_range_sec <= 60 ~ 60,
        .data$max_internal_standard_range_sec <= 90 ~ 90,
        TRUE ~ 120
      ),
      mean_drift_min = min(.data$mean_drift_sec, na.rm = TRUE),
      mean_drift_max = max(.data$mean_drift_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data$rt_block_num)
  
  manual_summary <- rt_df %>%
    filter(.data$rt_block == "needs_manual_check") %>%
    summarise(n_manual_check = n())
  
  method_summary <- data.frame(
    n_samples_input = length(cdf_files),
    n_samples_complete_internal_standards = nrow(complete),
    min_block_size = min_block_size,
    max_block_size = max_block_size,
    target_max_time_dist_sec = target_max_time_dist,
    achieved_worst_block_range_sec = partition$objective,
    n_blocks = length(partition$blocks),
    stringsAsFactors = FALSE
  )
  
  write.csv(rt_df, file.path(out_dir, "auto_block_raw_internal_standard_rt_by_sample.csv"), row.names = FALSE)
  write.csv(block_summary, file.path(out_dir, "auto_block_summary.csv"), row.names = FALSE)
  write.csv(manual_summary, file.path(out_dir, "auto_block_manual_check_summary.csv"), row.names = FALSE)
  write.csv(method_summary, file.path(out_dir, "auto_block_method_summary.csv"), row.names = FALSE)
  
  sample_assignment <- rt_df %>%
    arrange(is.na(.data$rt_block_num), .data$rt_block_num, .data$mean_drift_sec, .data$sample) %>%
    select(
      "sample", "rt_block", "cdf_file",
      "c24_rt", "chrysene_rt",
      "c24_drift_sec", "chrysene_drift_sec", "mean_drift_sec",
      "c24_found", "chrysene_found"
    )
  write.csv(sample_assignment, file.path(out_dir, "sample_block_assignment.csv"), row.names = FALSE)
  
  for (b in unique(sample_assignment$rt_block)) {
    block_files <- rt_df %>%
      filter(.data$rt_block == b) %>%
      arrange(.data$mean_drift_sec, .data$sample) %>%
      select("sample", "rt_block", "cdf_file")
    safe_b <- gsub("[^A-Za-z0-9_]+", "_", b)
    write.csv(block_files, file.path(out_dir, paste0("auto_block_sample_list_", safe_b, ".csv")), row.names = FALSE)
  
    block_dir <- file.path(block_folder_root, safe_b)
    dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(block_files, file.path(block_dir, "sample_list.csv"), row.names = FALSE)
  
    for (src in block_files$cdf_file) {
      dest <- file.path(block_dir, basename(src))
      if (copy_block_files) {
        if (file.exists(dest) || file.info(dest)$isdir %in% FALSE) {
          unlink(dest)
        }
        ok <- file.copy(src, dest, overwrite = FALSE)
        if (!isTRUE(ok)) {
          warning("Could not copy ", src, " to ", dest, "; leaving it in sample_list.csv only.")
        }
      } else {
        if (file.exists(dest)) next
        ok <- file.symlink(src, dest)
        if (!isTRUE(ok)) {
          warning("Could not create symlink for ", src, "; leaving it in sample_list.csv only.")
        }
      }
    }
  }
  
  p <- ggplot(rt_df, aes(x = c24_drift_sec, y = chrysene_drift_sec, color = rt_block)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(size = 2) +
    labs(x = "N-Tetracosane-d50 EIC RT drift (sec)",
         y = "Chrysene-d12 EIC RT drift (sec)",
         color = "Automatic RT block") +
    theme_classic(base_size = 12)
  ggsave(file.path(out_dir, "auto_block_raw_internal_standard_rt_blocks.png"), p, width = 7, height = 5, dpi = 200)
  
  print(method_summary)
  print(block_summary)
  print(manual_summary)
  message("Outputs written to: ", out_dir)
  
}

run_internal_standard_guided_sample_grouping <- function(config = list()) {
  cdf_dir <- Sys.getenv("CDF_DIR", unset = config$cdf_dir %||% "")
  preblock_dir <- Sys.getenv("PREBLOCK_DIR", unset = config$preblock_dir %||% "")
  requested_mode <- Sys.getenv("ALIGNMENT_MODE", unset = config$alignment_mode %||% "auto")
  threshold <- as.integer(Sys.getenv("ALIGNMENT_SINGLE_BLOCK_THRESHOLD", unset = as.character(config$single_group_threshold %||% 25L)))
  single_dist <- as.numeric(Sys.getenv("ALIGNMENT_SINGLE_BLOCK_TIME_DIST_SEC", unset = as.character(config$single_group_time_dist_sec %||% 60)))
  if (!nzchar(cdf_dir)) stop("CDF_DIR is required.")
  if (!nzchar(preblock_dir)) stop("PREBLOCK_DIR is required.")
  dir.create(preblock_dir, recursive = TRUE, showWarnings = FALSE)
  cdf_files <- discover_alignment_cdf_files(cdf_dir, minimum = 2L)
  selected <- select_alignment_mode(requested_mode, length(cdf_files), threshold)
  writeLines(selected, file.path(preblock_dir, "alignment_mode.txt"))
  reuse <- tolower(Sys.getenv("USE_EXISTING_PREBLOCK", unset = "false")) %in% c("1", "true", "yes") ||
    tolower(Sys.getenv("USE_EXISTING_BLOCK_FOLDERS", unset = "false")) %in% c("1", "true", "yes")
  if (!reuse) {
    if (identical(selected, "single")) {
      write_single_block_preblock(preblock_dir, cdf_files, single_dist)
    } else {
      build_internal_standard_sample_groups()
    }
  }
  list(
    selected_alignment_mode = selected,
    cdf_files = cdf_files,
    group_manifest = file.path(preblock_dir, "auto_block_summary.csv"),
    group_summary = file.path(preblock_dir, "auto_block_summary.csv"),
    sample_assignment = file.path(preblock_dir, "sample_block_assignment.csv"),
    preblock_dir = preblock_dir
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(as.character(x))) y else x
