rm(list = ls())
gc()


# 0. User settings

project_root <- "G:/ML"
metadata_3d_file <- file.path(project_root, "metadata_3D_validate.csv")
output_root <- file.path(
  project_root,
  "results_3D_depth_validation_v7_FINAL_2D_PC_MATCHED_DIRECTION_CONCORDANCE"
)

# Existing 2D placenta PCA results are used as the explicit cross-dimensional reference.
# No 2D PCA is refitted here. The script reads the already disease-aligned 2D loadings.
placenta2d_result_root <- file.path(
  project_root,
  "results_pairwise_PCA_2D_placenta_300_1800_DISEASE_ALIGNED_PROMINENT_EXTREMA"
)
placenta2d_pairwise_root <- file.path(placenta2d_result_root, "05_pairwise_PCA")

analysis_min <- 300
analysis_max <- 1800
n_pcs <- 3L
pca_center <- TRUE
pca_scale <- FALSE
random_seed <- 20260807L


nonfinite_action <- "stop"  # "stop" or "drop_pixels"


require_one_sample_per_group <- TRUE
require_one_file_per_group_layer <- TRUE


peak_detection_min <- analysis_min
peak_detection_max <- analysis_max
loading_peak_sg_window <- 7L
peak_min_fraction_absmax <- 0.20
peak_min_distance_cm1 <- 20

merge_peak_distance_cm1 <- 15
consensus_max_region_width_cm1 <- 45


feature_overlap_tolerance_cm1 <- 0


target_loading_buffer_cm1 <- 5

band_edge_channel_tolerance_cm1 <- 0.10


cross_dimensional_pc_match_min_abs_correlation <- 0.60


reference_max_pixels_per_layer_per_group <- 3000L


strong_support_fraction <- 0.80
moderate_support_fraction <- 0.60

# Near-zero threshold for Disease-Control depth deltas.
# For each band, tolerance = max(absolute floor,
#                               MAD multiplier * Control layer MAD,
#                               relative fraction * median(abs(Control layer mean))).
delta_tolerance_mad_multiplier <- 0.50
delta_tolerance_relative_fraction <- 0.01
delta_tolerance_absolute_floor <- 0

# Total-intensity sensitivity/QC.
run_depth_qc <- TRUE
qc_correlation_caution_threshold <- 0.80

# Optional figures.
save_layerwise_loading_plots <- FALSE
plot_dpi <- 600

# Direct MAT-file reader. This avoids hidden preprocessing in external functions.
auto_install_packages <- TRUE

set.seed(random_seed)


# 0.1 Groups and output directories

group_order <- c("Control", "SLE", "APS", "sPE", "FGR", "ICP", "GDM")
disease_order <- c("SLE", "APS", "sPE", "FGR", "ICP", "GDM")
control_group <- "Control"

plot_group_colors <- c(
  Control = "#1F3A6D",
  SLE     = "#D74F8B",
  APS     = "#FE7882",
  sPE     = "#FF9360",
  FGR     = "#FFBE4A",
  ICP     = "#8A79AA",
  GDM     = "#1FB0B0"
)

out_qc <- file.path(output_root, "01_QC")
out_direct <- file.path(output_root, "02_DIRECT_fixed_band_depth_validation")
out_reference <- file.path(output_root, "03_FIXED_DEPTH_REFERENCE_PCA")
out_layerwise <- file.path(output_root, "04_LAYERWISE_pairwise_PCA")
out_consensus <- file.path(output_root, "05_LAYERWISE_consensus_features")
out_integrated <- file.path(output_root, "06_INTEGRATED_2D_3D_evidence")
out_figures <- file.path(output_root, "07_FIGURES")
for (d in c(output_root, out_qc, out_direct, out_reference, out_layerwise,
            out_consensus, out_integrated, out_figures)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# 0.2 Packages

required_packages <- c("data.table", "ggplot2", "R.matlab")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  if (!isTRUE(auto_install_packages)) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# 0.3 Seven predefined 2D recurrent feature regions (FINAL)

feature_2d_dt <- data.table(
  feature_label = c(
    "998.8-1002.6",
    "1151.1-1154.9",
    "1509.8-1517.0",
    "1154.9-1177.5",
    "1245.0-1286.0",
    "1440.7-1448.0",
    "1657.0-1664.1"
  ),
  feature_min = c(
    998.8, 1151.1, 1509.8,
    1154.9, 1245.0, 1440.7, 1657.0
  ),
  feature_max = c(
    1002.6, 1154.9, 1517.0,
    1177.5, 1286.0, 1448.0, 1664.1
  ),
  expected_pc_label = c(
    "PC1", "PC1", "PC1",
    "PC2", "PC2", "PC2", "PC2"
  ),
  feature_class = c(
    rep("2D recurrent PC1 feature", 3),
    rep("2D recurrent PC2 feature", 4)
  ),
  assignment = rep(
    "2D-defined recurrent Raman feature; molecular attribution interpreted separately",
    7
  )
)

if (anyDuplicated(feature_2d_dt$feature_label)) {
  stop(
    "feature_2d_dt$feature_label contains duplicated labels: ",
    paste(
      unique(feature_2d_dt$feature_label[duplicated(feature_2d_dt$feature_label)]),
      collapse = "; "
    )
  )
}

feature_export <- copy(feature_2d_dt)
fwrite(
  feature_export,
  file.path(out_qc, "INPUT_2D_recurrent_feature_regions.csv"),
  bom = TRUE
)


# 1. Utility functions

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) ||
      all(trimws(as.character(x)) == "")) y else x
}

safe_group_name <- function(x) {
  x0 <- trimws(as.character(x))
  y <- tolower(x0)
  ifelse(y %in% c("control", "ctrl", "con", "healthy", "hc"), "Control",
  ifelse(y == "sle", "SLE",
  ifelse(y == "aps", "APS",
  ifelse(y %in% c("spe", "s_pe", "s-pe", "severe_pe", "severepe"), "sPE",
  ifelse(y == "fgr", "FGR",
  ifelse(y == "icp", "ICP",
  ifelse(y == "gdm", "GDM", x0)))))))
}

infer_layer_from_text <- function(x) {
  x <- basename(as.character(x))
  hit <- regexpr("(?i)(layer|z)[_ -]*[0-9]+", x, perl = TRUE)
  if (hit[1] < 0) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", regmatches(x, hit)))
}

clean_layer_value <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))
  out <- suppressWarnings(as.integer(as.character(x)))
  bad <- is.na(out)
  if (any(bad)) out[bad] <- vapply(x[bad], infer_layer_from_text, integer(1))
  out
}

normalize_name <- function(x) tolower(gsub("[^a-z0-9]", "", as.character(x)))

as_numeric_matrix <- function(x) {
  if (is.data.frame(x) || data.table::is.data.table(x)) x <- as.matrix(x)
  if (!is.numeric(x)) suppressWarnings(storage.mode(x) <- "double")
  as.matrix(x)
}

flatten_numeric_objects <- function(x, prefix = "") {
  out <- list()
  if (is.numeric(x) && !is.null(dim(x))) {
    out[[prefix]] <- x
    return(out)
  }
  if (is.numeric(x) && is.null(dim(x)) && length(x) > 1L) {
    out[[prefix]] <- x
    return(out)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      nm <- if (nzchar(prefix)) paste0(prefix, ".", nms[i]) else nms[i]
      out <- c(out, flatten_numeric_objects(x[[i]], nm))
    }
  }
  out
}

pick_named_object <- function(objects, requested_name, candidates, require_vector = FALSE) {
  req <- as.character(requested_name %||% "")
  if (nzchar(trimws(req))) {
    exact <- which(names(objects) == req)
    if (length(exact) == 1L) return(objects[[exact]])
    norm <- which(normalize_name(names(objects)) == normalize_name(req))
    if (length(norm) == 1L) return(objects[[norm]])
    stop("MAT variable not found: ", req)
  }
  normalized <- normalize_name(names(objects))
  for (cand in candidates) {
    idx <- which(normalized == normalize_name(cand))
    for (ii in idx) {
      obj <- objects[[ii]]
      vec_ok <- is.null(dim(obj)) || sum(dim(obj) > 1L) <= 1L
      if (!require_vector || vec_ok) return(obj)
    }
  }
  NULL
}

read_mat_safely <- function(file_path) {
  withCallingHandlers(
    R.matlab::readMat(file_path),
    warning = function(w) {
      if (grepl("convertUTF8", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

read_preprocessed_mat <- function(file_path, spectra_var = "", wavenumber_var = "") {
  raw <- tryCatch(
    read_mat_safely(file_path),
    error = function(e) {
      stop(
        "Unable to read MAT file: ", file_path,
        "\nThe file may be MATLAB v7.3/HDF5. Save the processed spectra as RDS/CSV or convert MAT format.",
        "\nOriginal error: ", conditionMessage(e)
      )
    }
  )

  objects <- flatten_numeric_objects(raw)
  if (length(objects) == 0L) stop("No numeric objects found in MAT file: ", file_path)

  wn_candidates <- c(
    "wn", "wavenumber", "wavenumbers", "wl", "raman_shift", "ramanshift",
    "axis", "xaxis", "spectral_axis", "shift"
  )
  spectra_candidates <- c(
    "X", "spectra", "spectrum", "data", "raman", "raman_data",
    "processed", "processed_spectra", "spectral_data", "cube"
  )

  wn_obj <- pick_named_object(
    objects,
    requested_name = wavenumber_var,
    candidates = wn_candidates,
    require_vector = TRUE
  )

  if (is.null(wn_obj)) {
    candidates <- objects[vapply(objects, function(z) {
      dz <- dim(z)
      vec_like <- is.null(dz) || sum(dz > 1L) <= 1L
      v <- as.numeric(z)
      vec_like && length(v) >= 10L && all(is.finite(v)) &&
        min(v) >= -100 && max(v) <= 5000 && length(unique(v)) == length(v)
    }, logical(1))]
    if (length(candidates) > 0L) {
      wn_obj <- candidates[[which.max(vapply(candidates, length, integer(1)))]]
    }
  }
  if (is.null(wn_obj)) stop("Unable to identify Raman axis in: ", file_path)
  wn <- as.numeric(wn_obj)

  X_obj <- pick_named_object(
    objects,
    requested_name = spectra_var,
    candidates = spectra_candidates,
    require_vector = FALSE
  )

  if (is.null(X_obj)) {
    candidates <- objects[vapply(objects, function(z) {
      dz <- dim(z)
      !is.null(dz) && length(dz) >= 2L && any(dz == length(wn))
    }, logical(1))]
    if (length(candidates) == 0L) {
      stop("Unable to identify spectra matrix in: ", file_path)
    }
    X_obj <- candidates[[which.max(vapply(candidates, length, integer(1)))]]
  }

  dims <- dim(X_obj)
  if (length(dims) == 2L) {
    X <- as_numeric_matrix(X_obj)
    if (ncol(X) == length(wn)) {
      # already pixels x wavenumbers
    } else if (nrow(X) == length(wn)) {
      X <- t(X)
    } else {
      stop("No spectra-matrix dimension matches Raman axis: ", file_path)
    }
  } else {
    spectral_dim <- which(dims == length(wn))
    if (length(spectral_dim) != 1L) {
      stop("Cannot uniquely identify spectral dimension in: ", file_path)
    }
    arr <- aperm(X_obj, c(setdiff(seq_along(dims), spectral_dim), spectral_dim))
    X <- matrix(arr, ncol = length(wn))
  }

  ord <- order(wn)
  wn <- wn[ord]
  X <- X[, ord, drop = FALSE]
  if (anyDuplicated(wn)) stop("Duplicated Raman channels in: ", file_path)

  keep <- is.finite(wn) & wn >= analysis_min & wn <= analysis_max
  X <- X[, keep, drop = FALSE]
  wn <- wn[keep]
  if (ncol(X) < 10L) stop("Too few Raman channels in analysis range: ", file_path)

  finite_rows <- rowSums(!is.finite(X)) == 0L
  n_bad <- sum(!finite_rows)
  if (n_bad > 0L) {
    if (identical(nonfinite_action, "stop")) {
      stop("Non-finite pixels found in ", file_path, ": ", n_bad)
    }
    if (identical(nonfinite_action, "drop_pixels")) {
      X <- X[finite_rows, , drop = FALSE]
    } else {
      stop("Unknown nonfinite_action: ", nonfinite_action)
    }
  }
  if (nrow(X) < 2L) stop("Too few valid pixels after QC: ", file_path)

  list(
    X = X,
    wn = wn,
    n_pixels_before = length(finite_rows),
    n_pixels_after = nrow(X),
    n_nonfinite_pixels = n_bad
  )
}

read_layer <- function(meta_row) {
  read_preprocessed_mat(
    file_path = as.character(meta_row$mat_path),
    spectra_var = as.character(meta_row$spectra_var %||% ""),
    wavenumber_var = as.character(meta_row$wavenumber_var %||% "")
  )
}

sample_rows_matrix <- function(X, max_n) {
  if (is.null(max_n) || is.infinite(max_n) || nrow(X) <= max_n) return(X)
  X[sample.int(nrow(X), as.integer(max_n), replace = FALSE), , drop = FALSE]
}

combine_mean_M2 <- function(n_a, mean_a, M2_a, n_b, mean_b, M2_b) {
  n_a <- as.double(n_a)
  n_b <- as.double(n_b)
  if (n_a == 0) return(list(n = n_b, mean = mean_b, M2 = M2_b))
  if (n_b == 0) return(list(n = n_a, mean = mean_a, M2 = M2_a))
  n_total <- n_a + n_b
  delta <- mean_b - mean_a
  mean_total <- mean_a + delta * (n_b / n_total)
  M2_total <- M2_a + M2_b + tcrossprod(delta) * ((n_a * n_b) / n_total)
  list(n = n_total, mean = mean_total, M2 = M2_total)
}

cohens_d_two_group <- function(x_control, x_disease) {
  x_control <- as.numeric(x_control)
  x_disease <- as.numeric(x_disease)
  x_control <- x_control[is.finite(x_control)]
  x_disease <- x_disease[is.finite(x_disease)]
  n_c <- length(x_control)
  n_d <- length(x_disease)
  if (n_c < 2L || n_d < 2L) return(NA_real_)
  pooled <- sqrt(
    ((n_c - 1) * var(x_control) + (n_d - 1) * var(x_disease)) /
      (n_c + n_d - 2)
  )
  if (!is.finite(pooled) || pooled <= 0) return(NA_real_)
  (mean(x_disease) - mean(x_control)) / pooled
}

fit_pairwise_pca <- function(X_control, X_disease, wn, pair_id) {
  X_control <- as.matrix(X_control)
  X_disease <- as.matrix(X_disease)
  if (ncol(X_control) != ncol(X_disease) || ncol(X_control) != length(wn)) {
    stop("PCA dimension mismatch: ", pair_id)
  }

  mean_c <- colMeans(X_control)
  centered_c_local <- sweep(X_control, 2L, mean_c, FUN = "-")
  M2_c <- crossprod(centered_c_local)

  mean_d <- colMeans(X_disease)
  centered_d_local <- sweep(X_disease, 2L, mean_d, FUN = "-")
  M2_d <- crossprod(centered_d_local)

  combined <- combine_mean_M2(
    nrow(X_control), mean_c, M2_c,
    nrow(X_disease), mean_d, M2_d
  )
  cov_full <- combined$M2 / (combined$n - 1)
  cov_full <- (cov_full + t(cov_full)) / 2

  keep_cols <- is.finite(diag(cov_full)) & diag(cov_full) > 0
  if (sum(keep_cols) < 10L) stop("Too few non-zero-variance channels: ", pair_id)

  cov_use <- cov_full[keep_cols, keep_cols, drop = FALSE]
  eig <- eigen(cov_use, symmetric = TRUE)
  eigenvalues_all <- pmax(eig$values, 0)
  if (sum(eigenvalues_all) <= 0) stop("Total PCA variance is zero: ", pair_id)

  n_keep <- min(n_pcs, ncol(eig$vectors))
  rotation_anchored <- eig$vectors[, seq_len(n_keep), drop = FALSE]
  colnames(rotation_anchored) <- paste0("PC", seq_len(n_keep))

  anchor_sign <- rep(1, n_keep)
  for (k in seq_len(n_keep)) {
    j <- which.max(abs(rotation_anchored[, k]))
    if (rotation_anchored[j, k] < 0) {
      rotation_anchored[, k] <- -rotation_anchored[, k]
      anchor_sign[k] <- -1
    }
  }

  center_use <- combined$mean[keep_cols]
  scores_c_anchored <- sweep(X_control[, keep_cols, drop = FALSE], 2L, center_use, FUN = "-") %*% rotation_anchored
  scores_d_anchored <- sweep(X_disease[, keep_cols, drop = FALSE], 2L, center_use, FUN = "-") %*% rotation_anchored
  colnames(scores_c_anchored) <- colnames(rotation_anchored)
  colnames(scores_d_anchored) <- colnames(rotation_anchored)

  disease_sign <- rep(1, n_keep)
  rotation_aligned <- rotation_anchored
  scores_c_aligned <- scores_c_anchored
  scores_d_aligned <- scores_d_anchored
  for (k in seq_len(n_keep)) {
    if (mean(scores_d_anchored[, k]) < mean(scores_c_anchored[, k])) {
      disease_sign[k] <- -1
      rotation_aligned[, k] <- -rotation_aligned[, k]
      scores_c_aligned[, k] <- -scores_c_aligned[, k]
      scores_d_aligned[, k] <- -scores_d_aligned[, k]
    }
  }

  var_exp <- eigenvalues_all[seq_len(n_keep)] / sum(eigenvalues_all)
  sep <- rbindlist(lapply(seq_len(n_keep), function(k) {
    pc <- paste0("PC", k)
    data.table(
      pair_id = pair_id,
      PC = pc,
      variance_explained = var_exp[k],
      variance_explained_percent = 100 * var_exp[k],
      anchor_sign = anchor_sign[k],
      disease_alignment_sign = disease_sign[k],
      control_mean_score = mean(scores_c_aligned[, k]),
      disease_mean_score = mean(scores_d_aligned[, k]),
      disease_minus_control_score = mean(scores_d_aligned[, k]) - mean(scores_c_aligned[, k]),
      cohens_d_disease_vs_control = cohens_d_two_group(scores_c_aligned[, k], scores_d_aligned[, k]),
      n_control_pixels = nrow(X_control),
      n_disease_pixels = nrow(X_disease)
    )
  }))

  diagnostics <- data.table(
    diagnostic = c(
      "max_abs_rotation_crossprod_minus_identity",
      "max_abs_offdiagonal_rotated_covariance",
      "minimum_retained_eigenvalue",
      "total_variance",
      "n_fit_pixels",
      "n_wavenumbers_retained"
    ),
    value = c(
      max(abs(crossprod(rotation_aligned) - diag(n_keep))),
      max(abs(crossprod(rotation_aligned, cov_use %*% rotation_aligned) -
                diag(diag(crossprod(rotation_aligned, cov_use %*% rotation_aligned))))),
      min(eigenvalues_all[seq_len(n_keep)]),
      sum(eigenvalues_all),
      nrow(X_control) + nrow(X_disease),
      sum(keep_cols)
    )
  )

  list(
    pair_id = pair_id,
    wn = wn[keep_cols],
    keep_cols = keep_cols,
    center = center_use,
    rotation_anchored = rotation_anchored,
    rotation_aligned = rotation_aligned,
    scores_control_aligned = scores_c_aligned,
    scores_disease_aligned = scores_d_aligned,
    variance_explained = var_exp,
    separation = sep,
    diagnostics = diagnostics,
    disease_alignment_sign = disease_sign
  )
}

project_to_reference <- function(X, model) {
  X_use <- X[, model$keep_cols, drop = FALSE]
  sweep(X_use, 2L, model$center, FUN = "-") %*% model$rotation_aligned
}

support_tier_from_fraction <- function(n_support, n_total) {
  if (!is.finite(n_support) || !is.finite(n_total) || n_total <= 0) return(NA_character_)
  frac <- n_support / n_total
  if (frac >= strong_support_fraction) return("strong_depth_support")
  if (frac >= moderate_support_fraction) return("moderate_depth_support")
  "weak_or_inconsistent_depth_support"
}

longest_consecutive_run <- function(x) {
  x <- sort(unique(as.integer(x[is.finite(x)])))
  if (length(x) == 0L) return(0L)
  id <- cumsum(c(TRUE, diff(x) != 1L))
  max(tabulate(id))
}

longest_condition_run <- function(layer, condition) {
  condition <- as.logical(condition)
  condition[is.na(condition)] <- FALSE
  longest_consecutive_run(as.integer(layer)[condition])
}

select_feature_indices <- function(wn, ff, edge_tol_cm1 = 0.30) {
  # Range-based selection is safer than hard-coding channel lists.
  # edge_tol_cm1 allows rounded feature boundaries (e.g. 998.8) to match the
  # actual measured channel (e.g. 998.773) without silently dropping it.
  which(
    wn >= ff$feature_min - edge_tol_cm1 &
      wn <= ff$feature_max + edge_tol_cm1
  )
}

feature_indices_by_window <- function(wn, ff, window_type = c("strict", "buffered"),
                                      edge_tol_cm1 = 0.30) {
  window_type <- match.arg(window_type)
  buffer <- if (window_type == "buffered") target_loading_buffer_cm1 else 0
  which(
    wn >= ff$feature_min - buffer - edge_tol_cm1 &
      wn <= ff$feature_max + buffer + edge_tol_cm1
  )
}

loading_table_from_model <- function(model, disease_group, layer = NA_integer_, source) {
  dt <- data.table(
    pair_id = model$pair_id,
    disease_group = disease_group,
    layer = layer,
    source = source,
    wavenumber = model$wn
  )
  for (k in seq_len(n_pcs)) {
    pc <- paste0("PC", k)
    loading <- model$rotation_aligned[, k]
    dt[, (paste0(pc, "_loading_disease_aligned")) := loading]
    # loading^2 is retained as an AUXILIARY audit metric only.
    dt[, (paste0(pc, "_loading_squared_contribution_percent")) :=
         100 * loading^2 / sum(loading^2)]
  }
  dt
}


smooth_loading_for_peak_detection <- function(x, window = 7L) {
  x <- as.numeric(x)
  if (length(x) < 5L) return(x)

  if (window == 5L) {
    coef <- c(-3, 12, 17, 12, -3) / 35
  } else if (window == 7L) {
    coef <- c(-2, 3, 6, 7, 6, 3, -2) / 21
  } else {
    stop("loading_peak_sg_window must be 5L or 7L.")
  }

  y <- as.numeric(stats::filter(x, filter = coef, sides = 2))
  y[!is.finite(y)] <- x[!is.finite(y)]
  y
}

detect_prominent_loading_extrema <- function(
    wavenumber,
    loading,
    smooth_window = 7L,
    min_fraction_absmax = 0.20,
    min_distance_cm1 = 20) {

  wavenumber <- as.numeric(wavenumber)
  loading_raw <- as.numeric(loading)

  if (length(wavenumber) != length(loading_raw)) {
    stop("wavenumber and loading lengths differ.")
  }

  n <- length(loading_raw)
  if (n < 5L || all(!is.finite(loading_raw))) {
    return(list(all_extrema = data.table(), prominent = data.table()))
  }

  loading_smooth <- smooth_loading_for_peak_detection(
    loading_raw,
    window = smooth_window
  )

  abs_max <- max(abs(loading_smooth), na.rm = TRUE)
  if (!is.finite(abs_max) || abs_max <= 0) {
    return(list(all_extrema = data.table(), prominent = data.table()))
  }

  half_window <- as.integer((smooth_window - 1L) / 2L)
  search_idx <- seq.int(
    max(2L, half_window + 1L),
    min(n - 1L, n - half_window)
  )

  is_pos_peak <- vapply(search_idx, function(i) {
    is.finite(loading_smooth[i]) &&
      loading_smooth[i] >= loading_smooth[i - 1L] &&
      loading_smooth[i] >= loading_smooth[i + 1L] &&
      loading_smooth[i] > 0
  }, logical(1))

  is_neg_trough <- vapply(search_idx, function(i) {
    is.finite(loading_smooth[i]) &&
      loading_smooth[i] <= loading_smooth[i - 1L] &&
      loading_smooth[i] <= loading_smooth[i + 1L] &&
      loading_smooth[i] < 0
  }, logical(1))

  candidate_idx <- unique(c(
    search_idx[is_pos_peak],
    search_idx[is_neg_trough]
  ))

  if (length(candidate_idx) == 0L) {
    return(list(all_extrema = data.table(), prominent = data.table()))
  }

  all_extrema <- data.table(
    index = candidate_idx,
    wavenumber = wavenumber[candidate_idx],
    loading_raw = loading_raw[candidate_idx],
    loading_smooth = loading_smooth[candidate_idx],
    abs_loading_smooth = abs(loading_smooth[candidate_idx]),
    relative_abs_loading = abs(loading_smooth[candidate_idx]) / abs_max,
    extremum_type = ifelse(
      loading_smooth[candidate_idx] > 0,
      "positive_peak",
      "negative_trough"
    )
  )

  all_extrema[, passes_relative_threshold :=
                relative_abs_loading >= min_fraction_absmax]
  setorder(all_extrema, -abs_loading_smooth)

  candidate_prominent <- all_extrema[
    passes_relative_threshold == TRUE &
      wavenumber >= peak_detection_min &
      wavenumber <= peak_detection_max
  ]

  # Greedy spacing by descending |loading|; NO fixed Top-N.
  selected_rows <- integer(0)
  if (nrow(candidate_prominent) > 0L) {
    for (rr in seq_len(nrow(candidate_prominent))) {
      wn_i <- candidate_prominent$wavenumber[rr]
      if (length(selected_rows) == 0L ||
          all(abs(wn_i - candidate_prominent$wavenumber[selected_rows]) >=
                min_distance_cm1)) {
        selected_rows <- c(selected_rows, rr)
      }
    }
  }

  prominent <- if (length(selected_rows) > 0L) {
    candidate_prominent[selected_rows]
  } else {
    candidate_prominent[0]
  }

  if (nrow(prominent) > 0L) {
    prominent[, rank_by_abs_loading := seq_len(.N)]
    setorder(prominent, wavenumber)
  }

  list(all_extrema = all_extrema, prominent = prominent)
}

extract_all_prominent_extrema <- function(loading_dt, delta_spectrum = NULL) {
  peaks <- rbindlist(lapply(seq_len(n_pcs), function(k) {
    pc <- paste0("PC", k)
    loading_col <- paste0(pc, "_loading_disease_aligned")
    contribution_col <- paste0(pc, "_loading_squared_contribution_percent")

    detected <- detect_prominent_loading_extrema(
      wavenumber = loading_dt$wavenumber,
      loading = loading_dt[[loading_col]],
      smooth_window = loading_peak_sg_window,
      min_fraction_absmax = peak_min_fraction_absmax,
      min_distance_cm1 = peak_min_distance_cm1
    )

    out <- detected$prominent
    if (nrow(out) == 0L) return(data.table())

    idx <- out$index
    out[, `:=`(
      pair_id = loading_dt$pair_id[1],
      disease_group = loading_dt$disease_group[1],
      layer = loading_dt$layer[1],
      source = loading_dt$source[1],
      PC = pc,
      loading_disease_aligned = loading_raw,
      loading_squared_contribution_percent = loading_dt[[contribution_col]][idx],
      loading_side = fifelse(
        loading_raw > 0,
        "disease_side_positive_loading",
        fifelse(loading_raw < 0, "control_side_negative_loading", "near_zero_loading")
      ),
      relative_threshold_used = peak_min_fraction_absmax,
      min_distance_cm1_used = peak_min_distance_cm1,
      SG_window_used = loading_peak_sg_window
    )]

    if (!is.null(delta_spectrum)) {
      delta <- delta_spectrum[idx]
      out[, disease_minus_control_intensity := delta]
      out[, intensity_direction := fifelse(
        delta > 0, "higher_in_disease",
        fifelse(delta < 0, "lower_in_disease", "no_mean_difference")
      )]
      out[, centroid_separation_contribution :=
            disease_minus_control_intensity * loading_disease_aligned]
      out[, separation_role := fifelse(
        centroid_separation_contribution > 0,
        "supports_disease_aligned_separation",
        fifelse(
          centroid_separation_contribution < 0,
          "opposes_disease_aligned_separation",
          "neutral"
        )
      )]
    }

    out
  }), fill = TRUE)

  peaks
}

pc_permutations_3 <- rbind(
  c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
  c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
)

match_layer_pcs_to_reference <- function(layer_loading_dt, reference_loading_dt) {
  merged <- merge(
    layer_loading_dt,
    reference_loading_dt,
    by = "wavenumber",
    suffixes = c("_layer", "_reference")
  )
  if (nrow(merged) < 10L) stop("Too few common channels for PC matching.")

  cor_mat <- matrix(
    NA_real_,
    nrow = n_pcs,
    ncol = n_pcs,
    dimnames = list(paste0("PC", 1:n_pcs), paste0("PC", 1:n_pcs))
  )

  for (i in seq_len(n_pcs)) {
    for (j in seq_len(n_pcs)) {
      cor_mat[i, j] <- suppressWarnings(cor(
        merged[[paste0("PC", i, "_loading_disease_aligned_layer")]],
        merged[[paste0("PC", j, "_loading_disease_aligned_reference")]],
        use = "complete.obs",
        method = "pearson"
      ))
    }
  }

  perm_scores <- apply(pc_permutations_3, 1, function(p) {
    sum(abs(cor_mat[cbind(seq_len(n_pcs), p)]), na.rm = TRUE)
  })
  best_perm <- pc_permutations_3[which.max(perm_scores), ]

  data.table(
    layer_PC = paste0("PC", seq_len(n_pcs)),
    reference_PC = paste0("PC", best_perm),
    signed_loading_correlation = cor_mat[cbind(seq_len(n_pcs), best_perm)],
    absolute_loading_correlation = abs(cor_mat[cbind(seq_len(n_pcs), best_perm)]),
    match_quality = fifelse(
      abs(cor_mat[cbind(seq_len(n_pcs), best_perm)]) >= 0.80, "strong",
      fifelse(
        abs(cor_mat[cbind(seq_len(n_pcs), best_perm)]) >= 0.60,
        "moderate",
        "weak"
      )
    ),
    total_assignment_score = max(perm_scores, na.rm = TRUE)
  )
}


read_existing_2d_loading <- function(disease_group) {
  f <- file.path(
    placenta2d_pairwise_root,
    paste0("Control_vs_", disease_group),
    "02_tables",
    "PCA_PC1_PC3_loadings_contributions_and_direction_audit.csv"
  )
  if (!file.exists(f)) {
    stop(
      "Cannot find existing 2D placenta loading table for ", disease_group,
      ":\n", f,
      "\nRun the finalized 2D pairwise PCA first or correct placenta2d_result_root."
    )
  }
  dt <- fread(f)
  req <- c(
    "wavenumber",
    paste0("PC", 1:n_pcs, "_loading_disease_aligned"),
    paste0("PC", 1:n_pcs, "_loading_squared_contribution_percent")
  )
  if (!all(req %in% names(dt))) {
    stop("2D loading table lacks required columns: ", f)
  }
  dt[, `:=`(
    pair_id = paste0("Control_vs_", disease_group, "_2D"),
    disease_group = disease_group,
    layer = NA_integer_,
    source = "existing_2D_placenta_PCA"
  )]
  dt
}

match_reference_pcs_to_2d <- function(reference_loading_dt, loading_2d_dt,
                                      axis_tolerance_cm1 = 0.50) {
  wn3 <- as.numeric(reference_loading_dt$wavenumber)
  wn2 <- as.numeric(loading_2d_dt$wavenumber)

  nearest2 <- vapply(wn3, function(x) {
    if (!is.finite(x)) return(NA_integer_)
    which.min(abs(wn2 - x))
  }, integer(1))
  axis_diff <- rep(NA_real_, length(wn3))
  ok_idx <- which(!is.na(nearest2))
  axis_diff[ok_idx] <- abs(wn3[ok_idx] - wn2[nearest2[ok_idx]])
  keep <- which(is.finite(axis_diff) & axis_diff <= axis_tolerance_cm1)

  min_required <- max(
    10L,
    floor(0.80 * min(nrow(reference_loading_dt), nrow(loading_2d_dt)))
  )
  if (length(keep) < min_required) {
    stop(
      "Too few Raman channels for 3D-reference <-> 2D PC matching: ",
      length(keep), " channels within ±", axis_tolerance_cm1, " cm^-1."
    )
  }

  cor_mat <- matrix(
    NA_real_, nrow = n_pcs, ncol = n_pcs,
    dimnames = list(paste0("3Dref_", 1:n_pcs), paste0("2D_", 1:n_pcs))
  )
  for (i in seq_len(n_pcs)) {
    x3 <- reference_loading_dt[[paste0("PC", i, "_loading_disease_aligned")]][keep]
    for (j in seq_len(n_pcs)) {
      x2 <- loading_2d_dt[[paste0("PC", j, "_loading_disease_aligned")]][nearest2[keep]]
      cor_mat[i, j] <- suppressWarnings(cor(
        x3, x2,
        use = "complete.obs",
        method = "pearson"
      ))
    }
  }

  perm_scores <- apply(pc_permutations_3, 1, function(p) {
    sum(abs(cor_mat[cbind(seq_len(n_pcs), p)]), na.rm = TRUE)
  })
  best_perm <- pc_permutations_3[which.max(perm_scores), ]
  matched_cor <- cor_mat[cbind(seq_len(n_pcs), best_perm)]

  data.table(
    reference_PC = paste0("PC", seq_len(n_pcs)),
    matched_2D_PC = paste0("PC", best_perm),
    signed_loading_correlation_3Dref_to_2D = matched_cor,
    absolute_loading_correlation_3Dref_to_2D = abs(matched_cor),
    match_quality_3Dref_to_2D = fifelse(
      abs(matched_cor) >= 0.80, "strong",
      fifelse(abs(matched_cor) >= 0.60, "moderate", "weak")
    ),
    total_assignment_score_3Dref_to_2D = max(perm_scores, na.rm = TRUE),
    n_axis_channels_matched = length(keep),
    max_axis_difference_cm1 = max(axis_diff[keep], na.rm = TRUE),
    axis_tolerance_cm1 = axis_tolerance_cm1
  )
}

build_2d_feature_direction_table <- function(loading_2d_list, feature_dt) {
  out <- list()
  edge_tol <- band_edge_channel_tolerance_cm1

  for (dg in names(loading_2d_list)) {
    dt <- loading_2d_list[[dg]]
    detected <- extract_all_prominent_extrema(dt, delta_spectrum = NULL)

    for (f in seq_len(nrow(feature_dt))) {
      ff <- feature_dt[f]
      hits <- detected[
        PC == ff$expected_pc_label &
          wavenumber >= ff$feature_min - edge_tol &
          wavenumber <= ff$feature_max + edge_tol
      ]

      if (nrow(hits) > 0L) {
        strongest <- hits[which.max(relative_abs_loading)]
        signs <- sign(hits$loading_disease_aligned)
        signs <- signs[is.finite(signs) & signs != 0]
        n_pos <- sum(signs > 0)
        n_neg <- sum(signs < 0)
        direction_pattern <- if (n_pos > 0L && n_neg == 0L) {
          "positive_only"
        } else if (n_neg > 0L && n_pos == 0L) {
          "negative_only"
        } else if (n_pos > 0L && n_neg > 0L) {
          "mixed"
        } else {
          "near_zero_or_undefined"
        }
        direction_sign <- if (direction_pattern == "positive_only") {
          1
        } else if (direction_pattern == "negative_only") {
          -1
        } else {
          NA_real_
        }

        out[[length(out) + 1L]] <- data.table(
          disease_group = dg,
          feature_label = ff$feature_label,
          source_2D_PC = ff$expected_pc_label,
          twoD_formal_prominent_support = TRUE,
          twoD_n_prominent_extrema_in_window = nrow(hits),
          twoD_n_positive_extrema = n_pos,
          twoD_n_negative_extrema = n_neg,
          twoD_direction_pattern = direction_pattern,
          twoD_direction_evaluable = is.finite(direction_sign),
          twoD_loading_sign = direction_sign,
          twoD_peak_wavenumber_strongest = strongest$wavenumber,
          twoD_extremum_type_strongest = strongest$extremum_type,
          twoD_loading_strongest = strongest$loading_disease_aligned,
          twoD_relative_abs_loading_strongest = strongest$relative_abs_loading
        )
      } else {
        # Do NOT invent a direction when the finalized 2D detector did not call
        # a formal prominent extremum for that disease/feature.
        out[[length(out) + 1L]] <- data.table(
          disease_group = dg,
          feature_label = ff$feature_label,
          source_2D_PC = ff$expected_pc_label,
          twoD_formal_prominent_support = FALSE,
          twoD_n_prominent_extrema_in_window = 0L,
          twoD_n_positive_extrema = 0L,
          twoD_n_negative_extrema = 0L,
          twoD_direction_pattern = "no_formal_extremum",
          twoD_direction_evaluable = FALSE,
          twoD_loading_sign = NA_real_,
          twoD_peak_wavenumber_strongest = NA_real_,
          twoD_extremum_type_strongest = NA_character_,
          twoD_loading_strongest = NA_real_,
          twoD_relative_abs_loading_strongest = NA_real_
        )
      }
    }
  }
  rbindlist(out, fill = TRUE)
}

safe_max_numeric <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

safe_median_numeric <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

# Target support is defined by the PRESENCE OF A PROMINENT SIGNED-LOADING
# EXTREMUM. PC identity is NOT assumed from the numeric PC label. Instead:
# layer PC -> fixed 3D reference PC -> matched source 2D PC.
# Loading-direction concordance is reported separately from feature recurrence.
summarize_target_support_one_loading <- function(
    loading_dt,
    feature_dt,
    reference_to_2d_mapping,
    feature_direction_dt,
    pc_mapping = NULL,
    window_type = c("strict", "buffered")) {

  window_type <- match.arg(window_type)
  buffer <- if (window_type == "buffered") target_loading_buffer_cm1 else 0
  edge_tol <- band_edge_channel_tolerance_cm1
  dg <- as.character(loading_dt$disease_group[1])

  detected_by_pc <- lapply(seq_len(n_pcs), function(k) {
    pc <- paste0("PC", k)
    detect_prominent_loading_extrema(
      wavenumber = loading_dt$wavenumber,
      loading = loading_dt[[paste0(pc, "_loading_disease_aligned")]],
      smooth_window = loading_peak_sg_window,
      min_fraction_absmax = peak_min_fraction_absmax,
      min_distance_cm1 = peak_min_distance_cm1
    )$prominent
  })

  out <- list()

  for (f in seq_len(nrow(feature_dt))) {
    ff <- feature_dt[f]
    source_2d_pc <- ff$expected_pc_label
    low <- ff$feature_min - buffer - edge_tol
    high <- ff$feature_max + buffer + edge_tol

    dir2d <- feature_direction_dt[
      disease_group == dg & feature_label == ff$feature_label
    ]
    twoD_sign <- if (nrow(dir2d) == 1L) dir2d$twoD_loading_sign[1] else NA_real_
    twoD_peak <- if (nrow(dir2d) == 1L) dir2d$twoD_peak_wavenumber_strongest[1] else NA_real_
    twoD_type <- if (nrow(dir2d) == 1L) dir2d$twoD_extremum_type_strongest[1] else NA_character_
    twoD_pattern <- if (nrow(dir2d) == 1L) dir2d$twoD_direction_pattern[1] else NA_character_
    twoD_direction_evaluable <- nrow(dir2d) == 1L && isTRUE(dir2d$twoD_direction_evaluable[1])
    twoD_formal <- nrow(dir2d) == 1L && isTRUE(dir2d$twoD_formal_prominent_support[1])

    for (k in seq_len(n_pcs)) {
      pc <- paste0("PC", k)
      reference_pc <- pc
      corr_layer_to_ref <- NA_real_

      if (!is.null(pc_mapping) && nrow(pc_mapping) > 0L) {
        reference_pc <- pc_mapping[layer_PC == pc, reference_PC][1]
        corr_layer_to_ref <- pc_mapping[
          layer_PC == pc,
          signed_loading_correlation
        ][1]
      }

      map_row <- reference_to_2d_mapping[reference_PC == reference_pc]
      matched_2d_pc <- if (nrow(map_row) == 1L) map_row$matched_2D_PC[1] else NA_character_
      corr_ref_to_2d <- if (nrow(map_row) == 1L) {
        map_row$signed_loading_correlation_3Dref_to_2D[1]
      } else NA_real_
      abs_corr_ref_to_2d <- if (nrow(map_row) == 1L) {
        map_row$absolute_loading_correlation_3Dref_to_2D[1]
      } else NA_real_
      cross_dimensional_match_evaluable <- is.finite(abs_corr_ref_to_2d) &&
        abs_corr_ref_to_2d >= cross_dimensional_pc_match_min_abs_correlation

      hits <- detected_by_pc[[k]][wavenumber >= low & wavenumber <= high]
      supported <- nrow(hits) > 0L
      strongest <- if (supported) hits[which.max(relative_abs_loading)] else NULL
      loading3d <- if (supported) strongest$loading_raw else NA_real_
      strongest_sign3d <- if (is.finite(loading3d)) sign(loading3d) else NA_real_

      hit_signs <- if (supported) sign(hits$loading_raw) else numeric(0)
      hit_signs <- hit_signs[is.finite(hit_signs) & hit_signs != 0]
      n_pos_3d <- sum(hit_signs > 0)
      n_neg_3d <- sum(hit_signs < 0)
      direction_pattern_3d <- if (!supported) {
        "no_formal_extremum"
      } else if (n_pos_3d > 0L && n_neg_3d == 0L) {
        "positive_only"
      } else if (n_neg_3d > 0L && n_pos_3d == 0L) {
        "negative_only"
      } else if (n_pos_3d > 0L && n_neg_3d > 0L) {
        "mixed"
      } else {
        "near_zero_or_undefined"
      }
      sign3d <- if (direction_pattern_3d == "positive_only") {
        1
      } else if (direction_pattern_3d == "negative_only") {
        -1
      } else {
        NA_real_
      }

      sign_match <- supported && twoD_direction_evaluable &&
        is.finite(twoD_sign) && is.finite(sign3d) && sign3d == twoD_sign
      strongest_sign_match <- supported && twoD_direction_evaluable &&
        is.finite(twoD_sign) && is.finite(strongest_sign3d) && strongest_sign3d == twoD_sign
      expected_2d_pc <- !is.na(matched_2d_pc) && matched_2d_pc == source_2d_pc

      out[[length(out) + 1L]] <- data.table(
        feature_label = ff$feature_label,
        feature_min = ff$feature_min,
        feature_max = ff$feature_max,
        expected_pc_label = source_2d_pc,
        pair_id = loading_dt$pair_id[1],
        disease_group = dg,
        layer = loading_dt$layer[1],
        source = loading_dt$source[1],
        window_type = window_type,
        layer_PC = pc,
        matched_reference_PC = reference_pc,
        matched_2D_PC = matched_2d_pc,
        signed_loading_correlation_to_reference = corr_layer_to_ref,
        signed_loading_correlation_reference_to_2D = corr_ref_to_2d,
        absolute_loading_correlation_reference_to_2D = abs_corr_ref_to_2d,
        cross_dimensional_PC_match_evaluable = cross_dimensional_match_evaluable,
        expected_reference_PC = expected_2d_pc,
        n_prominent_extrema_in_window = nrow(hits),
        supported = supported,
        strongest_relative_abs_loading = if (supported) strongest$relative_abs_loading else NA_real_,
        wavenumber_at_strongest_extremum = if (supported) strongest$wavenumber else NA_real_,
        extremum_type_at_strongest = if (supported) strongest$extremum_type else NA_character_,
        loading_at_strongest_extremum = loading3d,
        loading_sign_at_strongest_extremum = strongest_sign3d,
        n_positive_extrema_in_window = n_pos_3d,
        n_negative_extrema_in_window = n_neg_3d,
        loading_direction_pattern_in_window = direction_pattern_3d,
        loading_direction_sign_in_window = sign3d,
        loading_squared_contribution_percent_at_strongest = if (supported) {
          idx <- strongest$index
          loading_dt[[paste0(pc, "_loading_squared_contribution_percent")]][idx]
        } else NA_real_,
        twoD_formal_prominent_support = twoD_formal,
        twoD_peak_wavenumber = twoD_peak,
        twoD_extremum_type = twoD_type,
        twoD_direction_pattern = twoD_pattern,
        twoD_direction_evaluable = twoD_direction_evaluable,
        twoD_loading_sign = twoD_sign,
        loading_direction_concordant_with_2D = sign_match,
        strongest_extremum_direction_concordant_with_2D = strongest_sign_match,
        expected_PC_and_direction_concordant = expected_2d_pc && sign_match,
        expected_PC_and_direction_concordant_high_confidence =
          expected_2d_pc && sign_match && cross_dimensional_match_evaluable
      )
    }
  }

  rbindlist(out, fill = TRUE)
}

collapse_target_support_by_layer <- function(point_dt) {
  if (nrow(point_dt) == 0L) return(data.table())

  point_dt[, {
    strength <- strongest_relative_abs_loading
    any_finite <- any(is.finite(strength))
    j_any <- if (any_finite) which.max(replace(strength, !is.finite(strength), -Inf)) else NA_integer_

    exp_rows <- expected_reference_PC == TRUE
    exp_supported_rows <- which(exp_rows & supported)
    j_exp <- if (length(exp_supported_rows) > 0L) {
      z <- strongest_relative_abs_loading[exp_supported_rows]
      exp_supported_rows[which.max(replace(z, !is.finite(z), -Inf))]
    } else NA_integer_

    expected_support <- length(exp_supported_rows) > 0L
    direction_evaluable <- expected_support && any(
      exp_rows & supported & twoD_direction_evaluable &
        is.finite(loading_direction_sign_in_window),
      na.rm = TRUE
    )
    direction_concordant <- any(
      expected_PC_and_direction_concordant,
      na.rm = TRUE
    )
    direction_evaluable_high_confidence <- expected_support && any(
      exp_rows & supported & twoD_direction_evaluable &
        is.finite(loading_direction_sign_in_window) &
        cross_dimensional_PC_match_evaluable,
      na.rm = TRUE
    )
    direction_concordant_high_confidence <- any(
      expected_PC_and_direction_concordant_high_confidence,
      na.rm = TRUE
    )

    data.table(
      feature_min = feature_min[1],
      feature_max = feature_max[1],
      expected_pc_label = expected_pc_label[1],
      source = source[1],
      any_PC_support = any(supported, na.rm = TRUE),
      expected_reference_PC_support = expected_support,
      expected_reference_PC_direction_evaluable = direction_evaluable,
      expected_reference_PC_direction_concordant = direction_concordant,
      expected_reference_PC_direction_evaluable_high_confidence = direction_evaluable_high_confidence,
      expected_reference_PC_direction_concordant_high_confidence = direction_concordant_high_confidence,
      expected_3D_reference_PC = {
        z <- unique(matched_reference_PC[expected_reference_PC == TRUE])
        if (length(z)) z[1] else NA_character_
      },
      layer_PC_matched_to_2D_source_PC = if (is.finite(j_exp)) layer_PC[j_exp] else NA_character_,
      matched_2D_PC_at_expected_support = if (is.finite(j_exp)) matched_2D_PC[j_exp] else NA_character_,
      layer_PC_with_max_support = if (is.finite(j_any)) layer_PC[j_any] else NA_character_,
      matched_reference_PC_with_max_support = if (is.finite(j_any)) matched_reference_PC[j_any] else NA_character_,
      maximum_relative_abs_loading = safe_max_numeric(strongest_relative_abs_loading),
      wavenumber_at_max_support = if (is.finite(j_any)) wavenumber_at_strongest_extremum[j_any] else NA_real_,
      extremum_type_at_max_support = if (is.finite(j_any)) extremum_type_at_strongest[j_any] else NA_character_,
      loading_at_max_support = if (is.finite(j_any)) loading_at_strongest_extremum[j_any] else NA_real_,
      max_expected_reference_PC_relative_abs_loading = safe_max_numeric(
        strongest_relative_abs_loading[exp_rows]
      ),
      wavenumber_at_expected_PC_support = if (is.finite(j_exp)) wavenumber_at_strongest_extremum[j_exp] else NA_real_,
      extremum_type_at_expected_PC_support = if (is.finite(j_exp)) extremum_type_at_strongest[j_exp] else NA_character_,
      loading_at_expected_PC_support = if (is.finite(j_exp)) loading_at_strongest_extremum[j_exp] else NA_real_,
      loading_sign_at_expected_PC_support = if (is.finite(j_exp)) loading_direction_sign_in_window[j_exp] else NA_real_,
      loading_direction_pattern_at_expected_PC_support = if (is.finite(j_exp)) loading_direction_pattern_in_window[j_exp] else NA_character_,
      twoD_formal_prominent_support = any(twoD_formal_prominent_support, na.rm = TRUE),
      twoD_peak_wavenumber = {
        z <- unique(twoD_peak_wavenumber[is.finite(twoD_peak_wavenumber)])
        if (length(z)) z[1] else NA_real_
      },
      twoD_extremum_type = {
        z <- unique(twoD_extremum_type[!is.na(twoD_extremum_type)])
        if (length(z)) z[1] else NA_character_
      },
      twoD_direction_pattern = {
        z <- unique(twoD_direction_pattern[!is.na(twoD_direction_pattern)])
        if (length(z)) z[1] else NA_character_
      },
      twoD_loading_sign = {
        z <- unique(twoD_loading_sign[is.finite(twoD_loading_sign)])
        if (length(z)) z[1] else NA_real_
      }
    )
  }, by = .(feature_label, pair_id, disease_group, layer, window_type)]
}

summarize_target_support_across_layers <- function(layer_dt) {
  if (nrow(layer_dt) == 0L) return(data.table())

  layer_dt[, {
    n_valid <- uniqueN(layer)
    n_any <- sum(any_PC_support, na.rm = TRUE)
    n_expected <- sum(expected_reference_PC_support, na.rm = TRUE)
    n_dir_eval <- sum(expected_reference_PC_direction_evaluable, na.rm = TRUE)
    n_dir_same <- sum(expected_reference_PC_direction_concordant, na.rm = TRUE)
    n_dir_eval_hc <- sum(expected_reference_PC_direction_evaluable_high_confidence, na.rm = TRUE)
    n_dir_same_hc <- sum(expected_reference_PC_direction_concordant_high_confidence, na.rm = TRUE)

    data.table(
      feature_min = feature_min[1],
      feature_max = feature_max[1],
      expected_pc_label = expected_pc_label[1],
      source = source[1],
      n_layers_valid = n_valid,
      any_PC_layers_supported = n_any,
      any_PC_layer_fraction = n_any / n_valid,
      expected_reference_PC_layers_supported = n_expected,
      expected_reference_PC_layer_fraction = n_expected / n_valid,
      any_PC_support_tier = support_tier_from_fraction(n_any, n_valid),
      expected_reference_PC_support_tier = support_tier_from_fraction(n_expected, n_valid),
      longest_any_PC_support_run = longest_condition_run(layer, any_PC_support),
      longest_expected_reference_PC_support_run = longest_condition_run(
        layer,
        expected_reference_PC_support
      ),
      median_maximum_relative_abs_loading = safe_median_numeric(maximum_relative_abs_loading),
      median_expected_reference_PC_relative_abs_loading = safe_median_numeric(
        max_expected_reference_PC_relative_abs_loading
      ),
      loading_direction_layers_evaluable = n_dir_eval,
      loading_direction_layers_concordant = n_dir_same,
      loading_direction_concordance_fraction_among_evaluable = if (n_dir_eval > 0L) {
        n_dir_same / n_dir_eval
      } else NA_real_,
      loading_direction_concordance_fraction_all_layers = if (n_valid > 0L) {
        n_dir_same / n_valid
      } else NA_real_,
      loading_direction_layers_evaluable_high_confidence = n_dir_eval_hc,
      loading_direction_layers_concordant_high_confidence = n_dir_same_hc,
      loading_direction_concordance_fraction_high_confidence = if (n_dir_eval_hc > 0L) {
        n_dir_same_hc / n_dir_eval_hc
      } else NA_real_,
      longest_direction_concordant_run = longest_condition_run(
        layer,
        expected_reference_PC_direction_concordant
      )
    )
  }, by = .(feature_label, disease_group, window_type)]
}

# Recurrent regions are clustered separately by the PC identity matched to the 2D placenta PCA, using
# the same 15 cm^-1 adjacency rule as the revised 2D workflow.
merge_consensus_regions <- function(peak_dt) {
  if (nrow(peak_dt) == 0L) return(data.table())

  dt <- copy(peak_dt)
  dt <- dt[
    is.finite(wavenumber) &
      wavenumber >= peak_detection_min &
      wavenumber <= peak_detection_max &
      !is.na(matched_2D_PC)
  ]
  if (nrow(dt) == 0L) return(data.table())

  dt[, matched_2D_PC := as.character(matched_2D_PC)]
  setorder(dt, matched_2D_PC, wavenumber)

 
  dt[, region_number := {
    wn <- as.numeric(wavenumber)
    rid <- integer(.N)
    if (.N > 0L) {
      rid[1L] <- 1L
      current <- 1L
      region_start <- wn[1L]
      if (.N > 1L) {
        for (ii in 2L:.N) {
          gap <- wn[ii] - wn[ii - 1L]
          span_if_added <- wn[ii] - region_start
          if (!is.finite(gap) ||
              gap > merge_peak_distance_cm1 ||
              span_if_added > consensus_max_region_width_cm1) {
            current <- current + 1L
            region_start <- wn[ii]
          }
          rid[ii] <- current
        }
      }
    }
    rid
  }, by = matched_2D_PC]

  dt[, region_id := sprintf(
    "%s_R%02d",
    matched_2D_PC,
    region_number
  )]

  out <- dt[, .(
    region_min = min(wavenumber, na.rm = TRUE),
    region_max = max(wavenumber, na.rm = TRUE),
    region_center = mean(range(wavenumber, na.rm = TRUE)),
    dominant_wavenumber = wavenumber[which.max(relative_abs_loading)],
    n_hits = .N,
    n_disease_layer_pairs = uniqueN(paste(disease_group, layer, sep = "_")),
    n_diseases = uniqueN(disease_group),
    diseases = paste(sort(unique(disease_group)), collapse = ";"),
    n_unique_layer_numbers = uniqueN(layer),
    layers = paste(sort(unique(layer)), collapse = ";"),
    raw_layer_PCs = paste(sort(unique(PC)), collapse = ";"),
    n_positive = sum(extremum_type == "positive_peak", na.rm = TRUE),
    n_negative = sum(extremum_type == "negative_trough", na.rm = TRUE),
    max_relative_abs_loading = max(relative_abs_loading, na.rm = TRUE),
    median_relative_abs_loading = median(relative_abs_loading, na.rm = TRUE),
    disease_side_hits = sum(
      loading_side == "disease_side_positive_loading",
      na.rm = TRUE
    ),
    control_side_hits = sum(
      loading_side == "control_side_negative_loading",
      na.rm = TRUE
    ),
    higher_intensity_hits = sum(intensity_direction == "higher_in_disease", na.rm = TRUE),
    lower_intensity_hits = sum(intensity_direction == "lower_in_disease", na.rm = TRUE),
    separation_support_hits = sum(
      separation_role == "supports_disease_aligned_separation",
      na.rm = TRUE
    )
  ), by = .(matched_2D_PC, region_id)]

  out[, `:=`(
    region_width = region_max - region_min,
    loading_side_consistency = pmax(disease_side_hits, control_side_hits) / pmax(n_hits, 1),
    intensity_direction_consistency = pmax(higher_intensity_hits, lower_intensity_hits) / pmax(n_hits, 1),
    separation_support_fraction = separation_support_hits / pmax(n_hits, 1)
  )]

  out[, `:=`(
    overlaps_2D_feature = FALSE,
    overlapping_2D_features = NA_character_,
    nearest_2D_feature = NA_character_,
    distance_to_nearest_2D_center = NA_real_
  )]

  for (i in seq_len(nrow(out))) {
    candidate_features <- feature_2d_dt[
      expected_pc_label == out$matched_2D_PC[i]
    ]

    if (nrow(candidate_features) == 0L) next

    overlap <- which(
      out$region_max[i] >= candidate_features$feature_min - feature_overlap_tolerance_cm1 &
        out$region_min[i] <= candidate_features$feature_max + feature_overlap_tolerance_cm1
    )

    centers <- (candidate_features$feature_min + candidate_features$feature_max) / 2
    nearest <- which.min(abs(out$region_center[i] - centers))

    out$overlaps_2D_feature[i] <- length(overlap) > 0L
    out$overlapping_2D_features[i] <- if (length(overlap) > 0L) {
      paste(candidate_features$feature_label[overlap], collapse = ";")
    } else NA_character_
    out$nearest_2D_feature[i] <- candidate_features$feature_label[nearest]
    out$distance_to_nearest_2D_center[i] <- abs(
      out$region_center[i] - centers[nearest]
    )
  }

  out[, consensus_score :=
        n_disease_layer_pairs * 1.0 +
        n_diseases * 2.0 +
        median_relative_abs_loading * 2.0 +
        loading_side_consistency * 2.0 +
        separation_support_fraction * 2.0 +
        fifelse(overlaps_2D_feature, 5.0, 0.0)]

  setorder(
    out,
    matched_2D_PC,
    -n_diseases,
    -n_disease_layer_pairs,
    -max_relative_abs_loading,
    region_center
  )

  out
}


# 1.1 Read existing 2D placenta loadings and define disease-specific 2D directions

if (!dir.exists(placenta2d_pairwise_root)) {
  stop(
    "Cannot find finalized 2D placenta pairwise PCA results:\n",
    placenta2d_pairwise_root
  )
}

loading_2d_list <- setNames(lapply(disease_order, read_existing_2d_loading), disease_order)
twoD_feature_direction_dt <- build_2d_feature_direction_table(
  loading_2d_list,
  feature_2d_dt
)
fwrite(
  twoD_feature_direction_dt,
  file.path(out_qc, "INPUT_2D_feature_loading_direction_by_disease.csv"),
  bom = TRUE
)


# 2. Read and audit metadata

if (!file.exists(metadata_3d_file)) {
  stop("Cannot find metadata file: ", metadata_3d_file)
}
meta <- fread(metadata_3d_file, encoding = "UTF-8")
required_cols <- c("sample_id", "group", "mat_file")
if (!all(required_cols %in% names(meta))) {
  stop("metadata_3D_validate.csv must contain: ", paste(required_cols, collapse = ", "))
}
for (nm in c("spectra_var", "wavenumber_var")) {
  if (!nm %in% names(meta)) meta[, (nm) := ""]
}

meta[, `:=`(
  sample_id = as.character(sample_id),
  group = safe_group_name(group),
  mat_file = as.character(mat_file)
)]
meta <- meta[group %in% group_order]
meta[, mat_path := mat_file]
meta[!file.exists(mat_path), mat_path := file.path(project_root, mat_file)]
if (any(!file.exists(meta$mat_path))) {
  stop("Missing 3D files:\n", paste(meta[!file.exists(mat_path), mat_path], collapse = "\n"))
}

if ("layer" %in% names(meta)) {
  meta[, layer := clean_layer_value(layer)]
} else if ("z" %in% names(meta)) {
  meta[, layer := clean_layer_value(z)]
} else if ("z_layer" %in% names(meta)) {
  meta[, layer := clean_layer_value(z_layer)]
} else {
  meta[, layer := vapply(mat_file, infer_layer_from_text, integer(1))]
}
if (any(is.na(meta$layer))) stop("Unable to determine layer index for some files.")

meta[, group := factor(group, levels = group_order)]
setorder(meta, group, sample_id, layer)

sample_count <- meta[, .(n_samples = uniqueN(sample_id)), by = group]
if (isTRUE(require_one_sample_per_group) && any(sample_count$n_samples != 1L)) {
  stop(
    "This script is configured for one placenta per group. Problem groups:\n",
    paste(sample_count[n_samples != 1L, paste(group, n_samples, sep = ": ")], collapse = "\n")
  )
}

group_layer_count <- meta[, .N, by = .(group, layer)]
if (isTRUE(require_one_file_per_group_layer) && any(group_layer_count$N != 1L)) {
  stop(
    "Each group-layer must map to exactly one file. Problem rows:\n",
    paste(capture.output(print(group_layer_count[N != 1L])), collapse = "\n")
  )
}

fwrite(meta, file.path(out_qc, "metadata_3D_used.csv"), bom = TRUE)
metadata_audit <- meta[, .(
  n_files = .N,
  n_samples = uniqueN(sample_id),
  n_layers = uniqueN(layer),
  layers = paste(sort(unique(layer)), collapse = ";"),
  sample_ids = paste(sort(unique(sample_id)), collapse = ";")
), by = group]
fwrite(metadata_audit, file.path(out_qc, "metadata_3D_layer_audit.csv"), bom = TRUE)

available_diseases <- intersect(disease_order, as.character(unique(meta[group != control_group, group])))
if (length(available_diseases) == 0L) stop("No disease groups found in metadata.")
control_meta <- meta[group == control_group]
if (nrow(control_meta) == 0L) stop("No Control group found in metadata.")

available_pairs <- rbindlist(lapply(available_diseases, function(dg) {
  dl <- sort(unique(meta[group == dg, layer]))
  cl <- sort(unique(control_meta$layer))
  common <- intersect(cl, dl)
  data.table(disease_group = dg, layer = common)
}))
fwrite(available_pairs, file.path(out_qc, "available_Control_Disease_layer_pairs.csv"), bom = TRUE)


# 3. Read/QC audit and direct fixed-band depth profiles (primary evidence)

profile_list <- list()
single_channel_list <- list()
layer_qc_list <- list()
axis_audit_list <- list()
wn_reference <- NULL

for (i in seq_len(nrow(meta))) {
  row_i <- meta[i]
  message("Direct depth profile: ", as.character(row_i$group), " layer ", row_i$layer)
  obj <- read_layer(row_i)

  if (is.null(wn_reference)) {
    wn_reference <- obj$wn
  } else if (length(obj$wn) != length(wn_reference) ||
             max(abs(obj$wn - wn_reference), na.rm = TRUE) > 1e-8) {
    stop("Raman-axis mismatch: ", row_i$mat_path)
  }

  pixel_total_abs <- rowSums(abs(obj$X))
  valid_total <- is.finite(pixel_total_abs) & pixel_total_abs > 0

  layer_qc_list[[length(layer_qc_list) + 1L]] <- data.table(
    sample_id = as.character(row_i$sample_id),
    group = as.character(row_i$group),
    layer = as.integer(row_i$layer),
    n_pixels_before = obj$n_pixels_before,
    n_pixels_after = obj$n_pixels_after,
    n_nonfinite_pixels = obj$n_nonfinite_pixels,
    n_wavenumbers = length(obj$wn),
    layer_mean_spectrum_intensity = mean(rowMeans(obj$X)),
    layer_median_spectrum_intensity = median(rowMeans(obj$X)),
    layer_mean_sum_abs_intensity = mean(pixel_total_abs),
    layer_median_sum_abs_intensity = median(pixel_total_abs),
    negative_fraction = mean(obj$X < 0),
    zero_fraction = mean(obj$X == 0)
  )

  for (f in seq_len(nrow(feature_2d_dt))) {
    ff <- feature_2d_dt[f]
    idx <- select_feature_indices(obj$wn, ff)
    if (length(idx) == 0L) {
      warning("No measured channels selected for feature ", ff$feature_label,
              " in ", row_i$mat_path)
      next
    }
    channels <- obj$wn[idx]

    axis_audit_list[[length(axis_audit_list) + 1L]] <- data.table(
      sample_id = as.character(row_i$sample_id),
      group = as.character(row_i$group),
      layer = as.integer(row_i$layer),
      feature_label = ff$feature_label,
      expected_feature_range = ff$feature_label,
      actual_channels = paste(sprintf("%.3f", channels), collapse = ";"),
      n_actual_channels = length(idx),
      actual_min = min(channels),
      actual_max = max(channels)
    )

    band_matrix <- obj$X[, idx, drop = FALSE]
    band_pixel_mean <- rowMeans(band_matrix)
    band_pixel_sum <- rowSums(band_matrix)

    if (length(idx) >= 2L) {
      dx <- diff(channels)
      weights <- numeric(length(idx))
      weights[1] <- dx[1] / 2
      weights[length(idx)] <- dx[length(dx)] / 2
      if (length(idx) > 2L) {
        weights[2:(length(idx) - 1L)] <- (head(dx, -1L) + tail(dx, -1L)) / 2
      }
      band_pixel_auc <- as.numeric(band_matrix %*% weights)
    } else {
      band_pixel_auc <- band_pixel_mean
    }

    relative_band <- rep(NA_real_, nrow(obj$X))
    relative_band[valid_total] <- band_pixel_mean[valid_total] / pixel_total_abs[valid_total]

    profile_list[[length(profile_list) + 1L]] <- data.table(
      sample_id = as.character(row_i$sample_id),
      group = as.character(row_i$group),
      layer = as.integer(row_i$layer),
      feature_label = ff$feature_label,
      feature_min = ff$feature_min,
      feature_max = ff$feature_max,
      expected_pc_label = ff$expected_pc_label,
      n_pixels = nrow(obj$X),
      n_channels = length(idx),
      actual_channels = paste(sprintf("%.3f", channels), collapse = ";"),
      band_mean = mean(band_pixel_mean),
      band_median = median(band_pixel_mean),
      band_sd_pixel = sd(band_pixel_mean),
      band_q25 = as.numeric(quantile(band_pixel_mean, 0.25, names = FALSE)),
      band_q75 = as.numeric(quantile(band_pixel_mean, 0.75, names = FALSE)),
      band_sum_mean = mean(band_pixel_sum),
      band_auc_mean = mean(band_pixel_auc),
      band_relative_total_abs_mean = mean(relative_band, na.rm = TRUE),
      band_relative_total_abs_median = median(relative_band, na.rm = TRUE),
      layer_mean_sum_abs_intensity = mean(pixel_total_abs)
    )

    for (j in seq_along(idx)) {
      x <- band_matrix[, j]
      single_channel_list[[length(single_channel_list) + 1L]] <- data.table(
        sample_id = as.character(row_i$sample_id),
        group = as.character(row_i$group),
        layer = as.integer(row_i$layer),
        feature_label = ff$feature_label,
        actual_wavenumber = channels[j],
        channel_mean = mean(x),
        channel_median = median(x),
        channel_sd_pixel = sd(x),
        channel_q25 = as.numeric(quantile(x, 0.25, names = FALSE)),
        channel_q75 = as.numeric(quantile(x, 0.75, names = FALSE))
      )
    }
  }
  rm(obj)
  gc()
}

depth_profile <- rbindlist(profile_list, fill = TRUE)
single_channel_profile <- rbindlist(single_channel_list, fill = TRUE)
layer_qc <- rbindlist(layer_qc_list, fill = TRUE)
axis_audit <- rbindlist(axis_audit_list, fill = TRUE)

for (dt in list(depth_profile, single_channel_profile, layer_qc, axis_audit)) {
  if ("group" %in% names(dt)) dt[, group := factor(group, levels = group_order)]
}

fwrite(depth_profile, file.path(out_direct, "FIXED_BAND_depth_profiles_all_layers.csv"), bom = TRUE)
fwrite(single_channel_profile, file.path(out_direct, "FIXED_BAND_single_channel_depth_profiles.csv"), bom = TRUE)
fwrite(layer_qc, file.path(out_qc, "depth_layer_total_intensity_QC.csv"), bom = TRUE)
fwrite(axis_audit, file.path(out_qc, "FIXED_BAND_actual_channel_selection_AUDIT.csv"), bom = TRUE)

# Dynamic near-zero tolerance estimated from Control depth variation.
control_tolerance <- depth_profile[group == control_group, .(
  control_layer_mad = mad(band_mean, constant = 1.4826),
  control_abs_median = median(abs(band_mean)),
  n_control_layers = uniqueN(layer)
), by = feature_label]
control_tolerance[, delta_tolerance := pmax(
  delta_tolerance_absolute_floor,
  delta_tolerance_mad_multiplier * control_layer_mad,
  delta_tolerance_relative_fraction * control_abs_median,
  na.rm = TRUE
)]
control_tolerance[!is.finite(delta_tolerance), delta_tolerance := delta_tolerance_absolute_floor]
fwrite(control_tolerance, file.path(out_direct, "FIXED_BAND_delta_near_zero_tolerance.csv"), bom = TRUE)

control_depth <- depth_profile[group == control_group, .(
  feature_label,
  layer,
  control_sample_id = sample_id,
  control_band_mean = band_mean,
  control_band_median = band_median,
  control_band_auc_mean = band_auc_mean,
  control_band_relative_total_abs_mean = band_relative_total_abs_mean,
  control_layer_mean_sum_abs_intensity = layer_mean_sum_abs_intensity
)]

depth_delta <- merge(
  depth_profile[group != control_group],
  control_depth,
  by = c("feature_label", "layer"),
  all.x = TRUE,
  sort = FALSE
)
depth_delta <- merge(depth_delta, control_tolerance, by = "feature_label", all.x = TRUE)
depth_delta[, `:=`(
  delta_mean = band_mean - control_band_mean,
  delta_median = band_median - control_band_median,
  delta_auc_mean = band_auc_mean - control_band_auc_mean,
  delta_relative_total_abs_mean =
    band_relative_total_abs_mean - control_band_relative_total_abs_mean,
  total_intensity_delta =
    layer_mean_sum_abs_intensity - control_layer_mean_sum_abs_intensity
)]
depth_delta[, direction_mean := fifelse(
  delta_mean > delta_tolerance, "up_vs_control",
  fifelse(delta_mean < -delta_tolerance, "down_vs_control", "near_zero")
)]
setorder(depth_delta, feature_label, group, layer)
fwrite(depth_delta, file.path(out_direct, "FIXED_BAND_Disease_minus_Control_layerwise_delta.csv"), bom = TRUE)

depth_direction_summary <- depth_delta[, {
  n_valid <- sum(is.finite(delta_mean))
  n_up <- sum(direction_mean == "up_vs_control", na.rm = TRUE)
  n_down <- sum(direction_mean == "down_vs_control", na.rm = TRUE)
  n_near <- sum(direction_mean == "near_zero", na.rm = TRUE)
  dominant <- if (n_up > n_down) "up_vs_control" else if (n_down > n_up) {
    "down_vs_control"
  } else "tie_or_near_zero"
  dominant_n <- max(n_up, n_down)
  data.table(
    feature_min = feature_min[1],
    feature_max = feature_max[1],
    expected_pc_label = expected_pc_label[1],
    n_layers_valid = n_valid,
    n_layers_up = n_up,
    n_layers_down = n_down,
    n_layers_near_zero = n_near,
    dominant_direction = dominant,
    dominant_direction_layers = dominant_n,
    direction_consistency_all_layers = if (n_valid > 0) dominant_n / n_valid else NA_real_,
    direction_consistency_nonzero_layers = if ((n_up + n_down) > 0) {
      dominant_n / (n_up + n_down)
    } else NA_real_,
    median_delta_across_layers = median(delta_mean, na.rm = TRUE),
    mean_delta_across_layers = mean(delta_mean, na.rm = TRUE),
    median_abs_delta_across_layers = median(abs(delta_mean), na.rm = TRUE),
    longest_up_run = longest_condition_run(layer, direction_mean == "up_vs_control"),
    longest_down_run = longest_condition_run(layer, direction_mean == "down_vs_control"),
    longest_dominant_run = if (dominant == "up_vs_control") {
      longest_condition_run(layer, direction_mean == "up_vs_control")
    } else if (dominant == "down_vs_control") {
      longest_condition_run(layer, direction_mean == "down_vs_control")
    } else 0L,
    delta_tolerance = delta_tolerance[1]
  )
}, by = .(feature_label, group)]
depth_direction_summary[, depth_support_tier := mapply(
  support_tier_from_fraction,
  dominant_direction_layers,
  n_layers_valid,
  USE.NAMES = FALSE
)]
fwrite(depth_direction_summary, file.path(out_direct, "FIXED_BAND_depth_direction_consistency_SUMMARY.csv"), bom = TRUE)

# Robust raw-signal depth stability within each placenta.
depth_stability <- depth_profile[, .(
  n_layers = uniqueN(layer),
  mean_across_layers = mean(band_mean),
  median_across_layers = median(band_mean),
  sd_across_layers = sd(band_mean),
  mad_across_layers = mad(band_mean, constant = 1.4826),
  mean_abs_across_layers = mean(abs(band_mean)),
  robust_cv_mad = mad(band_mean, constant = 1) / median(abs(band_mean)),
  range_over_mean_abs = (max(band_mean) - min(band_mean)) / mean(abs(band_mean))
), by = .(group, sample_id, feature_label)]
depth_stability[!is.finite(robust_cv_mad), robust_cv_mad := NA_real_]
depth_stability[!is.finite(range_over_mean_abs), range_over_mean_abs := NA_real_]
fwrite(depth_stability, file.path(out_direct, "FIXED_BAND_depth_stability_ROBUST_METRICS.csv"), bom = TRUE)

# Total-intensity correlation is a caution flag, not automatic invalidation.
delta_total_correlation <- data.table()
if (isTRUE(run_depth_qc)) {
  delta_total_correlation <- depth_delta[, .(
    spearman_delta_feature_vs_total_intensity_delta = suppressWarnings(cor(
      delta_mean,
      total_intensity_delta,
      method = "spearman",
      use = "complete.obs"
    ))
  ), by = .(feature_label, group)]
  delta_total_correlation[, QC_caution :=
    abs(spearman_delta_feature_vs_total_intensity_delta) >=
      qc_correlation_caution_threshold]
  fwrite(delta_total_correlation,
         file.path(out_qc, "FIXED_BAND_delta_vs_total_intensity_QC.csv"), bom = TRUE)
}


# 4. Balanced all-common-layer fixed reference PCA for each disease

reference_models <- list()
reference_loading_list <- list()
reference_peak_list <- list()  # formal PROMINENT signed-loading extrema
reference_sep_list <- list()
reference_projection_list <- list()
reference_target_point_list <- list()
reference_target_layer_list <- list()
reference_to_2d_match_list <- list()

for (dg in available_diseases) {
  common_layers <- available_pairs[disease_group == dg, layer]
  if (length(common_layers) == 0L) next
  message("\nFixed reference PCA: Control vs ", dg)

  Xc_list <- list()
  Xd_list <- list()
  wn_ref <- NULL
  for (lz in common_layers) {
    c_row <- control_meta[layer == lz][1]
    d_row <- meta[group == dg & layer == lz][1]
    c_obj <- read_layer(c_row)
    d_obj <- read_layer(d_row)
    if (length(c_obj$wn) != length(d_obj$wn) ||
        max(abs(c_obj$wn - d_obj$wn)) > 1e-8) {
      stop("Raman-axis mismatch in reference PCA: Control vs ", dg, " layer ", lz)
    }
    if (is.null(wn_ref)) wn_ref <- c_obj$wn
    # Exact per-layer balance: both groups contribute the same number of pixels.
    n_take <- min(
      nrow(c_obj$X),
      nrow(d_obj$X),
      as.integer(reference_max_pixels_per_layer_per_group)
    )
    Xc_list[[as.character(lz)]] <- c_obj$X[
      sample.int(nrow(c_obj$X), n_take, replace = FALSE),
      , drop = FALSE
    ]
    Xd_list[[as.character(lz)]] <- d_obj$X[
      sample.int(nrow(d_obj$X), n_take, replace = FALSE),
      , drop = FALSE
    ]
    rm(c_obj, d_obj)
  }

  Xc_ref <- do.call(rbind, Xc_list)
  Xd_ref <- do.call(rbind, Xd_list)
  pair_id <- paste0("Control_vs_", dg, "_fixed_depth_reference")
  model <- fit_pairwise_pca(Xc_ref, Xd_ref, wn_ref, pair_id)
  reference_models[[dg]] <- model

  ref_dir <- file.path(out_reference, paste0("Control_vs_", dg))
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(model, file.path(ref_dir, "fixed_reference_PCA_model.rds"), compress = TRUE)
  fwrite(model$separation, file.path(ref_dir, "fixed_reference_PC_separation.csv"), bom = TRUE)
  fwrite(model$diagnostics, file.path(ref_dir, "fixed_reference_PCA_numerical_diagnostics.csv"), bom = TRUE)

  loading_dt <- loading_table_from_model(
    model,
    disease_group = dg,
    layer = NA_integer_,
    source = "fixed_depth_reference_PCA"
  )
  loading_dt[, layer := NA_integer_]
  fwrite(loading_dt, file.path(ref_dir, "fixed_reference_PC1_PC3_loadings.csv"), bom = TRUE)
  reference_loading_list[[dg]] <- loading_dt

  # Explicitly match this disease's 3D all-depth reference PCs to its existing
  # 2D placenta PC axes. PC numbers are not assumed to correspond a priori.
  ref_to_2d <- match_reference_pcs_to_2d(
    loading_dt,
    loading_2d_list[[dg]]
  )
  ref_to_2d[, disease_group := dg]
  fwrite(
    ref_to_2d,
    file.path(ref_dir, "fixed_reference_PC_match_to_2D_placenta.csv"),
    bom = TRUE
  )
  reference_to_2d_match_list[[dg]] <- ref_to_2d

  delta_ref <- colMeans(Xd_ref) - colMeans(Xc_ref)
  delta_ref <- delta_ref[model$keep_cols]
  peaks_ref <- extract_all_prominent_extrema(loading_dt, delta_ref)
  peaks_ref <- merge(
    peaks_ref,
    ref_to_2d[, .(
      PC = reference_PC,
      matched_2D_PC,
      signed_loading_correlation_3Dref_to_2D,
      absolute_loading_correlation_3Dref_to_2D,
      match_quality_3Dref_to_2D
    )],
    by = "PC",
    all.x = TRUE
  )
  fwrite(peaks_ref, file.path(ref_dir, "fixed_reference_PROMINENT_loading_extrema.csv"), bom = TRUE)
  reference_peak_list[[dg]] <- peaks_ref
  reference_sep_list[[dg]] <- model$separation[, disease_group := dg]

  target_points <- rbindlist(lapply(c("strict", "buffered"), function(w) {
    summarize_target_support_one_loading(
      loading_dt,
      feature_2d_dt,
      reference_to_2d_mapping = ref_to_2d,
      feature_direction_dt = twoD_feature_direction_dt,
      pc_mapping = NULL,
      window_type = w
    )
  }), fill = TRUE)
  target_layers <- collapse_target_support_by_layer(target_points)
  fwrite(target_points, file.path(ref_dir, "fixed_reference_TARGET_feature_support_points.csv"), bom = TRUE)
  fwrite(target_layers, file.path(ref_dir, "fixed_reference_TARGET_feature_support_summary.csv"), bom = TRUE)
  reference_target_point_list[[dg]] <- target_points
  reference_target_layer_list[[dg]] <- target_layers

  # Project every layer into this one fixed coordinate system.
  for (lz in common_layers) {
    c_row <- control_meta[layer == lz][1]
    d_row <- meta[group == dg & layer == lz][1]
    c_obj <- read_layer(c_row)
    d_obj <- read_layer(d_row)
    sc <- project_to_reference(c_obj$X, model)
    sd <- project_to_reference(d_obj$X, model)
    for (k in seq_len(n_pcs)) {
      reference_projection_list[[length(reference_projection_list) + 1L]] <- data.table(
        disease_group = dg,
        layer = lz,
        reference_PC = paste0("PC", k),
        matched_2D_PC = ref_to_2d[reference_PC == paste0("PC", k), matched_2D_PC][1],
        signed_loading_correlation_3Dref_to_2D = ref_to_2d[
          reference_PC == paste0("PC", k),
          signed_loading_correlation_3Dref_to_2D
        ][1],
        variance_explained_percent = 100 * model$variance_explained[k],
        control_mean_score = mean(sc[, k]),
        disease_mean_score = mean(sd[, k]),
        disease_minus_control_score = mean(sd[, k]) - mean(sc[, k]),
        cohens_d_disease_vs_control = cohens_d_two_group(sc[, k], sd[, k]),
        n_control_pixels = nrow(sc),
        n_disease_pixels = nrow(sd)
      )
    }
    rm(c_obj, d_obj, sc, sd)
  }
  rm(Xc_ref, Xd_ref, Xc_list, Xd_list)
  gc()
}

reference_loadings_all <- rbindlist(reference_loading_list, fill = TRUE)
reference_peaks_all <- rbindlist(reference_peak_list, fill = TRUE)
reference_separation_all <- rbindlist(reference_sep_list, fill = TRUE)
reference_projection <- rbindlist(reference_projection_list, fill = TRUE)
reference_target_points <- rbindlist(reference_target_point_list, fill = TRUE)
reference_target_summary <- rbindlist(reference_target_layer_list, fill = TRUE)
reference_to_2d_matches_all <- rbindlist(reference_to_2d_match_list, fill = TRUE)

fwrite(reference_to_2d_matches_all, file.path(out_reference, "ALL_fixed_reference_PC_matches_to_2D_placenta.csv"), bom = TRUE)
fwrite(reference_loadings_all, file.path(out_reference, "ALL_fixed_reference_PC_loadings.csv"), bom = TRUE)
fwrite(reference_peaks_all, file.path(out_reference, "ALL_fixed_reference_PROMINENT_loading_extrema.csv"), bom = TRUE)
fwrite(reference_separation_all, file.path(out_reference, "ALL_fixed_reference_PC_separation.csv"), bom = TRUE)
fwrite(reference_projection, file.path(out_reference, "FIXED_REFERENCE_PCA_layerwise_score_centroid_profiles.csv"), bom = TRUE)
fwrite(reference_target_points, file.path(out_reference, "ALL_fixed_reference_TARGET_feature_support_points.csv"), bom = TRUE)
fwrite(reference_target_summary, file.path(out_reference, "ALL_fixed_reference_TARGET_feature_support_summary.csv"), bom = TRUE)


# 5. Independent layer-wise PCA, reference-PC matching, and PROMINENT-extremum recurrence

layer_loading_list <- list()
layer_peak_list <- list()      # formal PROMINENT signed-loading extrema
layer_sep_list <- list()
pc_match_list <- list()
target_point_list <- list()
target_layer_list <- list()

for (dg in available_diseases) {
  model_ref <- reference_models[[dg]]
  ref_loading <- reference_loading_list[[dg]]
  ref_to_2d <- reference_to_2d_match_list[[dg]]
  common_layers <- available_pairs[disease_group == dg, layer]
  disease_dir <- file.path(out_layerwise, paste0("Control_vs_", dg))
  dir.create(disease_dir, recursive = TRUE, showWarnings = FALSE)

  for (lz in common_layers) {
    message("Layer-wise PCA: Control vs ", dg, " layer ", lz)
    c_row <- control_meta[layer == lz][1]
    d_row <- meta[group == dg & layer == lz][1]
    c_obj <- read_layer(c_row)
    d_obj <- read_layer(d_row)

    pair_id <- paste0("Control_vs_", dg, "_layer", sprintf("%02d", lz))
    model <- fit_pairwise_pca(c_obj$X, d_obj$X, c_obj$wn, pair_id)
    layer_dir <- file.path(disease_dir, paste0("layer", sprintf("%02d", lz)))
    dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)

    loading_dt <- loading_table_from_model(
      model,
      disease_group = dg,
      layer = lz,
      source = "independent_layerwise_PCA"
    )
    pc_map <- match_layer_pcs_to_reference(loading_dt, ref_loading)
    pc_map <- merge(
      pc_map,
      ref_to_2d[, .(
        reference_PC,
        matched_2D_PC,
        signed_loading_correlation_3Dref_to_2D,
        absolute_loading_correlation_3Dref_to_2D,
        match_quality_3Dref_to_2D
      )],
      by = "reference_PC",
      all.x = TRUE,
      sort = FALSE
    )
    pc_map[, `:=`(pair_id = pair_id, disease_group = dg, layer = lz)]

    delta_layer <- colMeans(d_obj$X) - colMeans(c_obj$X)
    delta_layer <- delta_layer[model$keep_cols]
    peaks <- extract_all_prominent_extrema(loading_dt, delta_layer)
    peaks <- merge(
      peaks,
      pc_map[, .(
        PC = layer_PC,
        matched_reference_PC = reference_PC,
        matched_2D_PC = matched_2D_PC,
        signed_loading_correlation_to_reference = signed_loading_correlation,
        absolute_loading_correlation_to_reference = absolute_loading_correlation,
        reference_PC_match_quality = match_quality,
        signed_loading_correlation_reference_to_2D = signed_loading_correlation_3Dref_to_2D,
        absolute_loading_correlation_reference_to_2D = absolute_loading_correlation_3Dref_to_2D,
        reference_to_2D_match_quality = match_quality_3Dref_to_2D
      )],
      by = "PC",
      all.x = TRUE
    )

    target_points <- rbindlist(lapply(c("strict", "buffered"), function(w) {
      summarize_target_support_one_loading(
        loading_dt,
        feature_2d_dt,
        reference_to_2d_mapping = ref_to_2d,
        feature_direction_dt = twoD_feature_direction_dt,
        pc_mapping = pc_map,
        window_type = w
      )
    }), fill = TRUE)
    target_layers <- collapse_target_support_by_layer(target_points)

    fwrite(model$separation, file.path(layer_dir, "PC_separation.csv"), bom = TRUE)
    fwrite(model$diagnostics, file.path(layer_dir, "PCA_numerical_diagnostics.csv"), bom = TRUE)
    fwrite(loading_dt, file.path(layer_dir, "PC1_PC3_loadings_contributions.csv"), bom = TRUE)
    fwrite(pc_map, file.path(layer_dir, "PC_axis_match_to_fixed_reference.csv"), bom = TRUE)
    fwrite(peaks, file.path(layer_dir, "PROMINENT_loading_extrema.csv"), bom = TRUE)
    fwrite(target_points, file.path(layer_dir, "TARGET_feature_support_points.csv"), bom = TRUE)
    fwrite(target_layers, file.path(layer_dir, "TARGET_feature_support_layer_summary.csv"), bom = TRUE)

    if (isTRUE(save_layerwise_loading_plots)) {
      plot_long <- rbindlist(lapply(seq_len(n_pcs), function(k) {
        pc <- paste0("PC", k)
        data.table(
          wavenumber = loading_dt$wavenumber,
          PC = pc,
          loading = loading_dt[[paste0(pc, "_loading_disease_aligned")]]
        )
      }))
      p <- ggplot(plot_long, aes(wavenumber, loading)) +
        geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
        geom_line(linewidth = 0.45) +
        facet_wrap(~PC, ncol = 1, scales = "free_y") +
        coord_cartesian(xlim = c(analysis_min, analysis_max)) +
        theme_classic(base_size = 10) +
        labs(
          title = pair_id,
          x = expression("Wavenumber (cm"^{-1}*")"),
          y = "Disease-aligned loading"
        )
      ggsave(file.path(layer_dir, "PC1_PC3_loading_spectra.png"), p,
             width = 7.5, height = 7.5, dpi = plot_dpi, bg = "white")
    }

    layer_loading_list[[pair_id]] <- loading_dt
    layer_peak_list[[pair_id]] <- peaks
    layer_sep_list[[pair_id]] <- model$separation[, `:=`(disease_group = dg, layer = lz)]
    pc_match_list[[pair_id]] <- pc_map
    target_point_list[[pair_id]] <- target_points
    target_layer_list[[pair_id]] <- target_layers

    rm(c_obj, d_obj, model)
    gc()
  }
}

layer_loadings_all <- rbindlist(layer_loading_list, fill = TRUE)
layer_peaks_all <- rbindlist(layer_peak_list, fill = TRUE)
layer_separation_all <- rbindlist(layer_sep_list, fill = TRUE)
pc_matches_all <- rbindlist(pc_match_list, fill = TRUE)
target_points_all <- rbindlist(target_point_list, fill = TRUE)
target_layers_all <- rbindlist(target_layer_list, fill = TRUE)
target_summary_all <- summarize_target_support_across_layers(target_layers_all)

# Defensive schema audit: ggplot2 requires uniquely named columns.
if (anyDuplicated(names(target_layers_all))) {
  stop(
    "target_layers_all has duplicate column names: ",
    paste(names(target_layers_all)[duplicated(names(target_layers_all))], collapse = ", ")
  )
}
if (anyDuplicated(names(target_summary_all))) {
  stop(
    "target_summary_all has duplicate column names: ",
    paste(names(target_summary_all)[duplicated(names(target_summary_all))], collapse = ", ")
  )
}

fwrite(layer_loadings_all, file.path(out_layerwise, "ALL_layerwise_PC_loadings_contributions.csv"), bom = TRUE)
fwrite(layer_peaks_all, file.path(out_layerwise, "ALL_layerwise_PROMINENT_loading_extrema.csv"), bom = TRUE)
fwrite(layer_separation_all, file.path(out_layerwise, "ALL_layerwise_PC_separation.csv"), bom = TRUE)
fwrite(pc_matches_all, file.path(out_layerwise, "ALL_layerwise_PC_axis_matches_to_fixed_reference.csv"), bom = TRUE)
fwrite(target_points_all, file.path(out_layerwise, "ALL_TARGET_feature_support_points.csv"), bom = TRUE)
fwrite(target_layers_all, file.path(out_layerwise, "ALL_TARGET_feature_support_layerwise.csv"), bom = TRUE)
fwrite(target_summary_all, file.path(out_layerwise, "ALL_TARGET_feature_support_across_layers_SUMMARY.csv"), bom = TRUE)


# 6. De novo layer-wise PROMINENT-extremum consensus regions and PC-specific 2D overlap

consensus_regions <- merge_consensus_regions(layer_peaks_all)
fwrite(consensus_regions, file.path(out_consensus, "ALL_3D_layerwise_consensus_feature_regions_RANKED.csv"), bom = TRUE)

# Explicit target overlap summary.
target_consensus_overlap <- rbindlist(lapply(seq_len(nrow(feature_2d_dt)), function(f) {
  ff <- feature_2d_dt[f]
  hits <- consensus_regions[
    matched_2D_PC == ff$expected_pc_label &
      region_max >= ff$feature_min - feature_overlap_tolerance_cm1 &
      region_min <= ff$feature_max + feature_overlap_tolerance_cm1
  ]
  if (nrow(hits) == 0L) {
    return(data.table(
      feature_label = ff$feature_label,
      n_consensus_regions_overlapping = 0L,
      best_consensus_region = NA_character_,
      best_consensus_score = NA_real_,
      n_diseases_in_best_region = 0L,
      n_disease_layer_pairs_in_best_region = 0L
    ))
  }
  best <- hits[which.max(consensus_score)]
  data.table(
    feature_label = ff$feature_label,
    n_consensus_regions_overlapping = nrow(hits),
    best_consensus_region = paste0(round(best$region_min, 1), "-", round(best$region_max, 1)),
    best_consensus_score = best$consensus_score,
    n_diseases_in_best_region = best$n_diseases,
    n_disease_layer_pairs_in_best_region = best$n_disease_layer_pairs
  )
}))
fwrite(target_consensus_overlap, file.path(out_consensus, "TARGET_2D_features_overlap_with_3D_consensus.csv"), bom = TRUE)


# 7. Integrated feature-level evidence summary

direct_feature_summary <- depth_direction_summary[, .(
  n_diseases_direct_strong = sum(depth_support_tier == "strong_depth_support", na.rm = TRUE),
  n_diseases_direct_moderate_or_strong = sum(
    depth_support_tier %in% c("strong_depth_support", "moderate_depth_support"),
    na.rm = TRUE
  ),
  median_direct_direction_consistency = median(direction_consistency_all_layers, na.rm = TRUE),
  median_abs_Disease_Control_delta = median(abs(median_delta_across_layers), na.rm = TRUE),
  major_direct_direction = {
    tt <- table(dominant_direction)
    if (length(tt) == 0L) NA_character_ else names(tt)[which.max(tt)]
  },
  n_diseases_major_direct_direction = {
    tt <- table(dominant_direction)
    if (length(tt) == 0L) 0L else as.integer(max(tt))
  }
), by = feature_label]

# Layer-wise PCA target support across diseases. PRIMARY inference uses the STRICT
# target interval and the PC identity matched to the fixed reference PCA.
pca_feature_summary <- target_summary_all[window_type == "strict", .(
  n_diseases_anyPC_strong = sum(any_PC_support_tier == "strong_depth_support", na.rm = TRUE),
  n_diseases_anyPC_moderate_or_strong = sum(
    any_PC_support_tier %in% c("strong_depth_support", "moderate_depth_support"),
    na.rm = TRUE
  ),
  n_diseases_expected_refPC_strong = sum(
    expected_reference_PC_support_tier == "strong_depth_support", na.rm = TRUE
  ),
  n_diseases_expected_refPC_moderate_or_strong = sum(
    expected_reference_PC_support_tier %in%
      c("strong_depth_support", "moderate_depth_support"),
    na.rm = TRUE
  ),
  median_anyPC_layer_fraction = median(any_PC_layer_fraction, na.rm = TRUE),
  median_expected_refPC_layer_fraction = median(
    expected_reference_PC_layer_fraction, na.rm = TRUE
  ),
  median_layerwise_max_relative_abs_loading = median(
    median_maximum_relative_abs_loading, na.rm = TRUE
  ),
  n_diseases_direction_evaluable = sum(loading_direction_layers_evaluable > 0, na.rm = TRUE),
  n_diseases_direction_concordance_ge80 = sum(
    loading_direction_concordance_fraction_among_evaluable >= 0.80,
    na.rm = TRUE
  ),
  n_diseases_direction_concordance_ge60 = sum(
    loading_direction_concordance_fraction_among_evaluable >= 0.60,
    na.rm = TRUE
  ),
  median_loading_direction_concordance_fraction = safe_median_numeric(
    loading_direction_concordance_fraction_among_evaluable
  ),
  median_loading_direction_concordance_fraction_high_confidence = safe_median_numeric(
    loading_direction_concordance_fraction_high_confidence
  )
), by = feature_label]

# Fixed reference PCA support: one summary per disease and feature.
# PRIMARY inference uses the STRICT feature interval.
reference_feature_summary <- reference_target_summary[window_type == "strict", .(
  n_diseases_fixed_reference_anyPC_supported = sum(any_PC_support, na.rm = TRUE),
  n_diseases_fixed_reference_expectedPC_supported = sum(
    expected_reference_PC_support, na.rm = TRUE
  ),
  median_fixed_reference_max_relative_abs_loading = median(
    maximum_relative_abs_loading, na.rm = TRUE
  ),
  n_diseases_fixed_reference_direction_evaluable = sum(
    expected_reference_PC_direction_evaluable, na.rm = TRUE
  ),
  n_diseases_fixed_reference_direction_concordant = sum(
    expected_reference_PC_direction_concordant, na.rm = TRUE
  ),
  fixed_reference_direction_concordance_fraction = if (
    sum(expected_reference_PC_direction_evaluable, na.rm = TRUE) > 0
  ) {
    sum(expected_reference_PC_direction_concordant, na.rm = TRUE) /
      sum(expected_reference_PC_direction_evaluable, na.rm = TRUE)
  } else NA_real_
), by = feature_label]

qc_feature_summary <- data.table(feature_label = feature_2d_dt$feature_label)
qc_feature_summary[, `:=`(n_disease_QC_cautions = 0L, QC_caution_groups = "")]
if (nrow(delta_total_correlation) > 0L) {
  tmp <- delta_total_correlation[, .(
    n_disease_QC_cautions = sum(QC_caution, na.rm = TRUE),
    QC_caution_groups = paste(as.character(group[QC_caution]), collapse = ";")
  ), by = feature_label]
  qc_feature_summary <- merge(
    qc_feature_summary[, .(feature_label)],
    tmp,
    by = "feature_label",
    all.x = TRUE
  )
  qc_feature_summary[is.na(n_disease_QC_cautions), n_disease_QC_cautions := 0L]
  qc_feature_summary[is.na(QC_caution_groups), QC_caution_groups := ""]
}

integrated <- merge(feature_export, direct_feature_summary, by = "feature_label", all.x = TRUE)
integrated <- merge(integrated, pca_feature_summary, by = "feature_label", all.x = TRUE)
integrated <- merge(integrated, reference_feature_summary, by = "feature_label", all.x = TRUE)
integrated <- merge(integrated, target_consensus_overlap, by = "feature_label", all.x = TRUE)
integrated <- merge(integrated, qc_feature_summary, by = "feature_label", all.x = TRUE)

# Count fields are zero when a feature has no support in a given evidence module.
count_columns <- c(
  "n_diseases_direct_strong",
  "n_diseases_direct_moderate_or_strong",
  "n_diseases_major_direct_direction",
  "n_diseases_anyPC_strong",
  "n_diseases_anyPC_moderate_or_strong",
  "n_diseases_expected_refPC_strong",
  "n_diseases_expected_refPC_moderate_or_strong",
  "n_diseases_fixed_reference_anyPC_supported",
  "n_diseases_fixed_reference_expectedPC_supported",
  "n_diseases_direction_evaluable",
  "n_diseases_direction_concordance_ge80",
  "n_diseases_direction_concordance_ge60",
  "n_diseases_fixed_reference_direction_evaluable",
  "n_diseases_fixed_reference_direction_concordant",
  "n_consensus_regions_overlapping",
  "n_diseases_in_best_region",
  "n_disease_layer_pairs_in_best_region",
  "n_disease_QC_cautions"
)
for (nm in intersect(count_columns, names(integrated))) {
  set(integrated, which(is.na(integrated[[nm]])), nm, 0)
}

n_available_diseases <- length(available_diseases)
majority_cutoff <- ceiling(n_available_diseases * 2 / 3)
half_cutoff <- ceiling(n_available_diseases / 2)

integrated[, final_3D_reproducibility_tier := fifelse(
  n_diseases_direct_moderate_or_strong >= majority_cutoff &
    n_diseases_expected_refPC_moderate_or_strong >= majority_cutoff &
    n_diseases_fixed_reference_expectedPC_supported >= majority_cutoff &
    n_disease_QC_cautions <= 2,
  "Tier 1: high 3D depth reproducibility with PC-specific PCA support",
  fifelse(
    n_diseases_direct_moderate_or_strong >= half_cutoff &
      n_diseases_expected_refPC_moderate_or_strong >= half_cutoff,
    "Tier 2: moderate or disease-heterogeneous 3D reproducibility",
    "Tier 3: weak/inconsistent 3D reproducibility"
  )
)]

setorder(integrated, feature_min)
fwrite(integrated, file.path(out_integrated, "INTEGRATED_2D_to_3D_feature_reproducibility_summary.csv"), bom = TRUE)


# 8. Figures

if (nrow(depth_profile) > 0L) {
  p <- ggplot(depth_profile, aes(layer, band_mean, color = group, group = group)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3) +
    facet_wrap(~feature_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = plot_group_colors, drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(depth_profile$layer))) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 7),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Fixed 2D Raman feature regions across 3D depth",
      x = "Z-stack layer",
      y = "Band mean intensity",
      color = NULL
    )
  ggsave(file.path(out_figures, "Fig3D_1_FIXED_BAND_depth_profiles.png"), p,
         width = 9.5, height = 8.2, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_1_FIXED_BAND_depth_profiles.pdf"), p,
         width = 9.5, height = 8.2, bg = "white")
}

if (nrow(depth_delta) > 0L) {
  p <- ggplot(depth_delta, aes(layer, delta_mean, color = group, group = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_ribbon(
      aes(ymin = -delta_tolerance, ymax = delta_tolerance),
      inherit.aes = TRUE,
      color = NA,
      fill = "grey90",
      alpha = 0.45
    ) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3) +
    facet_wrap(~feature_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = plot_group_colors[available_diseases], drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(depth_delta$layer))) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 7),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Disease-Control fixed-band depth deltas",
      subtitle = "Grey zone denotes the feature-specific near-zero tolerance",
      x = "Z-stack layer",
      y = "Disease - Control band mean",
      color = NULL
    )
  ggsave(file.path(out_figures, "Fig3D_2_FIXED_BAND_Disease_Control_delta.png"), p,
         width = 9.5, height = 8.2, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_2_FIXED_BAND_Disease_Control_delta.pdf"), p,
         width = 9.5, height = 8.2, bg = "white")
}

if (nrow(depth_direction_summary) > 0L) {
  plot_dt <- copy(depth_direction_summary)
  plot_dt[, label := paste0(n_layers_up, "↑/", n_layers_down, "↓/", n_layers_near_zero, "≈")]
  plot_dt[, feature_label := factor(feature_label, levels = rev(feature_2d_dt$feature_label))]
  plot_dt[, group := factor(group, levels = available_diseases)]
  p <- ggplot(plot_dt, aes(group, feature_label, fill = direction_consistency_all_layers)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient(low = "white", high = "#D7191C", limits = c(0, 1)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) +
    labs(
      title = "Fixed-band depth direction consistency",
      x = NULL,
      y = NULL,
      fill = "Consistency"
    )
  ggsave(file.path(out_figures, "Fig3D_3_FIXED_BAND_direction_consistency.png"), p,
         width = 8.5, height = 4.8, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_3_FIXED_BAND_direction_consistency.pdf"), p,
         width = 8.5, height = 4.8, bg = "white")
}

if (nrow(target_summary_all) > 0L) {
  plot_dt <- target_summary_all[window_type == "strict"]
  plot_dt[, feature_label := factor(feature_label, levels = rev(feature_2d_dt$feature_label))]
  plot_dt[, disease_group := factor(disease_group, levels = available_diseases)]
  p <- ggplot(plot_dt, aes(disease_group, feature_label, fill = any_PC_layer_fraction)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = paste0(any_PC_layers_supported, "/", n_layers_valid)), size = 3) +
    scale_fill_gradient(low = "white", high = "#D7191C", limits = c(0, 1)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) +
    labs(
      title = "Layer-wise prominent-loading support in any PC1-PC3",
      subtitle = "Strict 2D feature intervals; PROMINENT signed-loading extrema",
      x = NULL,
      y = NULL,
      fill = "Layer fraction"
    )
  ggsave(file.path(out_figures, "Fig3D_4_LAYERWISE_PCA_anyPC_target_support.png"), p,
         width = 8.5, height = 4.8, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_4_LAYERWISE_PCA_anyPC_target_support.pdf"), p,
         width = 8.5, height = 4.8, bg = "white")

  p2 <- ggplot(plot_dt, aes(disease_group, feature_label,
                            fill = expected_reference_PC_layer_fraction)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = paste0(expected_reference_PC_layers_supported, "/", n_layers_valid)), size = 3) +
    scale_fill_gradient(low = "white", high = "#2C7BB6", limits = c(0, 1)) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) +
    labs(
      title = "Reference-PC-matched prominent-loading support",
      subtitle = "Strict 2D intervals; layer PCs are matched one-to-one to the fixed reference PCA",
      x = NULL,
      y = NULL,
      fill = "Layer fraction"
    )
  ggsave(file.path(out_figures, "Fig3D_5_LAYERWISE_PCA_reference_matched_support.png"), p2,
         width = 8.5, height = 4.8, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_5_LAYERWISE_PCA_reference_matched_support.pdf"), p2,
         width = 8.5, height = 4.8, bg = "white")
}

if (nrow(reference_projection) > 0L) {
  p <- ggplot(reference_projection,
              aes(layer, disease_minus_control_score, color = disease_group,
                  group = disease_group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    facet_wrap(~reference_PC, scales = "free_y", ncol = 1) +
    scale_color_manual(values = plot_group_colors[available_diseases], drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(reference_projection$layer))) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 7),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = "Layer projection into fixed depth-reference PCA coordinates",
      x = "Z-stack layer",
      y = "Disease - Control score centroid",
      color = NULL
    )
  ggsave(file.path(out_figures, "Fig3D_6_FIXED_REFERENCE_score_centroid_depth_profiles.png"), p,
         width = 8.5, height = 8.5, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_6_FIXED_REFERENCE_score_centroid_depth_profiles.pdf"), p,
         width = 8.5, height = 8.5, bg = "white")
}

# Loading-direction concordance is displayed separately from recurrence.
if (nrow(target_summary_all) > 0L) {
  dir_dt <- target_summary_all[window_type == "strict"]
  dir_dt[, feature_label := factor(feature_label, levels = rev(feature_2d_dt$feature_label))]
  dir_dt[, disease_group := factor(disease_group, levels = available_diseases)]

  p_dir <- ggplot(
    dir_dt,
    aes(
      x = disease_group,
      y = feature_label,
      fill = loading_direction_concordance_fraction_high_confidence
    )
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = paste0(
        loading_direction_layers_concordant_high_confidence,
        "/",
        loading_direction_layers_evaluable_high_confidence
      )),
      size = 3
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#2C7BB6",
      limits = c(0, 1),
      na.value = "grey90"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) +
    labs(
      title = "2D-to-3D loading-direction concordance",
      subtitle = "Strict target interval; denominator = layers with expected-PC support, unidirectional formal 2D/3D loading patterns, and |3D-ref vs 2D PC loading correlation| >= 0.60",
      x = NULL,
      y = NULL,
      fill = "Direction\nconcordance"
    )

  ggsave(
    file.path(out_figures, "Fig3D_7_LOADING_direction_concordance.png"),
    p_dir, width = 8.5, height = 4.8, dpi = plot_dpi, bg = "white"
  )
  ggsave(
    file.path(out_figures, "Fig3D_7_LOADING_direction_concordance.pdf"),
    p_dir, width = 8.5, height = 4.8, bg = "white"
  )
}

if (nrow(consensus_regions) > 0L) {
  plot_dt <- head(consensus_regions, 30)

  # IMPORTANT:
  # Consensus regions are clustered separately by the PC identity matched to the 2D placenta PCA.
  # Therefore the same rounded Raman-region text can legitimately occur in
  # more than one reference PC (e.g. a ~1155 cm^-1 region in both PC1 and PC2).
  # Using the bare Raman label itself as factor levels can therefore create
  # duplicated levels and trigger: "factor level [...] is duplicated".
  plot_dt[, region_label := ifelse(
    region_width < 0.5,
    sprintf("%.1f", region_center),
    paste0(round(region_min, 1), "-", round(region_max, 1))
  )]

  # Make the plotting label PC-specific. This is also scientifically clearer
  # because 3D consensus regions are defined within 2D-matched PC identities.
  plot_dt[, plot_label := paste0(matched_2D_PC, " | ", region_label)]

  # Extra safeguard: if two distinct region_ids still round to the same text
  # within one PC, append region_id so every plotted row remains unique.
  dup_label <- duplicated(plot_dt$plot_label) | duplicated(plot_dt$plot_label, fromLast = TRUE)
  plot_dt[dup_label, plot_label := paste0(plot_label, " | ", region_id)]

  # factor levels MUST be unique.
  plot_dt[, plot_label := factor(
    plot_label,
    levels = rev(unique(as.character(plot_label)))
  )]

  p <- ggplot(
    plot_dt,
    aes(consensus_score, plot_label, fill = overlaps_2D_feature)
  ) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey65")) +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 8), legend.position = "bottom") +
    labs(
      title = "Layer-wise PROMINENT-loading depth-consensus Raman regions",
      x = "Consensus score",
      y = expression("Matched 2D PC | Raman region (cm"^{-1}*")"),
      fill = "Overlaps predefined 2D feature"
    )

  ggsave(file.path(out_figures, "Fig3D_8_LAYERWISE_consensus_feature_regions.png"), p,
         width = 8.5, height = 7, dpi = plot_dpi, bg = "white")
  ggsave(file.path(out_figures, "Fig3D_8_LAYERWISE_consensus_feature_regions.pdf"), p,
         width = 8.5, height = 7, bg = "white")
}


# 9. Method audit and interpretation notes

method_audit <- data.table(
  parameter = c(
    "analysis_range_cm1",
    "pca_method",
    "pca_center",
    "pca_scale",
    "extra_baseline_correction",
    "extra_smoothing",
    "extra_normalization",
    "nonfinite_action",
    "loading_peak_sg_window",
    "peak_min_fraction_absmax",
    "peak_min_distance_cm1",
    "merge_peak_distance_cm1",
    "consensus_max_region_width_cm1",
    "feature_overlap_tolerance_cm1",
    "band_edge_channel_tolerance_cm1",
    "reference_max_pixels_per_layer_per_group",
    "strong_support_fraction",
    "moderate_support_fraction",
    "delta_tolerance_mad_multiplier",
    "delta_tolerance_relative_fraction",
    "placenta2d_result_root",
    "cross_dimensional_PC_matching",
    "cross_dimensional_pc_match_min_abs_correlation",
    "loading_direction_concordance_rule",
    "n_available_diseases"
  ),
  value = c(
    paste0(analysis_min, "-", analysis_max),
    "covariance eigendecomposition equivalent to centered unscaled prcomp/sklearn PCA",
    as.character(pca_center),
    as.character(pca_scale),
    "FALSE",
    "FALSE",
    "FALSE",
    nonfinite_action,
    as.character(loading_peak_sg_window),
    as.character(peak_min_fraction_absmax),
    as.character(peak_min_distance_cm1),
    as.character(merge_peak_distance_cm1),
    as.character(consensus_max_region_width_cm1),
    as.character(feature_overlap_tolerance_cm1),
    as.character(band_edge_channel_tolerance_cm1),
    as.character(reference_max_pixels_per_layer_per_group),
    as.character(strong_support_fraction),
    as.character(moderate_support_fraction),
    as.character(delta_tolerance_mad_multiplier),
    as.character(delta_tolerance_relative_fraction),
    placenta2d_result_root,
    "one-to-one absolute loading-correlation assignment: 3D fixed reference PC1-PC3 <-> existing 2D PC1-PC3; layer PCs first matched to 3D reference PCs",
    as.character(cross_dimensional_pc_match_min_abs_correlation),
    "primary direction concordance requires a unidirectional formal 2D prominent-extremum pattern and a unidirectional 3D prominent-extremum pattern of the same disease-aligned sign within the strict 2D interval; mixed-sign windows are not forced into one direction; high-confidence metric additionally requires the cross-dimensional PC match threshold",
    as.character(n_available_diseases)
  )
)
fwrite(method_audit, file.path(out_qc, "METHOD_AUDIT_v7.csv"), bom = TRUE)

notes <- c(
  "3D placenta depth validation v7: explicit 2D-PC-matched PROMINENT-extrema depth validation with loading-direction concordance",
  "",
  "Primary evidence:",
  "- FIXED_BAND_depth_profiles_all_layers.csv",
  "- FIXED_BAND_Disease_minus_Control_layerwise_delta.csv",
  "- FIXED_BAND_depth_direction_consistency_SUMMARY.csv",
  "",
  "PCA evidence:",
  "- A balanced fixed reference PCA is fitted across all common depth layers for each disease pair.",
  "- Existing finalized 2D placenta disease-aligned PC1-PC3 loadings are read; the 2D PCA is not refitted.",
  "- For each disease, fixed 3D reference PC1-PC3 are matched one-to-one to the existing 2D PC1-PC3 axes by absolute loading correlation.",
  "- Every depth layer is projected into the same fixed 3D PC coordinates.",
  "- Independent layer-wise PCA is retained as unsupervised sensitivity analysis.",
  "- Layer-wise PC1-PC3 are matched one-to-one to fixed 3D reference PCs by loading correlation, then inherit the matched 2D PC identity.",
  "- Formal 3D PCA feature calling uses the same PROMINENT signed-loading-extremum rule as the revised 2D analysis.",
  "- SG7 is used for peak detection only; relative |loading| threshold = 20%; minimum peak distance = 20 cm^-1; no fixed Top-N.",
  "- Feature recurrence and loading-direction concordance are separate outputs; 3D loadings are never flipped merely to force agreement with 2D.",
  "- Mixed-sign prominent extrema within a target interval are reported as mixed and are not forced into a single direction.",
  "- loading^2 contribution is retained only as an auxiliary audit metric.",
  "",
  "Interpretation limits:",
  "- Z-stack layers are not independent biological replicates.",
  "- Positive disease-aligned loading is not equivalent to higher raw intensity.",
  "- Direct narrow-band intensity direction and PCA loading direction must be reported separately.",
  "- Full-field depth profiles do not prove blood-vessel-specific localization unless a registered vascular ROI is supplied.",
  "- Raman bands are overlapping chemical features and are not single-molecule proof."
)
writeLines(notes, file.path(output_root, "README_v7_interpretation_notes.txt"))

message("\nDone.")
message("Output folder: ", output_root)
message("Primary integrated table: ",
        file.path(out_integrated, "INTEGRATED_2D_to_3D_feature_reproducibility_summary.csv"))
