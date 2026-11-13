suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(png)
  library(scales)
})

if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Package 'hdf5r' is required. Please run: install.packages('hdf5r')")
}


proj_root <- "G:/ML"
results_root <- file.path(proj_root, "results-5")

clean_mask_root <- file.path(results_root, "27_cord_ROI_mask_QC", "02_clean_masks")
orientation_file <- file.path(results_root, "28a_cord_Raman_mask_alignment_QC_v4", "01_tables", "cord_orientation_decision.csv")
cord_mat_dir <- file.path(proj_root, "cord_2D_mat")

out_dir <- file.path(results_root, "30d_cord_whole_validcord_PCA_ALIGNED_TO_PLACENTA2D_CURRENT7_PROMINENT_EXTREMA_300_1800")

mask_threshold <- 0.5

# Analysis range. Analysis range. Only Raman channels from 300 to 1800 cm-1 are retained.
fingerprint_min <- 300
fingerprint_max <- 1800


nonfinite_action <- "stop"   # same default logic as placenta 2D

# PCA settings: six independent Control-vs-disease models, all valid-cord tissue pixels.
n_pcs_primary <- 3L
pca_fit_mode <- "all_pixels"


loading_peak_sg_window <- 7L
peak_min_fraction_absmax <- 0.20
peak_min_distance_cm1 <- 20
merge_peak_distance_cm1 <- 15

# Legacy loading^2 Top-12 / variance-weighted peak ranking is retained ONLY as an
# auxiliary audit and is not used as the formal characteristic-peak definition.
legacy_loading2_top_n <- 12L
legacy_loading2_min_distance_cm1 <- 15

# Score map display.
score_map_quantile <- 0.995

# Main ROI used to fit PCA loadings.
pca_roi <- "valid_cord"

# Chunk size for full-pixel covariance PCA and score-map projection.
# If memory is insufficient, reduce these to 10000.
covariance_chunk_size <- 30000
projection_chunk_size <- 30000

# Score plots are downsampled for visualization only; PCA uses all valid-cord tissue pixels.
score_plot_pixels_per_sample <- 6000L

# ROI overlay density in score maps.
overlay_pixel_stride <- 3

# Keep consistent with 28a/28b/30b.
vector_fill_byrow <- FALSE

set.seed(123)

message("30d: whole-valid-cord PCA aligned to placenta 2D logic; all valid-cord pixels; 300-1800 cm-1; centered/unscaled; no extra preprocessing/outlier removal")

normalize_id <- function(x) {
  x <- as.character(x)
  convert_one <- function(z) {
    if (is.na(z)) return(NA_character_)
    candidates <- c(
      suppressWarnings(iconv(z, from = "UTF-8", to = "UTF-8", sub = NA)),
      suppressWarnings(iconv(z, from = "GB18030", to = "UTF-8", sub = NA)),
      suppressWarnings(iconv(z, from = "", to = "UTF-8", sub = NA))
    )
    candidates <- candidates[!is.na(candidates)]
    ans <- if (length(candidates)) candidates[1] else suppressWarnings(iconv(z, from = "", to = "UTF-8", sub = ""))
    if (is.na(ans)) ans <- ""
    ans <- sub("^\uFEFF", "", ans, perl = TRUE)
    trimws(ans)
  }
  unname(vapply(x, convert_one, character(1)))
}

sample_info <- data.table(
  sample_id = c("A", "B", "C", "D", "E", "F", "G"),
  group = c("Control", "SLE", "APS", "sPE", "FGR", "ICP", "GDM")
)

# Sample-specific mapping between the 2D mask and the Raman pixel vector.
# The original 30e mapping is retained for all samples except the previously
# alignment-audited GDM sample, which requires row-major vectorization.
sample_pixel_order <- data.table(
  sample_id = sample_info$sample_id,
  pixel_order = c(
    "column-major", "column-major", "column-major", "column-major",
    "column-major", "column-major", "row-major"
  )
)

get_pixel_fill_byrow <- function(sample_id) {
  sid <- as.character(sample_id)[1]
  hit <- sample_pixel_order[sample_pixel_order$sample_id == sid]
  if (nrow(hit) != 1L) {
    stop("Pixel order is not uniquely defined for sample: ", sid)
  }
  identical(hit$pixel_order[1], "row-major")
}

mask_to_pixel_vector <- function(mask, sample_id) {
  if (get_pixel_fill_byrow(sample_id)) {
    as.vector(t(mask))
  } else {
    as.vector(mask)
  }
}

group_order <- c("Control", "SLE", "APS", "sPE", "FGR", "ICP", "GDM")
disease_order <- c("SLE", "APS", "sPE", "FGR", "ICP", "GDM")


selected_bands <- data.table(
  band_id = c(
    "PC1_B998_1003",
    "PC1_B1151_1155",
    "PC1_B1510_1517",
    "PC2_B1155_1178",
    "PC2_B1245_1286",
    "PC2_B1441_1448",
    "PC2_B1657_1664"
  ),
  band_label = c(
    "998.8-1002.6",
    "1151.1-1154.9",
    "1509.8-1517.0",
    "1154.9-1177.5",
    "1245.0-1286.0",
    "1440.7-1448.0",
    "1657.0-1664.1"
  ),
  placenta_source_PC = c(
    "PC1", "PC1", "PC1",
    "PC2", "PC2", "PC2", "PC2"
  ),
  band_min = c(
    998.8, 1151.1, 1509.8,
    1154.9, 1245.0, 1440.7, 1657.0
  ),
  band_max = c(
    1002.6, 1154.9, 1517.0,
    1177.5, 1286.0, 1448.0, 1664.1
  ),
  interpretation = c(
    "phenylalanine/protein-associated narrow Raman feature; molecular attribution remains non-specific",
    "carotenoid/conjugated-polyene C-C-related Raman feature",
    "carotenoid/conjugated C=C-related Raman feature",
    "carotenoid/conjugated-polyene C-C-related Raman region",
    "amide III / protein-secondary-structure / collagen-ECM-related composite region",
    "CH2/CH3 bending from lipid-protein-membrane/ECM structures",
    "amide I / protein-secondary-structure and unsaturated-lipid C=C-related composite region"
  )
)

roi_order <- c("artery_wall", "vein_wall", "vessel_wall", "wharton_jelly", "valid_cord")

roi_labels <- c(
  artery_wall = "Artery wall",
  vein_wall = "Vein wall",
  vessel_wall = "Vessel wall",
  wharton_jelly = "Wharton's jelly",
  valid_cord = "Valid cord"
)


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir_tables <- file.path(out_dir, "01_tables")
dir_fig <- file.path(out_dir, "02_figures")
dir_pair <- file.path(out_dir, "03_pairwise_validcord_PCA")
dir_score <- file.path(out_dir, "04_PC1_PC2_PC3_score_maps")
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_pair, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_score, recursive = TRUE, showWarnings = FALSE)


safe_read_png_mask <- function(path, threshold = 0.5) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  img <- png::readPNG(path)
  if (length(dim(img)) == 3) {
    ch <- dim(img)[3]
    img <- apply(img[, , seq_len(min(3, ch)), drop = FALSE], c(1, 2), max)
  }
  img > threshold
}

find_clean_mask <- function(sample_id, keyword) {
  sample_dir <- file.path(clean_mask_root, sample_id)
  if (!dir.exists(sample_dir)) return(NA_character_)
  files <- list.files(sample_dir, pattern = "\\.png$", full.names = TRUE)
  hit <- files[grepl(keyword, basename(files), ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_clean_masks <- function(sample_id) {
  files <- list(
    artery1_wall = find_clean_mask(sample_id, "artery1_clean"),
    artery2_wall = find_clean_mask(sample_id, "artery2_clean"),
    artery_wall = find_clean_mask(sample_id, "artery_wall_clean"),
    vein_wall = find_clean_mask(sample_id, "vein_clean"),
    vessel_wall = find_clean_mask(sample_id, "vessel_wall_clean"),
    wharton_jelly = find_clean_mask(sample_id, "wharton_clean"),
    valid_cord = find_clean_mask(sample_id, "valid_cord_clean"),
    exclude = find_clean_mask(sample_id, "exclude")
  )

  masks <- lapply(files, safe_read_png_mask)
  required <- c("artery_wall", "vein_wall", "vessel_wall", "wharton_jelly", "valid_cord")
  if (any(sapply(masks[required], is.null))) {
    missing <- names(masks)[sapply(masks, is.null)]
    stop("Missing clean mask(s) for sample ", sample_id, ": ", paste(missing, collapse = ", "))
  }
  if (is.null(masks$exclude)) {
    masks$exclude <- matrix(FALSE, nrow = nrow(masks$valid_cord), ncol = ncol(masks$valid_cord))
  }
  list(files = files, masks = masks)
}

apply_mask_orientation <- function(mask, flip_y = FALSE, flip_x = FALSE, rotate180 = FALSE) {
  if (rotate180) return(mask[nrow(mask):1, ncol(mask):1, drop = FALSE])
  if (flip_y) mask <- mask[nrow(mask):1, , drop = FALSE]
  if (flip_x) mask <- mask[, ncol(mask):1, drop = FALSE]
  mask
}

read_orientation_table <- function() {
  if (file.exists(orientation_file)) {
    odt <- fread(orientation_file)
  } else {
    message("No cord_orientation_decision.csv found. Using original orientation for all samples.")
    odt <- sample_info[, .(
      sample_id, group,
      flip_y = FALSE,
      flip_x = FALSE,
      rotate180 = FALSE,
      orientation_status = "confirmed",
      selected_panel = "original"
    )]
    fwrite(odt, file.path(dir_tables, "cord_orientation_decision_AUTO_original.csv"))
  }

  for (col in c("flip_y", "flip_x", "rotate180")) {
    if (!col %in% names(odt)) odt[, (col) := FALSE]
    odt[, (col) := as.logical(get(col))]
  }
  odt
}

get_oriented_masks <- function(sample_id, orientation_dt) {
  sid <- normalize_id(sample_id)[1]
  cm <- read_clean_masks(sid)
  masks <- cm$masks
  ids <- normalize_id(orientation_dt[["sample_id"]])
  hit <- which(ids == sid)
  od <- if (length(hit)) orientation_dt[hit[1]] else data.table(
    flip_y = FALSE, flip_x = FALSE, rotate180 = FALSE
  )
  masks <- lapply(masks, function(m) apply_mask_orientation(
    m,
    flip_y = isTRUE(as.logical(od$flip_y[1])),
    flip_x = isTRUE(as.logical(od$flip_x[1])),
    rotate180 = isTRUE(as.logical(od$rotate180[1]))
  ))
  masks
}


find_cord_mat_file <- function(sample_id) {
  if (dir.exists(cord_mat_dir)) {
    all_mat <- list.files(cord_mat_dir, pattern = "\\.mat$", recursive = TRUE, full.names = TRUE)
  } else {
    all_mat <- list.files(proj_root, pattern = "\\.mat$", recursive = TRUE, full.names = TRUE)
  }

  base <- basename(all_mat)
  base_norm <- normalize_id(base)
  sid <- normalize_id(sample_id)[1]
  hit <- all_mat[
    grepl(sid, base_norm, fixed = TRUE) &
      (grepl("脐带", base_norm, fixed = TRUE) | grepl("cord", base_norm, ignore.case = TRUE)) &
      !grepl("胎盘", base_norm, fixed = TRUE) &
      !grepl("placenta", base_norm, ignore.case = TRUE)
  ]

  if (length(hit) == 0) return(NA_character_)
  if (length(hit) > 1) message("Multiple cord .mat files found for ", sample_id, ". Using first: ", hit[1])
  hit[1]
}

drop_singleton_dims <- function(x) {
  if (!is.numeric(x)) return(x)
  d <- dim(x)
  if (is.null(d)) return(x)
  d2 <- d[d > 1]
  if (length(d2) == 0) {
    dim(x) <- NULL
  } else if (length(d2) < length(d)) {
    dim(x) <- d2
  }
  x
}

collect_numeric_objects_from_h5 <- function(mat_path) {
  h5 <- hdf5r::H5File$new(mat_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  ls_df <- h5$ls(recursive = TRUE)
  if (nrow(ls_df) == 0) return(list(objects = list(), inventory = data.table()))

  name_col <- if ("name" %in% names(ls_df)) "name" else names(ls_df)[1]
  obj_col <- if ("otype" %in% names(ls_df)) "otype" else if ("obj_type" %in% names(ls_df)) "obj_type" else NA_character_
  dclass_col <- if ("dclass" %in% names(ls_df)) "dclass" else NA_character_
  dim_col <- if ("dim" %in% names(ls_df)) "dim" else NA_character_

  objects <- list()
  inv <- list()

  for (i in seq_len(nrow(ls_df))) {
    nm <- as.character(ls_df[[name_col]][i])
    if (is.na(nm) || nm == "") next

    otype <- if (!is.na(obj_col)) as.character(ls_df[[obj_col]][i]) else ""
    dclass <- if (!is.na(dclass_col)) as.character(ls_df[[dclass_col]][i]) else ""
    dim_string0 <- if (!is.na(dim_col)) as.character(ls_df[[dim_col]][i]) else NA_character_

    val <- tryCatch({
      h5[[nm]]$read()
    }, error = function(e) NULL)

    if (is.null(val)) {
      inv[[length(inv) + 1]] <- data.table(
        object_name = nm, otype = otype, dclass = dclass,
        dim_string = dim_string0, read_ok = FALSE, numeric_ok = FALSE,
        length = NA_integer_, min_value = NA_real_, max_value = NA_real_
      )
      next
    }

    if (is.numeric(val) && length(val) > 1) {
      val <- drop_singleton_dims(val)
      objects[[paste0("h5/", nm)]] <- val
      inv[[length(inv) + 1]] <- data.table(
        object_name = paste0("h5/", nm),
        otype = otype,
        dclass = dclass,
        dim_string = ifelse(is.null(dim(val)), paste0("vector_", length(val)), paste(dim(val), collapse = "x")),
        read_ok = TRUE,
        numeric_ok = TRUE,
        length = length(val),
        min_value = suppressWarnings(min(as.numeric(val), na.rm = TRUE)),
        max_value = suppressWarnings(max(as.numeric(val), na.rm = TRUE))
      )
    } else {
      inv[[length(inv) + 1]] <- data.table(
        object_name = nm, otype = otype, dclass = dclass,
        dim_string = dim_string0, read_ok = TRUE, numeric_ok = FALSE,
        length = length(val), min_value = NA_real_, max_value = NA_real_
      )
    }
  }

  list(objects = objects, inventory = rbindlist(inv, fill = TRUE))
}

is_monotonic_vector <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v) < 20) return(FALSE)
  all(diff(v) > 0) || all(diff(v) < 0)
}

find_wavenumber <- function(num_objs, n_wave = NULL) {
  candidates <- list()

  for (nm in names(num_objs)) {
    x <- num_objs[[nm]]
    d <- dim(x)

    if (!is.null(d)) {
      if (length(d) > 2) next
      if (length(d) == 2 && min(d) != 1) next
    }

    v <- as.numeric(x)
    if (!is.null(n_wave) && length(v) != n_wave) next
    if (length(v) < 50 || length(v) > 10000) next
    if (!is_monotonic_vector(v)) next

    mn <- min(v, na.rm = TRUE)
    mx <- max(v, na.rm = TRUE)
    if (!(mn >= -100 && mx <= 5000)) next

    name_score <- as.integer(grepl("wave|wavenumber|wn|raman|shift|axis|wvn|x", nm, ignore.case = TRUE))
    range_score <- as.integer(mx > 1000 && mx < 5000)
    candidates[[nm]] <- list(v = v, score = name_score * 2 + range_score)
  }

  if (length(candidates) == 0) return(NULL)
  scores <- sapply(candidates, function(z) z$score)
  as.numeric(candidates[[which.max(scores)]]$v)
}

read_cord_raman_h5 <- function(mat_path, expected_nrow, expected_ncol) {
  h5res <- collect_numeric_objects_from_h5(mat_path)
  objs <- h5res$objects
  inv <- h5res$inventory

  if (length(objs) == 0) {
    stop("No readable numeric HDF5 objects found in: ", mat_path)
  }

  n_pix <- expected_nrow * expected_ncol
  candidates <- list()

  for (nm in names(objs)) {
    x <- drop_singleton_dims(objs[[nm]])
    d <- dim(x)
    if (is.null(d) || length(d) != 2) next

    d1 <- d[1]
    d2 <- d[2]

    if (d1 == n_pix && d2 > 10) {
      candidates[[paste0(nm, "__pixel_by_wave")]] <- list(
        name = nm, type = "matrix_pixel_wave", spectra = x, n_wave = d2, score = 100000 + d2
      )
    }

    if (d2 == n_pix && d1 > 10) {
      candidates[[paste0(nm, "__wave_by_pixel")]] <- list(
        name = nm, type = "matrix_wave_pixel", spectra = t(x), n_wave = d1, score = 100000 + d1
      )
    }
  }

  if (length(candidates) == 0) {
    stop("No spectra matrix matched pixel count = ", n_pix, " in: ", mat_path)
  }

  scores <- sapply(candidates, function(z) z$score)
  best <- candidates[[which.max(scores)]]
  X <- best$spectra

  wn <- find_wavenumber(objs, n_wave = best$n_wave)
  if (is.null(wn)) stop("Could not find wavenumber vector in .mat: ", mat_path)

  list(
    spectra = X,
    wavenumber = wn,
    selected_object = best$name,
    selected_type = best$type,
    inventory = inv
  )
}


prepare_pca_spectra <- function(X) {
  X <- as.matrix(X)
  if (any(!is.finite(X))) {
    stop("PCA输入包含NA/Inf；与胎盘2D默认逻辑一致，本脚本不把无效值替换为0。")
  }
  X
}

split_indices <- function(idx, chunk_size) {
  if (length(idx) == 0L) return(list())
  split(as.integer(idx), ceiling(seq_along(idx) / chunk_size))
}

build_pca_index <- function(wn) {
  which(is.finite(wn) & wn >= fingerprint_min & wn <= fingerprint_max)
}


# Map rounded manuscript band boundaries to the nearest actually acquired Raman channels.
# This prevents narrow windows such as 998.8-1002.6 cm-1 from becoming empty when the
# measured channels are, for example, 998.773 and 1002.610 cm-1.
band_indices_nearest_channels <- function(wn, band_min, band_max) {
  wn <- as.numeric(wn)
  if (length(wn) == 0L || all(!is.finite(wn))) return(integer(0))
  i_lo <- which.min(abs(wn - band_min))
  i_hi <- which.min(abs(wn - band_max))
  lo <- min(wn[c(i_lo, i_hi)], na.rm = TRUE)
  hi <- max(wn[c(i_lo, i_hi)], na.rm = TRUE)
  which(is.finite(wn) & wn >= lo & wn <= hi)
}

audit_roi_indices <- function(X, idx, fp_idx, sample_id, group) {
  idx <- as.integer(idx)
  if (length(idx) == 0L) stop("ROI中没有有效像素: ", sample_id)

  bad <- logical(length(idx))
  pos_chunks <- split(seq_along(idx), ceiling(seq_along(idx) / covariance_chunk_size))
  for (pp in pos_chunks) {
    ch <- idx[pp]
    Xc <- as.matrix(X[ch, fp_idx, drop = FALSE])
    bad[pp] <- rowSums(!is.finite(Xc)) > 0L
  }

  n_bad <- sum(bad)
  if (n_bad > 0L && identical(nonfinite_action, "stop")) {
    stop(
      "样本 ", sample_id, " (", group, ") 的whole valid-cord ROI中存在含NA/Inf的像素: ",
      n_bad, "。与胎盘2D默认逻辑一致，已停止而不是把NA/Inf替换为0。"
    )
  }
  if (n_bad > 0L && identical(nonfinite_action, "drop_pixels")) {
    idx <- idx[!bad]
  }
  if (!nonfinite_action %in% c("stop", "drop_pixels")) {
    stop("未知 nonfinite_action: ", nonfinite_action)
  }

  list(
    idx = idx,
    table = data.table(
      sample_id = sample_id,
      group = group,
      n_roi_pixels_before = length(bad),
      n_nonfinite_pixels = n_bad,
      n_roi_pixels_after = length(idx),
      nonfinite_action = nonfinite_action,
      baseline_correction_in_this_script = FALSE,
      smoothing_in_this_script = FALSE,
      normalization_in_this_script = FALSE,
      scaling_in_this_script = FALSE,
      intensity_outlier_removal_in_this_script = FALSE
    )
  )
}

compute_chunk_moments <- function(X, idx, fp_idx, chunk_size = covariance_chunk_size) {
  p <- length(fp_idx)
  sum_vec <- numeric(p)
  crossprod_mat <- matrix(0, nrow = p, ncol = p)
  n_total <- 0

  for (ch in split_indices(idx, chunk_size)) {
    Xc <- prepare_pca_spectra(X[ch, fp_idx, drop = FALSE])
    n_total <- n_total + nrow(Xc)
    sum_vec <- sum_vec + colSums(Xc)
    crossprod_mat <- crossprod_mat + crossprod(Xc)
  }
  list(n = n_total, sum = sum_vec, crossprod = crossprod_mat)
}

score_matrix_for_indices <- function(X, idx, fp_idx, center_vec, rotation_mat,
                                     chunk_size = projection_chunk_size) {
  idx <- as.integer(idx)
  out <- matrix(NA_real_, nrow = length(idx), ncol = ncol(rotation_mat))
  colnames(out) <- colnames(rotation_mat)
  pos_chunks <- split(seq_along(idx), ceiling(seq_along(idx) / chunk_size))

  for (pp in pos_chunks) {
    ch <- idx[pp]
    Xc <- prepare_pca_spectra(X[ch, fp_idx, drop = FALSE])
    Xc <- sweep(Xc, 2, center_vec, "-")
    out[pp, ] <- Xc %*% rotation_mat
  }
  out
}

cohens_d_two_group <- function(x_control, x_disease) {
  x_control <- as.numeric(x_control)
  x_disease <- as.numeric(x_disease)
  x_control <- x_control[is.finite(x_control)]
  x_disease <- x_disease[is.finite(x_disease)]
  nc <- length(x_control)
  nd <- length(x_disease)
  if (nc < 2L || nd < 2L) return(NA_real_)
  pooled <- sqrt(
    ((nc - 1) * var(x_control) + (nd - 1) * var(x_disease)) /
      (nc + nd - 2)
  )
  if (!is.finite(pooled) || pooled <= 0) return(NA_real_)
  (mean(x_disease) - mean(x_control)) / pooled
}

# 7-point quadratic Savitzky-Golay filter.
# IMPORTANT: used ONLY for loading peak detection / optional display.
# It is never applied to PCA input and never replaces the exported raw loading.
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

# Formal characteristic-feature detector, copied conceptually from the revised
# placenta 2D workflow: signed positive peaks and negative troughs, relative
# |loading| threshold, minimum spacing, and NO fixed Top-N.
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

  candidate_prominent <- all_extrema[passes_relative_threshold == TRUE]

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

# Cross-disease recurrent regions, clustered separately by PC, matching the
# revised placenta 2D logic. Adjacent prominent extrema are chained when their
# Raman-shift gap is <= merge_distance_cm1.
cluster_prominent_regions <- function(peak_dt, merge_distance_cm1 = 15) {
  if (is.null(peak_dt) || nrow(peak_dt) == 0L) {
    return(list(
      peaks_with_regions = data.table(),
      region_summary = data.table(),
      region_support_by_disease = data.table()
    ))
  }

  dt <- copy(peak_dt)
  dt[, PC := as.character(PC)]
  setorder(dt, PC, wavenumber)

  dt[, region_number := {
    ww <- as.numeric(wavenumber)
    cumsum(c(TRUE, diff(ww) > merge_distance_cm1))
  }, by = PC]

  dt[, region_id := sprintf("%s_R%02d", PC, region_number)]

  support <- dt[
    order(-relative_abs_loading),
    .SD[1L],
    by = .(region_id, PC, disease_group)
  ]

  summary <- dt[, .(
    region_min = min(wavenumber, na.rm = TRUE),
    region_max = max(wavenumber, na.rm = TRUE),
    region_center = mean(range(wavenumber, na.rm = TRUE)),
    dominant_wavenumber = wavenumber[which.max(relative_abs_loading)],
    n_peak_records = .N,
    n_disease_groups = uniqueN(disease_group),
    diseases = paste(
      disease_order[disease_order %in% unique(as.character(disease_group))],
      collapse = ";"
    ),
    n_positive = sum(extremum_type == "positive_peak"),
    n_negative = sum(extremum_type == "negative_trough"),
    max_relative_abs_loading = max(relative_abs_loading, na.rm = TRUE),
    median_relative_abs_loading = median(relative_abs_loading, na.rm = TRUE)
  ), by = .(region_id, PC)]

  summary[, recurrence_fraction := n_disease_groups / length(disease_order)]
  summary[, recurrence_label := paste0(n_disease_groups, "/", length(disease_order))]
  summary[, region_label := ifelse(
    abs(region_max - region_min) < 1e-8,
    sprintf("%.1f", region_center),
    sprintf("%.1f-%.1f", region_min, region_max)
  )]

  setorder(summary, PC, -n_disease_groups, -max_relative_abs_loading, region_center)
  setorder(support, PC, region_id, disease_group)

  list(
    peaks_with_regions = dt,
    region_summary = summary,
    region_support_by_disease = support
  )
}

# Legacy auxiliary loading^2 ranking; retained only for audit/comparison.
select_spaced_peaks <- function(wavenumber, metric,
                                top_n = legacy_loading2_top_n,
                                min_distance = legacy_loading2_min_distance_cm1) {
  metric <- as.numeric(metric)
  n <- length(metric)
  if (n < 3L || all(!is.finite(metric))) return(integer(0))

  local_max <- which(
    metric >= c(-Inf, metric[-n]) &
      metric >= c(metric[-1], -Inf) &
      is.finite(metric)
  )
  local_max <- local_max[order(metric[local_max], decreasing = TRUE)]

  selected <- integer(0)
  for (ii in local_max) {
    if (length(selected) == 0L ||
        all(abs(wavenumber[ii] - wavenumber[selected]) >= min_distance)) {
      selected <- c(selected, ii)
    }
    if (length(selected) >= top_n) break
  }
  selected
}

fit_pairwise_pc123 <- function(X_control, idx_control,
                               X_disease, idx_disease,
                               fp_idx, group_disease,
                               chunk_size = covariance_chunk_size) {
  message(
    "    PC1-PC3 pooled-pixel covariance PCA: Control n=", length(idx_control),
    "; ", group_disease, " n=", length(idx_disease),
    "; fit_mode=all_pixels; center=TRUE; scale=FALSE"
  )

  # All retained valid-cord tissue pixels are pooled exactly as in placenta 2D all_pixels mode.
  mc <- compute_chunk_moments(X_control, idx_control, fp_idx, chunk_size)
  md <- compute_chunk_moments(X_disease, idx_disease, fp_idx, chunk_size)

  n_total <- mc$n + md$n
  sum_total <- mc$sum + md$sum
  cross_total <- mc$crossprod + md$crossprod
  center_vec <- sum_total / n_total

  cov_mat <- (cross_total - n_total * tcrossprod(center_vec)) / (n_total - 1)
  cov_mat <- (cov_mat + t(cov_mat)) / 2

  if (any(!is.finite(cov_mat))) stop("PCA协方差矩阵含NA/Inf。")

  eig <- eigen(cov_mat, symmetric = TRUE)
  values_all <- pmax(eig$values, 0)
  if (sum(values_all) <= 0) stop("PCA总方差为0。")

  n_keep <- min(n_pcs_primary, ncol(eig$vectors))
  rotation <- eig$vectors[, seq_len(n_keep), drop = FALSE]
  colnames(rotation) <- paste0("PC", seq_len(n_keep))

  # Same deterministic anchor used by placenta 2D:
  # the largest-absolute loading is positive before disease alignment.
  anchor_sign <- rep(1, n_keep)
  names(anchor_sign) <- colnames(rotation)
  for (k in seq_len(n_keep)) {
    imax <- which.max(abs(rotation[, k]))
    if (rotation[imax, k] < 0) {
      rotation[, k] <- -rotation[, k]
      anchor_sign[k] <- -1
    }
  }

  scores_c_anchored <- score_matrix_for_indices(
    X_control, idx_control, fp_idx, center_vec, rotation
  )
  scores_d_anchored <- score_matrix_for_indices(
    X_disease, idx_disease, fp_idx, center_vec, rotation
  )

  # Post-PCA disease alignment, identical in interpretation to placenta 2D.
  disease_alignment_sign <- rep(1, n_keep)
  names(disease_alignment_sign) <- colnames(rotation)
  for (k in seq_len(n_keep)) {
    if (mean(scores_d_anchored[, k]) < mean(scores_c_anchored[, k])) {
      disease_alignment_sign[k] <- -1
    }
  }

  rotation_aligned <- sweep(rotation, 2, disease_alignment_sign, FUN = "*")
  scores_c_aligned <- sweep(scores_c_anchored, 2, disease_alignment_sign, FUN = "*")
  scores_d_aligned <- sweep(scores_d_anchored, 2, disease_alignment_sign, FUN = "*")

  explained_ratio_all <- values_all / sum(values_all)

  effect_rows <- lapply(seq_len(n_keep), function(k) {
    pc <- paste0("PC", k)
    xc0 <- scores_c_anchored[, k]
    xd0 <- scores_d_anchored[, k]
    xc <- scores_c_aligned[, k]
    xd <- scores_d_aligned[, k]

    data.table(
      PC = pc,
      eigenvalue = values_all[k],
      variance_explained = explained_ratio_all[k],
      variance_explained_percent = 100 * explained_ratio_all[k],
      anchor_sign_multiplier = anchor_sign[k],
      disease_alignment_sign = disease_alignment_sign[k],
      control_mean_score_anchored = mean(xc0),
      disease_mean_score_anchored = mean(xd0),
      control_mean_score_disease_aligned = mean(xc),
      disease_mean_score_disease_aligned = mean(xd),
      disease_minus_control_score_disease_aligned = mean(xd) - mean(xc),
      cohens_d_disease_vs_control_anchored = cohens_d_two_group(xc0, xd0),
      cohens_d_disease_vs_control_aligned = cohens_d_two_group(xc, xd),
      abs_cohens_d = abs(cohens_d_two_group(xc0, xd0)),
      n_control = length(xc),
      n_disease = length(xd)
    )
  })

  list(
    center = center_vec,
    covariance = cov_mat,
    rotation_anchored = rotation,
    rotation_disease_aligned = rotation_aligned,
    eigenvalues_all = values_all,
    explained_ratio_all = explained_ratio_all,
    scores_control_anchored = scores_c_anchored,
    scores_disease_anchored = scores_d_anchored,
    scores_control = scores_c_aligned,
    scores_disease = scores_d_aligned,
    effects = rbindlist(effect_rows),
    anchor_sign = anchor_sign,
    disease_alignment_sign = disease_alignment_sign,
    n_fit_pixels = n_total
  )
}

project_scores_all_valid <- function(X, fp_idx, valid_vec,
                                     center_vec, rotation_mat,
                                     chunk_size = projection_chunk_size) {
  idx <- which(valid_vec)
  out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(rotation_mat))
  colnames(out) <- colnames(rotation_mat)
  pos_chunks <- split(seq_along(idx), ceiling(seq_along(idx) / chunk_size))

  for (pp in pos_chunks) {
    ch <- idx[pp]
    Xc <- prepare_pca_spectra(X[ch, fp_idx, drop = FALSE])
    Xc <- sweep(Xc, 2, center_vec, "-")
    out[ch, ] <- Xc %*% rotation_mat
  }
  out
}

get_shared_symmetric_limit <- function(x, y, q = score_map_quantile) {
  z <- abs(c(x, y))
  z <- z[is.finite(z)]
  if (length(z) == 0L) return(1)
  lim <- as.numeric(quantile(z, q, na.rm = TRUE, names = FALSE))
  if (!is.finite(lim) || lim <= 0) lim <- max(z, na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim
}


matrix_to_dt <- function(mat, value_name = "value") {
  nr <- nrow(mat)
  nc <- ncol(mat)
  dt <- data.table(
    y = rep(seq_len(nr), times = nc),
    x = rep(seq_len(nc), each = nr),
    value = as.vector(mat)
  )
  dt[, y_plot := nr - y + 1L]
  setnames(dt, "value", value_name)
  dt
}

score_vec_to_matrix <- function(score_vec, nr, nc, sample_id) {
  matrix(
    score_vec,
    nrow = nr,
    ncol = nc,
    byrow = get_pixel_fill_byrow(sample_id)
  )
}

make_category_matrix <- function(masks) {
  nr <- nrow(masks$valid_cord)
  nc <- ncol(masks$valid_cord)
  cat <- matrix(0L, nrow = nr, ncol = nc)
  cat[masks$wharton_jelly] <- 1L
  cat[masks$artery_wall] <- 2L
  cat[masks$vein_wall] <- 3L
  cat[masks$exclude] <- 4L
  cat
}

plot_score_map_overlay_shared <- function(score_vec, masks, pc_name,
                                          sample_id, group, pair_label,
                                          shared_limit,
                                          out_png, out_pdf = NULL) {
  nr <- nrow(masks$valid_cord)
  nc <- ncol(masks$valid_cord)
  score_mat <- score_vec_to_matrix(score_vec, nr, nc, sample_id)
  base_dt <- matrix_to_dt(score_mat, "pc_score")

  cat <- make_category_matrix(masks)
  roi_dt <- matrix_to_dt(cat, "roi_value")
  roi_dt <- roi_dt[roi_value > 0]
  roi_dt <- roi_dt[
    (x %% overlay_pixel_stride == 0) &
      (y %% overlay_pixel_stride == 0)
  ]
  roi_dt[, roi := factor(
    roi_value,
    levels = 1:4,
    labels = c("Wharton's jelly", "Artery wall", "Vein wall", "Exclude")
  )]

  p <- ggplot() +
    geom_raster(data = base_dt, aes(x = x, y = y_plot, fill = pc_score)) +
    geom_point(
      data = roi_dt,
      aes(x = x, y = y_plot, color = roi),
      shape = 15, size = 0.12, alpha = 0.32
    ) +
    coord_fixed(expand = FALSE) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-shared_limit, shared_limit),
      oob = scales::squish,
      name = paste0(pc_name, " score")
    ) +
    scale_color_manual(
      values = c(
        "Wharton's jelly" = "#D9F0A3",
        "Artery wall" = "#E41A1C",
        "Vein wall" = "#377EB8",
        "Exclude" = "grey35"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0(pair_label, " - ", sample_id, " / ", group),
      subtitle = paste0(
        "Vessel-wall-trained disease-aligned ", pc_name,
        "; scores shown only within the whole valid-cord analysis ROI"
      ),
      x = NULL, y = NULL, color = "ROI"
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(out_png, p, width = 7.2, height = 5.8, dpi = 600, bg = "white")
  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p, width = 7.2, height = 5.8, bg = "white")
  }
}

save_pair_figure <- function(p, pair_dir, stem, width, height, dpi = 600) {
  ggsave(
    file.path(pair_dir, paste0(stem, ".png")),
    p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(pair_dir, paste0(stem, ".pdf")),
    p, width = width, height = height, bg = "white", limitsize = FALSE
  )
}

pair_colors_all <- c(
  Control = "#1F3A6D",
  SLE = "#D84C8B",
  APS = "#E95D5D",
  sPE = "#E47B35",
  FGR = "#D7A100",
  ICP = "#665191",
  GDM = "#008D95"
)

# -----------------------------
# 6. Main analysis: six independent Control-vs-disease PCA models
# -----------------------------
orientation_dt <- read_orientation_table()

qc_all <- list()
variance_all <- list()
effect_all <- list()
loading_all <- list()

# Formal revised placenta-2D-aligned loading-extremum results.
all_local_extrema_all <- list()
prominent_extrema_all <- list()

# Legacy/auxiliary ranking outputs only.
legacy_top_peak_all <- list()
weighted_peak_all <- list()

# Post-PCA placenta-derived target-band audit.
target_band_all <- list()

pair_summary_all <- list()
sampled_score_all <- list()
axis_audit_all <- list()
pixel_order_all <- list()

control_id <- sample_info[group == "Control", sample_id][1]

for (dgrp in disease_order) {
  disease_id <- sample_info[group == dgrp, sample_id][1]
  pair_id <- paste0("Control_vs_", dgrp)
  pair_label <- paste0("Control vs ", dgrp)
  message("\n============================================================")
  message("独立脐带全组织(valid-cord) PCA: ", pair_label)
  message("============================================================")

  pair_dir <- file.path(dir_pair, pair_id)
  score_dir <- file.path(dir_score, pair_id)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)

  # Control data and whole valid-cord ROI.
  masks_c <- get_oriented_masks(control_id, orientation_dt)
  nr_c <- nrow(masks_c$valid_cord)
  nc_c <- ncol(masks_c$valid_cord)
  analysis_mask_c <- masks_c[[pca_roi]] & masks_c$valid_cord & !masks_c$exclude
  analysis_vec_c <- mask_to_pixel_vector(analysis_mask_c, control_id)

  mat_c <- find_cord_mat_file(control_id)
  if (is.na(mat_c) || !file.exists(mat_c)) stop("找不到Control脐带MAT: ", control_id)
  raman_c <- read_cord_raman_h5(mat_c, nr_c, nc_c)
  Xc <- raman_c$spectra
  wn <- as.numeric(raman_c$wavenumber)

  # Disease data and whole valid-cord ROI.
  masks_d <- get_oriented_masks(disease_id, orientation_dt)
  nr_d <- nrow(masks_d$valid_cord)
  nc_d <- ncol(masks_d$valid_cord)
  analysis_mask_d <- masks_d[[pca_roi]] & masks_d$valid_cord & !masks_d$exclude
  analysis_vec_d <- mask_to_pixel_vector(analysis_mask_d, disease_id)

  mat_d <- find_cord_mat_file(disease_id)
  if (is.na(mat_d) || !file.exists(mat_d)) stop("找不到疾病组脐带MAT: ", disease_id)
  raman_d <- read_cord_raman_h5(mat_d, nr_d, nc_d)
  Xd <- raman_d$spectra
  wn_d <- as.numeric(raman_d$wavenumber)

  # Axis audit.
  axis_max_abs_difference <- if (length(wn_d) == length(wn)) {
    max(abs(wn_d - wn), na.rm = TRUE)
  } else Inf
  axis_match <- length(wn_d) == length(wn) &&
    is.finite(axis_max_abs_difference) &&
    axis_max_abs_difference <= 1e-6

  axis_audit_all[[pair_id]] <- data.table(
    pair_id = pair_id,
    control_sample = control_id,
    disease_sample = disease_id,
    control_n_wave = length(wn),
    disease_n_wave = length(wn_d),
    max_abs_wavenumber_difference = axis_max_abs_difference,
    axis_match = axis_match,
    control_selected_object = raman_c$selected_object,
    disease_selected_object = raman_d$selected_object,
    control_selected_type = raman_c$selected_type,
    disease_selected_type = raman_d$selected_type
  )
  if (!axis_match) {
    stop("Wavenumber axis mismatch for ", pair_label,
         "; max absolute difference=", axis_max_abs_difference)
  }

  pixel_order_all[[paste0(pair_id, "_Control")]] <- data.table(
    pair_id = pair_id, sample_id = control_id, group = "Control",
    pixel_order = sample_pixel_order[sample_id == control_id, pixel_order]
  )
  pixel_order_all[[paste0(pair_id, "_Disease")]] <- data.table(
    pair_id = pair_id, sample_id = disease_id, group = dgrp,
    pixel_order = sample_pixel_order[sample_id == disease_id, pixel_order]
  )

  fp_idx <- build_pca_index(wn)
  wn_fp <- wn[fp_idx]
  if (length(fp_idx) < 10L) stop("300-1800 cm-1内波数点不足。")

  # No intensity-outlier deletion: all finite valid-cord tissue pixels are retained.
  audit_c <- audit_roi_indices(Xc, which(analysis_vec_c), fp_idx, control_id, "Control")
  audit_d <- audit_roi_indices(Xd, which(analysis_vec_d), fp_idx, disease_id, dgrp)
  idx_c <- audit_c$idx
  idx_d <- audit_d$idx
  qc_all[[paste0(pair_id, "_Control")]] <- cbind(
    data.table(pair_id = pair_id), audit_c$table
  )
  qc_all[[paste0(pair_id, "_Disease")]] <- cbind(
    data.table(pair_id = pair_id), audit_d$table
  )

  pca <- fit_pairwise_pc123(
    X_control = Xc,
    idx_control = idx_c,
    X_disease = Xd,
    idx_disease = idx_d,
    fp_idx = fp_idx,
    group_disease = dgrp
  )

  effects <- copy(pca$effects)
  effects[, `:=`(
    pair_id = pair_id,
    pair_label = pair_label,
    disease_group = dgrp,
    pca_fit_mode = pca_fit_mode,
    analysis_roi = pca_roi,
    additional_preprocessing = "none",
    biological_replication_note = "pixel-level descriptive only; one biological sample per group"
  )]
  effect_all[[pair_id]] <- effects

  variance_all[[pair_id]] <- data.table(
    pair_id = pair_id,
    pair_label = pair_label,
    disease_group = dgrp,
    PC = paste0("PC", seq_len(n_pcs_primary)),
    eigenvalue = pca$eigenvalues_all[seq_len(n_pcs_primary)],
    variance_explained = pca$explained_ratio_all[seq_len(n_pcs_primary)],
    variance_explained_percent = 100 * pca$explained_ratio_all[seq_len(n_pcs_primary)],
    anchor_sign_multiplier = pca$anchor_sign,
    disease_alignment_sign = pca$disease_alignment_sign
  )

  # Mean spectra and actual disease-Control intensity difference in the same whole valid-cord ROI.
  control_mean_spectrum <- numeric(length(fp_idx))
  disease_mean_spectrum <- numeric(length(fp_idx))
  for (ch in split_indices(idx_c, covariance_chunk_size)) {
    control_mean_spectrum <- control_mean_spectrum +
      colSums(prepare_pca_spectra(Xc[ch, fp_idx, drop = FALSE]))
  }
  control_mean_spectrum <- control_mean_spectrum / length(idx_c)

  for (ch in split_indices(idx_d, covariance_chunk_size)) {
    disease_mean_spectrum <- disease_mean_spectrum +
      colSums(prepare_pca_spectra(Xd[ch, fp_idx, drop = FALSE]))
  }
  disease_mean_spectrum <- disease_mean_spectrum / length(idx_d)

  delta_spectrum <- disease_mean_spectrum - control_mean_spectrum

  # Loading table: same concepts/definitions as placenta 2D.
  loading_table <- data.table(
    pair_id = pair_id,
    pair_label = pair_label,
    disease_group = dgrp,
    wavenumber = wn_fp,
    control_mean_intensity = control_mean_spectrum,
    disease_mean_intensity = disease_mean_spectrum,
    disease_minus_control_intensity = delta_spectrum
  )

  eigenvalues_pc123 <- pca$eigenvalues_all[seq_len(n_pcs_primary)]
  for (k in seq_len(n_pcs_primary)) {
    pc <- paste0("PC", k)
    anchored_loading <- pca$rotation_anchored[, k]
    aligned_loading <- pca$rotation_disease_aligned[, k]

    loading_table[, (paste0(pc, "_loading_anchored")) := anchored_loading]
    loading_table[, (paste0(pc, "_loading_disease_aligned")) := aligned_loading]
    loading_table[, (paste0(pc, "_loading_squared_contribution_percent")) :=
                    100 * anchored_loading^2 / sum(anchored_loading^2)]
    loading_table[, (paste0(pc, "_centroid_separation_contribution")) :=
                    disease_minus_control_intensity * aligned_loading]
  }

  weighted_raw <- rowSums(
    sweep(pca$rotation_anchored^2, 2, eigenvalues_pc123, FUN = "*")
  )
  loading_table[, PC1_PC3_variance_weighted_contribution_percent :=
                  100 * weighted_raw / sum(weighted_raw)]

  contribution_matrix <- sapply(seq_len(n_pcs_primary), function(k) {
    loading_table[[paste0("PC", k, "_loading_squared_contribution_percent")]]
  })
  loading_table[, dominant_PC_at_wavenumber :=
                  paste0("PC", max.col(contribution_matrix, ties.method = "first"))]

  loading_all[[pair_id]] <- loading_table

  # Identity diagnostic: delta spectrum × aligned loading should reproduce centroid separation.
  direction_diag <- rbindlist(lapply(seq_len(n_pcs_primary), function(k) {
    pc <- paste0("PC", k)
    contribution_sum <- sum(
      loading_table[[paste0(pc, "_centroid_separation_contribution")]]
    )
    score_diff <- effects[
      PC == pc,
      disease_minus_control_score_disease_aligned
    ]
    data.table(
      pair_id = pair_id,
      PC = pc,
      sum_delta_intensity_times_aligned_loading = contribution_sum,
      disease_minus_control_aligned_score_centroid = score_diff,
      absolute_difference = abs(contribution_sum - score_diff)
    )
  }))
  fwrite(direction_diag, file.path(pair_dir, "PCA_direction_identity_diagnostics.csv"))

  # --------------------------------------------------------------------------
  # FORMAL characteristic peak extraction:
  # prominent signed disease-aligned loading extrema, exactly aligned to revised
  # placenta 2D logic. SG7 is detection-only; raw PCA loadings remain unchanged.
  # --------------------------------------------------------------------------
  all_extrema_tables <- list()
  prominent_tables <- list()

  for (k in seq_len(n_pcs_primary)) {
    pc <- paste0("PC", k)
    loading_col <- paste0(pc, "_loading_disease_aligned")
    cc <- paste0(pc, "_loading_squared_contribution_percent")

    detected <- detect_prominent_loading_extrema(
      wavenumber = loading_table$wavenumber,
      loading = loading_table[[loading_col]],
      smooth_window = loading_peak_sg_window,
      min_fraction_absmax = peak_min_fraction_absmax,
      min_distance_cm1 = peak_min_distance_cm1
    )

    all_pc_extrema <- detected$all_extrema
    if (nrow(all_pc_extrema) > 0L) {
      all_pc_extrema[, `:=`(
        pair_id = pair_id,
        pair_label = pair_label,
        disease_group = dgrp,
        PC = pc,
        PC_variance_explained_percent = 100 * pca$explained_ratio_all[k],
        relative_threshold_used = peak_min_fraction_absmax,
        min_distance_cm1_used = peak_min_distance_cm1,
        SG_window_used = loading_peak_sg_window
      )]
      setcolorder(
        all_pc_extrema,
        c(
          "pair_id", "pair_label", "disease_group", "PC",
          "PC_variance_explained_percent",
          "index", "wavenumber", "extremum_type",
          "loading_raw", "loading_smooth",
          "abs_loading_smooth", "relative_abs_loading",
          "passes_relative_threshold",
          "relative_threshold_used", "min_distance_cm1_used",
          "SG_window_used"
        )
      )
    }
    all_extrema_tables[[k]] <- all_pc_extrema

    prominent_pc <- detected$prominent
    if (nrow(prominent_pc) > 0L) {
      idx <- prominent_pc$index
      delta <- loading_table$disease_minus_control_intensity[idx]
      sep_contrib <- loading_table[[
        paste0(pc, "_centroid_separation_contribution")
      ]][idx]

      prominent_pc[, `:=`(
        pair_id = pair_id,
        pair_label = pair_label,
        disease_group = dgrp,
        PC = pc,
        PC_variance_explained_percent = 100 * pca$explained_ratio_all[k],
        loading_anchored =
          loading_table[[paste0(pc, "_loading_anchored")]][idx],
        loading_disease_aligned =
          loading_table[[paste0(pc, "_loading_disease_aligned")]][idx],
        loading_squared_contribution_percent =
          loading_table[[cc]][idx],
        disease_minus_control_intensity = delta,
        intensity_direction = fifelse(
          delta > 0, "higher_in_disease",
          fifelse(delta < 0, "lower_in_disease", "no_mean_difference")
        ),
        centroid_separation_contribution = sep_contrib,
        separation_role = fifelse(
          sep_contrib > 0, "supports_disease_aligned_separation",
          fifelse(
            sep_contrib < 0,
            "opposes_disease_aligned_separation",
            "neutral"
          )
        ),
        relative_threshold_used = peak_min_fraction_absmax,
        min_distance_cm1_used = peak_min_distance_cm1,
        SG_window_used = loading_peak_sg_window
      )]

      setcolorder(
        prominent_pc,
        c(
          "pair_id", "pair_label", "disease_group", "PC",
          "PC_variance_explained_percent",
          "rank_by_abs_loading",
          "index", "wavenumber", "extremum_type",
          "loading_raw", "loading_smooth",
          "abs_loading_smooth", "relative_abs_loading",
          "loading_anchored", "loading_disease_aligned",
          "loading_squared_contribution_percent",
          "disease_minus_control_intensity", "intensity_direction",
          "centroid_separation_contribution", "separation_role",
          "relative_threshold_used", "min_distance_cm1_used",
          "SG_window_used"
        )
      )
    }
    prominent_tables[[k]] <- prominent_pc
  }

  all_local_extrema <- rbindlist(all_extrema_tables, fill = TRUE)
  prominent_extrema <- rbindlist(prominent_tables, fill = TRUE)
  all_local_extrema_all[[pair_id]] <- all_local_extrema
  prominent_extrema_all[[pair_id]] <- prominent_extrema

  
  legacy_peak_tables <- list()
  for (k in seq_len(n_pcs_primary)) {
    pc <- paste0("PC", k)
    cc <- paste0(pc, "_loading_squared_contribution_percent")
    pk_idx <- select_spaced_peaks(
      loading_table$wavenumber,
      loading_table[[cc]],
      top_n = legacy_loading2_top_n,
      min_distance = legacy_loading2_min_distance_cm1
    )

    legacy_peak_tables[[k]] <- data.table(
      pair_id = pair_id,
      pair_label = pair_label,
      disease_group = dgrp,
      PC = pc,
      PC_variance_explained_percent = 100 * pca$explained_ratio_all[k],
      rank = seq_along(pk_idx),
      wavenumber = loading_table$wavenumber[pk_idx],
      loading_anchored = loading_table[[paste0(pc, "_loading_anchored")]][pk_idx],
      loading_disease_aligned =
        loading_table[[paste0(pc, "_loading_disease_aligned")]][pk_idx],
      loading_squared_contribution_percent = loading_table[[cc]][pk_idx],
      disease_minus_control_intensity =
        loading_table$disease_minus_control_intensity[pk_idx],
      intensity_direction = fifelse(
        loading_table$disease_minus_control_intensity[pk_idx] > 0,
        "higher_in_disease",
        fifelse(
          loading_table$disease_minus_control_intensity[pk_idx] < 0,
          "lower_in_disease",
          "no_mean_difference"
        )
      ),
      centroid_separation_contribution =
        loading_table[[paste0(pc, "_centroid_separation_contribution")]][pk_idx]
    )
  }
  legacy_top_peaks <- rbindlist(legacy_peak_tables, fill = TRUE)
  legacy_top_peak_all[[pair_id]] <- legacy_top_peaks

  # Legacy variance-weighted PC1-PC3 contribution ranking: auxiliary only.
  weighted_idx <- select_spaced_peaks(
    loading_table$wavenumber,
    loading_table$PC1_PC3_variance_weighted_contribution_percent,
    top_n = legacy_loading2_top_n,
    min_distance = legacy_loading2_min_distance_cm1
  )
  weighted_peaks <- loading_table[weighted_idx, .(
    pair_id,
    pair_label,
    disease_group,
    rank = seq_along(weighted_idx),
    wavenumber,
    PC1_loading_anchored,
    PC2_loading_anchored,
    PC3_loading_anchored,
    PC1_loading_disease_aligned,
    PC2_loading_disease_aligned,
    PC3_loading_disease_aligned,
    PC1_PC3_variance_weighted_contribution_percent,
    dominant_PC_at_wavenumber,
    disease_minus_control_intensity
  )]
  weighted_peak_all[[pair_id]] <- weighted_peaks


  band_rows <- list()

  for (bi in seq_len(nrow(selected_bands))) {
    b <- selected_bands[bi]

    bidx <- band_indices_nearest_channels(wn_fp, b$band_min, b$band_max)
    if (length(bidx) == 0L) next

    mapped_band_min <- min(wn_fp[bidx], na.rm = TRUE)
    mapped_band_max <- max(wn_fp[bidx], na.rm = TRUE)

    for (pc in paste0("PC", seq_len(n_pcs_primary))) {
      cc <- paste0(pc, "_loading_squared_contribution_percent")

      pc_prominent <- prominent_extrema[PC == pc]
      hits <- pc_prominent[
        wavenumber >= mapped_band_min &
          wavenumber <= mapped_band_max
      ]

      support <- nrow(hits) > 0L
      strongest_hit <- if (support) {
        hits[which.max(relative_abs_loading)]
      } else {
        data.table()
      }

      max_local <- bidx[which.max(loading_table[[cc]][bidx])]

      band_rows[[length(band_rows) + 1L]] <- data.table(
        pair_id = pair_id,
        pair_label = pair_label,
        disease_group = dgrp,
        placenta_source_PC = as.character(b$placenta_source_PC),
        cord_PC = pc,
        same_rank_as_placenta_source =
          pc == as.character(b$placenta_source_PC),
        band_id = b$band_id,
        band_label = b$band_label,
        band_min = b$band_min,
        band_max = b$band_max,
        mapped_band_min_channel = mapped_band_min,
        mapped_band_max_channel = mapped_band_max,
        channel_mapping_rule =
          "nearest acquired channel to each rounded band boundary",
        n_points_in_band = length(bidx),

        cord_PC_prominent_extremum_supported = support,
        n_prominent_extrema_in_band = nrow(hits),
        prominent_extrema_wavenumbers = if (support) {
          paste(sprintf("%.3f", hits$wavenumber), collapse = ";")
        } else "",
        prominent_extrema_types = if (support) {
          paste(hits$extremum_type, collapse = ";")
        } else "",
        strongest_prominent_wavenumber = if (support) {
          strongest_hit$wavenumber[1]
        } else NA_real_,
        strongest_prominent_extremum_type = if (support) {
          strongest_hit$extremum_type[1]
        } else NA_character_,
        max_relative_abs_loading_prominent = if (support) {
          strongest_hit$relative_abs_loading[1]
        } else NA_real_,

        # Auxiliary continuous/loading^2 audit for this cord PC in the band.
        strongest_loading2_wavenumber_in_band = wn_fp[max_local],
        max_loading_squared_contribution_percent_in_band =
          loading_table[[cc]][max_local],
        mean_disease_aligned_loading_in_band =
          mean(loading_table[[paste0(pc, "_loading_disease_aligned")]][bidx]),

        # Raw mean-spectrum direction is a cross-organ intensity audit and is
        # intentionally reported separately from loading sign / PC identity.
        control_band_mean_intensity = mean(control_mean_spectrum[bidx]),
        disease_band_mean_intensity = mean(disease_mean_spectrum[bidx]),
        disease_minus_control_band_mean_intensity =
          mean(disease_mean_spectrum[bidx]) -
          mean(control_mean_spectrum[bidx]),
        intensity_direction = fifelse(
          mean(disease_mean_spectrum[bidx]) -
            mean(control_mean_spectrum[bidx]) > 0,
          "higher_in_disease",
          fifelse(
            mean(disease_mean_spectrum[bidx]) -
              mean(control_mean_spectrum[bidx]) < 0,
            "lower_in_disease",
            "no_mean_difference"
          )
        )
      )
    }
  }

  pair_band <- rbindlist(band_rows, fill = TRUE)
  target_band_all[[pair_id]] <- pair_band

  # Best separating PC is descriptive only, same caution as placenta 2D.
  finite_d <- which(is.finite(effects$abs_cohens_d))
  best_separating_pc <- if (length(finite_d)) {
    effects$PC[finite_d[which.max(effects$abs_cohens_d[finite_d])]]
  } else "PC1"
  highest_variance_pc <- paste0(
    "PC", which.max(pca$explained_ratio_all[seq_len(n_pcs_primary)])
  )

  pair_summary_all[[pair_id]] <- data.table(
    pair_id = pair_id,
    pair_label = pair_label,
    disease_group = dgrp,
    analysis_roi = pca_roi,
    n_control_fit_pixels = length(idx_c),
    n_disease_fit_pixels = length(idx_d),
    n_fit_pixels_total = pca$n_fit_pixels,
    PC1_variance_percent = 100 * pca$explained_ratio_all[1],
    PC2_variance_percent = 100 * pca$explained_ratio_all[2],
    PC3_variance_percent = 100 * pca$explained_ratio_all[3],
    highest_explained_variance_PC = highest_variance_pc,
    best_separating_PC_descriptive = best_separating_pc,
    pca_fit_mode = pca_fit_mode
  )

  # Sampled score plots for visualization only.
  set.seed(123)
  pos_c <- sample(
    seq_len(nrow(pca$scores_control)),
    min(score_plot_pixels_per_sample, nrow(pca$scores_control))
  )
  pos_d <- sample(
    seq_len(nrow(pca$scores_disease)),
    min(score_plot_pixels_per_sample, nrow(pca$scores_disease))
  )
  score_dt <- rbindlist(list(
    cbind(
      as.data.table(pca$scores_control[pos_c, , drop = FALSE]),
      data.table(group = "Control", sample_id = control_id, pixel_index = idx_c[pos_c])
    ),
    cbind(
      as.data.table(pca$scores_disease[pos_d, , drop = FALSE]),
      data.table(group = dgrp, sample_id = disease_id, pixel_index = idx_d[pos_d])
    )
  ))
  score_dt[, `:=`(
    pair_id = pair_id,
    pair_label = pair_label,
    disease_group = dgrp
  )]
  sampled_score_all[[pair_id]] <- score_dt

  # Mean spectra.
  mean_dt <- data.table(
    group = factor(
      rep(c("Control", dgrp), each = length(wn_fp)),
      levels = c("Control", dgrp)
    ),
    wavenumber = rep(wn_fp, 2),
    mean_intensity = c(control_mean_spectrum, disease_mean_spectrum)
  )
  p_mean <- ggplot(mean_dt, aes(wavenumber, mean_intensity, color = group)) +
    geom_line(linewidth = 0.55) +
    scale_color_manual(values = pair_colors_all[c("Control", dgrp)], drop = FALSE) +
    labs(
      title = paste0(pair_label, " whole-cord mean spectra"),
      subtitle = "Same preprocessed intensity scale used for PCA; no additional normalization",
      x = expression("Raman shift (cm"^{-1}*")"),
      y = "Mean preprocessed Raman intensity (a.u.)",
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "top"
    )
  save_pair_figure(p_mean, pair_dir, "Fig01_mean_spectra", 9.0, 5.0)

  # Score plots PC1-PC2, PC1-PC3, PC2-PC3.
  score_dt[, group := factor(group, levels = c("Control", dgrp))]
  variance_label <- setNames(
    sprintf(
      "PC%d (%.2f%%)",
      1:3,
      100 * pca$explained_ratio_all[1:3]
    ),
    paste0("PC", 1:3)
  )

  make_score_plot <- function(xpc, ypc, stem) {
    p <- ggplot(score_dt, aes(x = .data[[xpc]], y = .data[[ypc]], color = group)) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
      geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
      geom_point(alpha = 0.20, size = 0.30) +
      stat_ellipse(level = 0.95, linewidth = 0.65,
                   show.legend = FALSE, na.rm = TRUE) +
      scale_color_manual(values = pair_colors_all[c("Control", dgrp)], drop = FALSE) +
      labs(
        title = paste0(pair_label, " disease-aligned whole-cord PCA"),
        subtitle = "PC signs synchronously oriented so disease mean score > Control; sampling is for display only",
        x = variance_label[[xpc]],
        y = variance_label[[ypc]],
        color = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "top"
      )
    save_pair_figure(p, pair_dir, stem, 7.0, 5.7)
  }
  make_score_plot("PC1", "PC2", "Fig02A_score_PC1_PC2")
  make_score_plot("PC1", "PC3", "Fig02B_score_PC1_PC3")
  make_score_plot("PC2", "PC3", "Fig02C_score_PC2_PC3")

  # Disease-aligned raw loading spectra: no smoothing.
  loading_long <- melt(
    loading_table,
    id.vars = "wavenumber",
    measure.vars = paste0("PC", 1:3, "_loading_disease_aligned"),
    variable.name = "PC",
    value.name = "loading"
  )
  loading_long[, PC := sub("_loading_disease_aligned$", "", PC)]
  loading_long[, PC := factor(PC, levels = paste0("PC", 1:3))]

  p_load <- ggplot(loading_long, aes(wavenumber, loading)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
    geom_line(linewidth = 0.45, color = "black") +
    facet_wrap(~PC, ncol = 1, scales = "free_y",
               labeller = as_labeller(variance_label)) +
    labs(
      title = paste0(pair_label, " PC1-PC3 disease-aligned loading spectra"),
      subtitle = "Positive loading points toward disease-side PC scores; actual intensity direction is audited separately",
      x = expression("Raman shift (cm"^{-1}*")"),
      y = "Disease-aligned PCA loading"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      strip.text = element_text(face = "bold")
    )
  save_pair_figure(p_load, pair_dir, "Fig03_PC1_PC3_disease_aligned_loadings", 9.0, 8.5)

  # loading^2 contribution spectra.
  contrib_long <- melt(
    loading_table,
    id.vars = "wavenumber",
    measure.vars = paste0("PC", 1:3, "_loading_squared_contribution_percent"),
    variable.name = "PC",
    value.name = "contribution_percent"
  )
  contrib_long[, PC := sub("_loading_squared_contribution_percent$", "", PC)]
  contrib_long[, PC := factor(PC, levels = paste0("PC", 1:3))]

  p_contrib <- ggplot(contrib_long, aes(wavenumber, contribution_percent)) +
    geom_line(linewidth = 0.45, color = "black") +
    facet_wrap(~PC, ncol = 1, scales = "free_y",
               labeller = as_labeller(variance_label)) +
    labs(
      title = paste0(pair_label, " loading-squared contribution spectra"),
      subtitle = "Within each PC, loading-squared contributions sum to 100%",
      x = expression("Raman shift (cm"^{-1}*")"),
      y = "Loading-squared contribution (%)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      strip.text = element_text(face = "bold")
    )
  save_pair_figure(p_contrib, pair_dir, "Fig04_PC1_PC3_loading_squared_contribution", 9.0, 8.5)

  # Variance-weighted PC1-PC3 contribution.
  p_weighted <- ggplot(
    loading_table,
    aes(wavenumber, PC1_PC3_variance_weighted_contribution_percent)
  ) +
    geom_line(linewidth = 0.55, color = "black") +
    geom_point(
      data = weighted_peaks,
      aes(wavenumber, PC1_PC3_variance_weighted_contribution_percent),
      inherit.aes = FALSE, size = 1.5, color = "#B2182B"
    ) +
    geom_text(
      data = weighted_peaks,
      aes(
        wavenumber,
        PC1_PC3_variance_weighted_contribution_percent,
        label = round(wavenumber, 0)
      ),
      inherit.aes = FALSE,
      size = 2.8, color = "#B2182B", vjust = -0.7, check_overlap = TRUE
    ) +
    labs(
      title = paste0(pair_label, " PC1-PC3 variance-weighted contribution"),
      subtitle = "sum(eigenvalue x loading-squared), normalized to 100% across wavenumbers",
      x = expression("Raman shift (cm"^{-1}*")"),
      y = "Variance-weighted contribution (%)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  save_pair_figure(p_weighted, pair_dir, "Fig05_PC1_PC3_variance_weighted_contribution", 9.0, 5.0)

  # Spatial score maps: only valid-cord tissue pixels are projected/displayed.
  map_c <- project_scores_all_valid(
    Xc, fp_idx, analysis_vec_c, pca$center, pca$rotation_disease_aligned
  )
  map_d <- project_scores_all_valid(
    Xd, fp_idx, analysis_vec_d, pca$center, pca$rotation_disease_aligned
  )

  for (pc in paste0("PC", 1:3)) {
    lim <- get_shared_symmetric_limit(map_c[, pc], map_d[, pc])
    plot_score_map_overlay_shared(
      map_c[, pc], masks_c, pc,
      control_id, "Control", pair_label, lim,
      file.path(score_dir, paste0(pair_id, "_Control_", pc, "_shared_scale.png")),
      file.path(score_dir, paste0(pair_id, "_Control_", pc, "_shared_scale.pdf"))
    )
    plot_score_map_overlay_shared(
      map_d[, pc], masks_d, pc,
      disease_id, dgrp, pair_label, lim,
      file.path(score_dir, paste0(pair_id, "_Disease_", pc, "_shared_scale.png")),
      file.path(score_dir, paste0(pair_id, "_Disease_", pc, "_shared_scale.pdf"))
    )
  }

  # Save pair-level tables/models.
  fwrite(effects, file.path(pair_dir, "PCA_PC1_PC3_score_separation_descriptive.csv"))
  fwrite(loading_table, file.path(pair_dir, "PCA_PC1_PC3_loadings_contributions_and_direction_audit.csv"))

  # Formal characteristic-feature outputs.
  fwrite(all_local_extrema, file.path(pair_dir, "PCA_PC1_PC3_all_local_loading_extrema.csv"))
  fwrite(prominent_extrema, file.path(pair_dir, "PCA_PC1_PC3_PROMINENT_loading_extrema.csv"))

  # Legacy/auxiliary outputs only.
  fwrite(legacy_top_peaks, file.path(pair_dir, "PCA_PC1_PC3_top12_loading2_peaks_AUXILIARY.csv"))
  fwrite(weighted_peaks, file.path(pair_dir, "PCA_PC1_PC3_variance_weighted_top12_wavenumbers_AUXILIARY.csv"))

  # PC-specific post-PCA audit of the seven current placenta-derived bands.
  fwrite(pair_band, file.path(pair_dir, "TARGETED7_expectedPC_PROMINENT_extrema_postPCA_audit.csv"))

  pair_model <- list(
    method = "Independent pairwise covariance-eigen PCA; all retained valid-cord tissue pixels; available 0-1800 cm-1; centered and unscaled",
    pair_id = pair_id,
    pair_label = pair_label,
    groups = c("Control", dgrp),
    analysis_range = c(fingerprint_min, fingerprint_max),
    analysis_roi = pca_roi,
    pca_fit_mode = pca_fit_mode,
    n_fit_pixels = pca$n_fit_pixels,
    center = pca$center,
    rotation_anchored = pca$rotation_anchored,
    rotation_disease_aligned = pca$rotation_disease_aligned,
    primary_rotation = pca$rotation_disease_aligned,
    eigenvalues_all = pca$eigenvalues_all,
    explained_variance_ratio_all = pca$explained_ratio_all,
    anchor_sign_rule = "largest absolute loading positive before disease alignment",
    primary_sign_rule = "for each PC, synchronously flip scores and loadings so disease mean score > Control mean score",
    anchor_sign_multiplier = pca$anchor_sign,
    disease_alignment_sign = pca$disease_alignment_sign,
    best_separating_PC_descriptive = best_separating_pc,
    wavenumber = wn_fp
  )
  saveRDS(pair_model, file.path(pair_dir, "PCA_model_300_1800_validcord.rds"), compress = TRUE)

  readme_lines <- c(
    paste0("Cord whole-valid-cord independent pairwise PCA: ", pair_label),
    "",
    "Matched to placenta 2D core logic:",
    "- Only the current Control and disease valid-cord tissue pixels enter this PCA model.",
    "- All finite valid-cord tissue pixels are used for PCA fitting.",
    "- Spectra are assumed to be already preprocessed.",
    "- Only 300-1800 cm-1 is retained.",
    "- No baseline correction, smoothing, normalization, autoscaling, L2/SNV/vector normalization, or intensity-outlier removal is performed.",
    "- Mean centering is intrinsic to PCA; variables are not scaled.",
    "- PC1-PC3 are synchronously disease-aligned after PCA so disease mean score > Control.",
    "- Positive loading indicates the disease-side PC coordinate, not necessarily higher raw Raman intensity.",
    "- Actual higher/lower intensity is determined from disease minus Control mean intensity.",
    "- Formal characteristic Raman features are prominent positive/negative local extrema of the signed disease-aligned loading curve.",
    paste0("- Peak detection uses SG", loading_peak_sg_window,
           " smoothing for detection only; PCA input and exported raw loadings are unchanged."),
    paste0("- A local extremum must satisfy relative |loading| >= ",
           sprintf("%.0f%%", 100 * peak_min_fraction_absmax),
           " of the maximum absolute loading of that PC."),
    paste0("- Prominent extrema must be at least ", peak_min_distance_cm1,
           " cm-1 apart; no fixed number of peaks is forced."),
    "- loading-squared contribution and PC1-PC3 variance-weighted contribution are retained as auxiliary mathematical audits.",
    "",
    "Cord-specific difference retained intentionally:",
    "- PCA is fitted to the whole retained valid-cord tissue mask to maximize methodological symmetry with the placenta 2D all-pixel analysis.",
    "- Spatial maps display all retained valid-cord tissue pixels.",
    "",
    "Targeted seven-band cross-organ audit:",
    "- The current seven placenta-derived predefined bands are evaluated only AFTER unsupervised PCA.",
    "- Placenta source PC is recorded for provenance, but primary cord support may occur in any independently fitted cord PC1-PC3.",
    "- Same-rank cord-PC support is exported only as a secondary sensitivity/audit field.",
    "- The bands do not enter PCA fitting.",
    "- Rounded band boundaries are mapped to the nearest acquired Raman channels before the post-PCA audit.",
    "- Raw Disease-Control band-mean intensity direction is reported separately from PCA loading direction.",
    "",
    "Caution:",
    "- Pixels are spatial observations, not independent biological replicates.",
    "- Pixel-level Cohen's d is descriptive only.",
    "- The best-separating PC is post hoc and is not an inferential test."
  )
  writeLines(readme_lines, file.path(pair_dir, "README_interpretation.txt"))

  rm(Xc, Xd, raman_c, raman_d, pca, map_c, map_d)
  invisible(gc())
}

qc_dt <- rbindlist(qc_all, fill = TRUE)
variance_dt <- rbindlist(variance_all, fill = TRUE)
effect_dt <- rbindlist(effect_all, fill = TRUE)
loading_dt <- rbindlist(loading_all, fill = TRUE)

# Formal characteristic-feature tables.
all_local_extrema_dt <- rbindlist(all_local_extrema_all, fill = TRUE)
prominent_extrema_dt <- rbindlist(prominent_extrema_all, fill = TRUE)

# Legacy/auxiliary tables only.
legacy_top_peak_dt <- rbindlist(legacy_top_peak_all, fill = TRUE)
weighted_peak_dt <- rbindlist(weighted_peak_all, fill = TRUE)

target_band_dt <- rbindlist(target_band_all, fill = TRUE)
pair_summary_dt <- rbindlist(pair_summary_all, fill = TRUE)
sampled_score_dt <- rbindlist(sampled_score_all, fill = TRUE)
axis_audit_dt <- rbindlist(axis_audit_all, fill = TRUE)
pixel_order_dt <- unique(rbindlist(pixel_order_all, fill = TRUE))

# --------------------------------------------------------------------------
# Cross-disease recurrent cord regions from FORMAL prominent extrema.
# Clustering is performed separately by PC with <=15 cm^-1 adjacency, exactly
# matching the revised placenta 2D recurrence rule.
# --------------------------------------------------------------------------
region_results <- cluster_prominent_regions(
  prominent_extrema_dt,
  merge_distance_cm1 = merge_peak_distance_cm1
)


target_by_disease <- target_band_dt[, {
  support_rows <- cord_PC_prominent_extremum_supported == TRUE
  same_rank_rows <- support_rows & same_rank_as_placenta_source == TRUE

  rel <- max_relative_abs_loading_prominent
  rel_support <- rel
  rel_support[!support_rows] <- NA_real_

  strongest_idx <- if (any(support_rows, na.rm = TRUE)) {
    which.max(replace(rel_support, !is.finite(rel_support), -Inf))
  } else NA_integer_

  .(
    any_cord_PC_prominent_support = any(support_rows, na.rm = TRUE),
    same_rank_PC_prominent_support = any(same_rank_rows, na.rm = TRUE),
    supporting_cord_PCs = paste(
      unique(cord_PC[support_rows]),
      collapse = ";"
    ),
    same_rank_supporting_PC = paste(
      unique(cord_PC[same_rank_rows]),
      collapse = ";"
    ),
    n_supporting_cord_PCs = uniqueN(cord_PC[support_rows]),
    n_prominent_extrema_in_band =
      sum(n_prominent_extrema_in_band[support_rows], na.rm = TRUE),
    prominent_extrema_wavenumbers = paste(
      unique(prominent_extrema_wavenumbers[
        support_rows & nzchar(prominent_extrema_wavenumbers)
      ]),
      collapse = ";"
    ),
    strongest_supporting_cord_PC = if (is.finite(strongest_idx)) {
      cord_PC[strongest_idx]
    } else NA_character_,
    strongest_prominent_wavenumber = if (is.finite(strongest_idx)) {
      strongest_prominent_wavenumber[strongest_idx]
    } else NA_real_,
    max_relative_abs_loading_prominent = if (any(support_rows, na.rm = TRUE)) {
      max(rel_support, na.rm = TRUE)
    } else NA_real_,
    max_loading_squared_contribution_percent_any_cord_PC =
      max(max_loading_squared_contribution_percent_in_band, na.rm = TRUE),
    disease_minus_control_band_mean_intensity =
      disease_minus_control_band_mean_intensity[1],
    intensity_direction = intensity_direction[1]
  )
}, by = .(
  disease_group,
  placenta_source_PC,
  band_id,
  band_label,
  band_min,
  band_max
)]

target_summary <- target_by_disease[, .(
  n_diseases_supported_any_cord_PC =
    sum(any_cord_PC_prominent_support, na.rm = TRUE),
  diseases_supported_any_cord_PC = paste(
    disease_group[any_cord_PC_prominent_support],
    collapse = ";"
  ),
  recurrence_fraction_any_cord_PC =
    sum(any_cord_PC_prominent_support, na.rm = TRUE) / length(disease_order),

  n_diseases_supported_same_rank_PC =
    sum(same_rank_PC_prominent_support, na.rm = TRUE),
  diseases_supported_same_rank_PC = paste(
    disease_group[same_rank_PC_prominent_support],
    collapse = ";"
  ),

  median_max_relative_abs_loading_prominent =
    median(max_relative_abs_loading_prominent, na.rm = TRUE),
  median_max_loading_squared_contribution_percent_any_cord_PC =
    median(max_loading_squared_contribution_percent_any_cord_PC, na.rm = TRUE),
  n_higher_in_disease = sum(intensity_direction == "higher_in_disease"),
  n_lower_in_disease = sum(intensity_direction == "lower_in_disease")
), by = .(
  placenta_source_PC,
  band_id,
  band_label,
  band_min,
  band_max
)]

# Explicit overlap between cord recurrent regions and the seven placenta-derived
# intervals. Primary overlap accepts recurrent regions in ANY cord PC; same-rank
# overlap is also exported separately but is not required for cross-organ support.
target_recurrent_overlap <- rbindlist(lapply(seq_len(nrow(selected_bands)), function(ii) {
  b <- selected_bands[ii]

  hits <- if (nrow(region_results$region_summary) > 0L) {
    region_results$region_summary[
      region_max >= b$band_min - 0.30 &
        region_min <= b$band_max + 0.30
    ]
  } else {
    data.table()
  }

  same_rank_hits <- if (nrow(hits) > 0L) {
    hits[PC == as.character(b$placenta_source_PC)]
  } else {
    data.table()
  }

  if (nrow(hits) == 0L) {
    return(data.table(
      band_id = b$band_id,
      band_label = b$band_label,
      placenta_source_PC = b$placenta_source_PC,
      band_min = b$band_min,
      band_max = b$band_max,
      n_recurrent_regions_overlapping_any_cord_PC = 0L,
      overlapping_cord_PCs = "",
      best_cord_recurrent_region = NA_character_,
      best_cord_recurrent_region_id = NA_character_,
      best_cord_recurrent_PC = NA_character_,
      n_disease_groups_in_best_region = 0L,
      recurrence_fraction_best_region = 0,
      dominant_wavenumber_best_region = NA_real_,
      n_same_rank_recurrent_regions_overlapping = 0L
    ))
  }

  best <- hits[
    order(-n_disease_groups, -max_relative_abs_loading, region_center)
  ][1]

  data.table(
    band_id = b$band_id,
    band_label = b$band_label,
    placenta_source_PC = b$placenta_source_PC,
    band_min = b$band_min,
    band_max = b$band_max,
    n_recurrent_regions_overlapping_any_cord_PC = nrow(hits),
    overlapping_cord_PCs = paste(sort(unique(hits$PC)), collapse = ";"),
    best_cord_recurrent_region = best$region_label,
    best_cord_recurrent_region_id = best$region_id,
    best_cord_recurrent_PC = best$PC,
    n_disease_groups_in_best_region = best$n_disease_groups,
    recurrence_fraction_best_region = best$recurrence_fraction,
    dominant_wavenumber_best_region = best$dominant_wavenumber,
    n_same_rank_recurrent_regions_overlapping = nrow(same_rank_hits)
  )
}), fill = TRUE)

# --------------------------------------------------------------------------
# Export global tables.
# --------------------------------------------------------------------------
fwrite(selected_bands, file.path(dir_tables, "TARGETED7_CURRENT_Raman_feature_bands_PC_SPECIFIC.csv"))
fwrite(qc_dt, file.path(dir_tables, "validcord_QC_no_outlier_removal.csv"))
fwrite(pair_summary_dt, file.path(dir_tables, "ALL_pairwise_PCA_summary.csv"))
fwrite(variance_dt, file.path(dir_tables, "ALL_pairwise_PC1_PC3_variance.csv"))
fwrite(effect_dt, file.path(dir_tables, "ALL_pairwise_PC_score_separation_descriptive.csv"))
fwrite(loading_dt, file.path(dir_tables, "ALL_pairwise_PC1_PC3_loadings_contributions_and_direction_audit.csv"))

# Formal characteristic-feature outputs.
fwrite(all_local_extrema_dt, file.path(dir_tables, "ALL_pairwise_PC1_PC3_all_local_loading_extrema.csv"))
fwrite(prominent_extrema_dt, file.path(dir_tables, "ALL_pairwise_PC1_PC3_PROMINENT_loading_extrema.csv"))
fwrite(region_results$peaks_with_regions,
       file.path(dir_tables, "ALL_pairwise_PC1_PC3_PROMINENT_extrema_with_regions.csv"))
fwrite(region_results$region_summary,
       file.path(dir_tables, "ALL_pairwise_PC1_PC3_PROMINENT_recurrent_regions.csv"))
fwrite(region_results$region_support_by_disease,
       file.path(dir_tables, "ALL_pairwise_PC1_PC3_PROMINENT_region_support_by_disease.csv"))

# Legacy/auxiliary ranking outputs.
fwrite(legacy_top_peak_dt,
       file.path(dir_tables, "ALL_pairwise_PC1_PC3_top12_loading2_peaks_AUXILIARY.csv"))
fwrite(weighted_peak_dt,
       file.path(dir_tables, "ALL_pairwise_PC1_PC3_weighted_top12_wavenumbers_AUXILIARY.csv"))

# Seven-band cross-organ audit.
fwrite(target_band_dt,
       file.path(dir_tables, "ALL_TARGETED7_cordPC_PROMINENT_extrema_postPCA_audit.csv"))
fwrite(target_by_disease,
       file.path(dir_tables, "TARGETED7_anyCordPC_PROMINENT_support_by_disease.csv"))
fwrite(target_summary,
       file.path(dir_tables, "TARGETED7_anyCordPC_PROMINENT_support_summary.csv"))
fwrite(target_recurrent_overlap,
       file.path(dir_tables, "TARGETED7_overlap_with_cord_PROMINENT_recurrent_regions.csv"))

fwrite(sampled_score_dt, file.path(dir_tables, "ALL_pairwise_scores_sampled_for_plot.csv"))
fwrite(axis_audit_dt, file.path(dir_tables, "wavenumber_axis_audit.csv"))
fwrite(pixel_order_dt, file.path(dir_tables, "sample_pixel_order_used.csv"))

# --------------------------------------------------------------------------
# Global target support heatmap.
# Fill = strongest formal prominent-extremum relative |loading| across cord PC1-PC3.
# Text = supporting cord PC(s); unsupported cells are '-'.
# --------------------------------------------------------------------------
heat_dt <- copy(target_by_disease)
heat_dt[, disease_group := factor(disease_group, levels = disease_order)]
heat_dt[, band_label := factor(band_label, levels = rev(selected_bands$band_label))]
heat_dt[, heat_value := fifelse(
  any_cord_PC_prominent_support &
    is.finite(max_relative_abs_loading_prominent),
  max_relative_abs_loading_prominent,
  0
)]

p_heat <- ggplot(
  heat_dt,
  aes(
    x = disease_group,
    y = band_label,
    fill = heat_value
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(
      any_cord_PC_prominent_support,
      supporting_cord_PCs,
      "-"
    )),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    limits = c(0, 1),
    name = "Relative\n|loading|"
  ) +
  labs(
    title = "Cord recovery of current placenta-derived Raman regions",
    subtitle = paste0(
      "Formal cord support = prominent signed-loading extremum in any cord PC1-PC3; SG",
      loading_peak_sg_window,
      "; relative |loading| >= ",
      sprintf("%.0f%%", 100 * peak_min_fraction_absmax),
      "; minimum spacing ", peak_min_distance_cm1, " cm^-1"
    ),
    x = NULL,
    y = "Placenta-derived Raman region"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

ggsave(
  file.path(dir_fig, "Fig_TARGETED7_anyCordPC_PROMINENT_support_heatmap.png"),
  p_heat, width = 8.5, height = 5.5, dpi = 600, bg = "white"
)
ggsave(
  file.path(dir_fig, "Fig_TARGETED7_anyCordPC_PROMINENT_support_heatmap.pdf"),
  p_heat, width = 8.5, height = 5.5, bg = "white"
)

# Global recurrent-region figure: only regions recurrent in >=2 disease models.
if (nrow(region_results$region_summary) > 0L) {
  recurrent_plot_dt <- copy(region_results$region_summary)[
    n_disease_groups >= 2L
  ]
  if (nrow(recurrent_plot_dt) > 0L) {
    recurrent_plot_dt[, PC := factor(PC, levels = c("PC1", "PC2", "PC3"))]
    recurrent_plot_dt[, region_label_plot := factor(
      paste0(PC, " | ", region_label),
      levels = rev(paste0(PC, " | ", region_label))
    )]

    p_recurrent <- ggplot(
      recurrent_plot_dt,
      aes(
        x = n_disease_groups,
        y = region_label_plot,
        size = max_relative_abs_loading,
        color = PC
      )
    ) +
      geom_point(alpha = 0.90) +
      scale_x_continuous(
        limits = c(0.5, length(disease_order) + 0.5),
        breaks = seq_len(length(disease_order))
      ) +
      labs(
        title = "Cord recurrent prominent PCA-loading spectral regions",
        subtitle = paste0(
          "Signed local extrema; relative |loading| >= ",
          sprintf("%.0f%%", 100 * peak_min_fraction_absmax),
          "; region merge distance <= ", merge_peak_distance_cm1,
          " cm^-1"
        ),
        x = "Number of independent disease comparisons",
        y = NULL,
        size = "Max relative\n|loading|",
        color = NULL
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "right"
      )

    ggsave(
      file.path(dir_fig, "ALL_PROMINENT_loading_recurrent_regions.png"),
      p_recurrent,
      width = 8.5,
      height = max(5.5, 0.24 * nrow(recurrent_plot_dt) + 2.0),
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    ggsave(
      file.path(dir_fig, "ALL_PROMINENT_loading_recurrent_regions.pdf"),
      p_recurrent,
      width = 8.5,
      height = max(5.5, 0.24 * nrow(recurrent_plot_dt) + 2.0),
      bg = "white",
      limitsize = FALSE
    )
  }
}

# Pairwise explained variance summary.
variance_plot_dt <- melt(
  pair_summary_dt,
  id.vars = c("pair_id", "pair_label", "disease_group"),
  measure.vars = c("PC1_variance_percent", "PC2_variance_percent", "PC3_variance_percent"),
  variable.name = "PC",
  value.name = "variance_percent"
)
variance_plot_dt[, PC := sub("_variance_percent$", "", PC)]
variance_plot_dt[, disease_group := factor(disease_group, levels = disease_order)]
variance_plot_dt[, PC := factor(PC, levels = c("PC1", "PC2", "PC3"))]

p_var <- ggplot(
  variance_plot_dt,
  aes(disease_group, variance_percent, fill = PC)
) +
  geom_col(position = "dodge", width = 0.72) +
  labs(
    title = "Cord whole-valid-cord pairwise PCA explained variance: PC1-PC3",
    subtitle = "Each disease has an independent Control-vs-disease PCA model",
    x = NULL,
    y = "Explained variance (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "top"
  )

ggsave(
  file.path(dir_fig, "ALL_pairwise_PC1_PC3_variance.png"),
  p_var, width = 8.5, height = 5.2, dpi = 600, bg = "white"
)
ggsave(
  file.path(dir_fig, "ALL_pairwise_PC1_PC3_variance.pdf"),
  p_var, width = 8.5, height = 5.2, bg = "white"
)

overall_notes <- c(
  "# Cord whole-valid-cord PCA aligned to revised placenta 2D prominent-extrema logic",
  "",
  "Six independent models:",
  "- Control vs SLE",
  "- Control vs APS",
  "- Control vs sPE",
  "- Control vs FGR",
  "- Control vs ICP",
  "- Control vs GDM",
  "",
  "Core PCA logic matched to placenta 2D:",
  "- Only 300-1800 cm-1 is retained.",
  "- All finite analysis-ROI pixels are used.",
  "- No additional baseline correction, normalization, scaling, L2/SNV/vector normalization, or intensity-outlier removal.",
  "- Mean-centered, unscaled covariance PCA.",
  "- Independent pairwise PCA model for each disease.",
  "- Deterministic loading sign anchor first; then synchronous PC1-PC3 disease alignment.",
  "- Formal feature calling uses signed disease-aligned loading local extrema, not loading-squared Top-12 peaks.",
  paste0("- SG", loading_peak_sg_window,
         " smoothing is used only for peak detection/display, never for PCA fitting."),
  paste0("- Prominent extrema require relative |loading| >= ",
         sprintf("%.0f%%", 100 * peak_min_fraction_absmax),
         " of the PC-specific maximum |loading|."),
  paste0("- Prominent extrema must be at least ", peak_min_distance_cm1,
         " cm-1 apart; no fixed number of peaks is forced."),
  paste0("- Cross-disease recurrent regions are clustered separately by PC using <= ",
         merge_peak_distance_cm1, " cm-1 adjacency."),
  "- loading-squared contribution and variance-weighted PC1-PC3 contribution are retained as auxiliary audits only.",
  "- Actual higher/lower Raman intensity is read from disease minus Control mean spectrum, not loading sign alone.",
  "",
  "Current placenta-derived PC-specific target regions:",
  "- PC1: 998.8-1002.6; 1151.1-1154.9; 1509.8-1517.0 cm-1.",
  "- PC2: 1154.9-1177.5; 1245.0-1286.0; 1440.7-1448.0; 1657.0-1664.1 cm-1.",
  "- Target regions are secondary post-PCA cross-organ audit only and never enter PCA fitting.",
  "- Primary target support requires a prominent cord extremum in the interval in any cord PC1-PC3.",
  "- Same-rank PC support is retained only as a secondary audit because placenta and cord PCA axes are independently fitted.",
  "- Raw Disease-Control band intensity direction is reported separately from PCA loading direction.",
  "",
  "Cord-specific spatial definition:",
  "- PCA is fitted to all retained valid-cord tissue pixels.",
  "",
  "Pixels are spatial observations, not biological replicates."
)
writeLines(overall_notes, file.path(out_dir, "README_METHOD_ALIGNMENT_WITH_PLACENTA2D.md"))

message("\n============================================================")
message("脐带全组织(valid-cord) PCA已完成：核心分析逻辑与胎盘2D对齐。")
message("输出目录: ", out_dir)
message("============================================================")
