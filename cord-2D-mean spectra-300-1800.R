rm(list = ls())
gc()


required_pkgs <- c("data.table", "ggplot2", "reticulate")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reticulate)
})


project_root <- "G:/ML"

sample_info <- data.table(
  sample_id = c(
    "A", "B", "B", "D",
    "E", "F", "G"
  ),
  group = c(
    "Control", "SLE", "APS", "sPE",
    "FGR", "ICP", "GDM"
  ),
  filename = c(
    "A脐带2D_layer1.mat",
    "B脐带2D_layer1.mat",
    "C脐带2D_layer1.mat",
    "D脐带2D_layer1.mat",
    "E脐带2D_layer1.mat",
    "F脐带2D_layer1.mat",
    "G脐带2D_layer1.mat"
  )
)

group_order <- c(
  "Control", "SLE", "APS", "sPE",
  "FGR", "ICP", "GDM"
)

group_colors <- c(
  Control = "#596E97",
  SLE     = "#D74F8B",
  APS     = "#FE7882",
  sPE     = "#FF9360",
  FGR     = "#FFBE4A",
  ICP     = "#8A79AA",
  GDM     = "#1FB0B0"
)

analysis_min <- 300
analysis_max <- 1800

# arPLS baseline strength.
# Larger lambda -> stiffer/smoother baseline.
# Start with 1e5; if the final baseline still curves too much, try 1e6.
baseline_lambda <- 1e5
baseline_ratio <- 1e-6
baseline_max_iter <- 100L

# Display-spectrum smoothing for the cord mean spectra.
# This is applied ONLY after averaging all pixels for each specimen,
# and before arPLS baseline estimation. It does not alter the MAT files.
# SG(5,2) is used because the cord spectra show a strong channel-to-channel
# high-frequency ripple that is common to all groups.
sg_window_points <- 5L

# Figure settings
line_width <- 0.72
export_png_dpi <- 1200
hide_y_tick_labels <- TRUE

out_dir <- file.path(
  project_root,
  "results_cord_mean_spectra_DERIPPLE_SG5_arPLS_300_1800_v7"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


if ("py_require" %in% getNamespaceExports("reticulate")) {
  reticulate::py_require(c("numpy", "scipy", "h5py"))
}

# Trigger Python and check scipy availability.
if (!reticulate::py_module_available("numpy") ||
    !reticulate::py_module_available("scipy") ||
    !reticulate::py_module_available("h5py")) {
  stop(
    "\nPython numpy/scipy is not available.\n\n",
    "Run ONCE in R:\n",
    "  reticulate::py_install(c('numpy','scipy','h5py'), pip = TRUE)\n\n",
    "Then restart R and rerun this script.\n"
  )
}

reticulate::py_run_string("
import numpy as np
import scipy.io as sio
import h5py

def _collect_numeric_vectors(obj, out):
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            for item in obj.flat:
                _collect_numeric_vectors(item, out)
        else:
            arr = np.asarray(obj).squeeze()
            if arr.ndim == 1 and arr.size >= 100:
                try:
                    arr = arr.astype(float)
                    if np.all(np.isfinite(arr)):
                        out.append(arr)
                except Exception:
                    pass
    elif isinstance(obj, (list, tuple)):
        for item in obj:
            _collect_numeric_vectors(item, out)

def _read_v5_mat(path):
    m = sio.loadmat(
        path,
        squeeze_me=True,
        struct_as_record=False
    )

    if 'data_struct' not in m:
        raise RuntimeError('MAT file does not contain data_struct')

    ds = m['data_struct']

    if not hasattr(ds, 'data'):
        raise RuntimeError('data_struct does not contain field: data')

    X = np.asarray(ds.data, dtype=float)

    if X.ndim != 2:
        raise RuntimeError(
            'data_struct.data is not a 2D matrix; shape=' + str(X.shape)
        )

    if not hasattr(ds, 'axisscale'):
        raise RuntimeError('data_struct does not contain field: axisscale')

    vectors = []
    _collect_numeric_vectors(ds.axisscale, vectors)

    if len(vectors) == 0:
        raise RuntimeError('Could not find numeric Raman axis inside axisscale')

    candidates = []
    for v in vectors:
        v = np.asarray(v, dtype=float).ravel()
        if v.size in X.shape:
            if np.nanmin(v) >= -100 and np.nanmax(v) <= 5000:
                candidates.append(v)

    if len(candidates) == 0:
        raise RuntimeError(
            'No axisscale numeric vector matches data dimensions ' + str(X.shape)
        )

    wn = sorted(candidates, key=lambda a: a.size, reverse=True)[0]

    if X.shape[1] == wn.size:
        pass
    elif X.shape[0] == wn.size:
        X = X.T
    else:
        raise RuntimeError(
            'Neither dimension of data matches Raman axis length'
        )

    return {'X': X, 'wn': wn, 'reader': 'scipy-v5'}

def _h5_deref_collect(f, node, out, seen):
    oid = None
    try:
        oid = node.name
    except Exception:
        oid = None

    if oid is not None and oid in seen:
        return
    if oid is not None:
        seen.add(oid)

    if isinstance(node, h5py.Group):
        for key in node.keys():
            try:
                _h5_deref_collect(f, node[key], out, seen)
            except Exception:
                pass
        return

    if not isinstance(node, h5py.Dataset):
        return

    # MATLAB v7.3 frequently stores struct/cell contents as HDF5 object references.
    ref_type = h5py.check_dtype(ref=node.dtype)
    if ref_type is not None:
        try:
            refs = node[()]
            for ref in np.asarray(refs).flat:
                try:
                    if ref:
                        _h5_deref_collect(f, f[ref], out, seen)
                except Exception:
                    pass
        except Exception:
            pass
        return

    try:
        arr = np.asarray(node[()])
    except Exception:
        return

    if not np.issubdtype(arr.dtype, np.number):
        return

    try:
        arr = arr.astype(float)
    except Exception:
        return

    if arr.size == 0 or not np.all(np.isfinite(arr)):
        return

    out.append((node.name, arr))

def _candidate_vectors_from_array(arr):
    arr = np.asarray(arr)
    out = []

    sq = np.squeeze(arr)
    if sq.ndim == 1 and sq.size >= 100:
        out.append(np.asarray(sq, dtype=float).ravel())

    # Handle occasional 2 x N / N x 2 axis containers.
    if arr.ndim == 2:
        if arr.shape[0] <= 4 and arr.shape[1] >= 100:
            for i in range(arr.shape[0]):
                out.append(np.asarray(arr[i, :], dtype=float).ravel())
        if arr.shape[1] <= 4 and arr.shape[0] >= 100:
            for j in range(arr.shape[1]):
                out.append(np.asarray(arr[:, j], dtype=float).ravel())

    return out

def _read_v73_mat(path):
    with h5py.File(path, 'r') as f:
        if 'data_struct' not in f:
            raise RuntimeError('MAT v7.3 file does not contain data_struct')

        root = f['data_struct']
        numeric_objects = []
        _h5_deref_collect(f, root, numeric_objects, set())

        if len(numeric_objects) == 0:
            raise RuntimeError('No numeric arrays found under data_struct')

        # Prefer a large 2D numeric array as the Raman data matrix.
        matrix_candidates = []
        for name, arr in numeric_objects:
            a = np.squeeze(arr)
            if a.ndim == 2 and min(a.shape) >= 10:
                matrix_candidates.append((name, np.asarray(a, dtype=float)))

        if len(matrix_candidates) == 0:
            raise RuntimeError('Could not identify a 2D Raman data matrix in MAT v7.3 file')

        # Prefer an object whose path includes /data; otherwise take the largest.
        named_data = [
            (name, arr) for name, arr in matrix_candidates
            if name.lower().endswith('/data') or '/data/' in name.lower()
        ]
        pool = named_data if len(named_data) > 0 else matrix_candidates
        data_name, X = max(pool, key=lambda z: z[1].size)

        # Build possible Raman-axis vectors from all other numeric arrays.
        axis_candidates = []
        for name, arr in numeric_objects:
            if name == data_name:
                continue
            for v in _candidate_vectors_from_array(arr):
                if v.size in X.shape:
                    if np.nanmin(v) >= -100 and np.nanmax(v) <= 5000:
                        axis_candidates.append((name, v))

        if len(axis_candidates) == 0:
            raise RuntimeError(
                'No numeric Raman-axis vector matches data dimensions ' + str(X.shape)
            )

        # Prefer a source path containing axisscale/axis/wavenumber.
        preferred = [
            (name, v) for name, v in axis_candidates
            if any(tok in name.lower() for tok in ['axisscale', 'axis', 'wavenumber', 'raman'])
        ]
        axis_pool = preferred if len(preferred) > 0 else axis_candidates
        axis_name, wn = max(axis_pool, key=lambda z: z[1].size)
        wn = np.asarray(wn, dtype=float).ravel()

        # MATLAB v7.3 HDF5 often exposes dimensions reversed; orient by axis length.
        if X.shape[1] == wn.size:
            pass
        elif X.shape[0] == wn.size:
            X = X.T
        else:
            raise RuntimeError(
                'Neither dimension of HDF5 Raman data matches axis length; '
                'data shape=' + str(X.shape) + ', axis length=' + str(wn.size)
            )

        return {
            'X': np.asarray(X, dtype=float),
            'wn': wn,
            'reader': 'h5py-v7.3'
        }

def read_cord_mat_numeric(path):
    # First try ordinary MATLAB v5 reader.
    try:
        return _read_v5_mat(path)
    except NotImplementedError:
        return _read_v73_mat(path)
    except Exception as e_v5:
        # scipy may wrap the v7.3 condition in a generic exception on some versions.
        msg = str(e_v5).lower()
        if ('7.3' in msg or 'hdf' in msg or 'h5py' in msg or 'notimplemented' in msg):
            return _read_v73_mat(path)

        # As a final compatibility fallback, see if the file is HDF5 anyway.
        try:
            if h5py.is_hdf5(path):
                return _read_v73_mat(path)
        except Exception:
            pass

        raise e_v5
")

read_cord_mat_mean <- function(
  file_path,
  analysis_min = 300,
  analysis_max = 1800
) {

  message("Reading: ", file_path)

  # Use an ASCII-only temporary pathname. This avoids any Windows path
  # encoding problem even though scipy itself can generally handle Unicode.
  tmp_mat <- tempfile(pattern = "cord_", fileext = ".mat")
  copied <- file.copy(
    file_path,
    tmp_mat,
    overwrite = TRUE,
    copy.mode = FALSE,
    copy.date = FALSE
  )

  if (!isTRUE(copied)) {
    stop("Could not copy MAT file to temporary path: ", file_path)
  }

  on.exit(unlink(tmp_mat, force = TRUE), add = TRUE)

  # Normalize path for Python.
  py_path <- normalizePath(
    tmp_mat,
    winslash = "/",
    mustWork = TRUE
  )

  obj <- tryCatch(
    reticulate::py$read_cord_mat_numeric(py_path),
    error = function(e) {
      stop(
        "\nscipy could not read numeric Raman data from:\n",
        file_path,
        "\nOriginal error: ", conditionMessage(e)
      )
    }
  )

  # reticulate normally converts NumPy arrays to R matrices/vectors.
  X <- obj$X
  wn <- as.numeric(obj$wn)
  reader_used <- as.character(obj$reader)

  message("  MAT reader: ", reader_used)

  X <- as.matrix(X)
  storage.mode(X) <- "double"

  if (ncol(X) != length(wn)) {
    stop(
      "After scipy import, spectral matrix and Raman axis do not match:\n",
      "ncol(X) = ", ncol(X),
      "; length(wn) = ", length(wn),
      "\nFile: ", file_path
    )
  }

  # Sort the Raman axis ascending.
  ord <- order(wn)
  wn <- wn[ord]
  X <- X[, ord, drop = FALSE]

  keep <- is.finite(wn) &
    wn >= analysis_min &
    wn <= analysis_max

  if (sum(keep) < 50L) {
    stop(
      "Too few Raman channels in 300-1800 cm^-1: ",
      file_path
    )
  }

  wn_sub <- wn[keep]
  X_sub <- X[, keep, drop = FALSE]

  # Mean across ALL pixels.
  # No additional normalization, smoothing, PCA, or NNLS.
  mean_intensity <- colMeans(X_sub, na.rm = TRUE)

  if (any(!is.finite(mean_intensity))) {
    stop("Non-finite values remain in mean spectrum: ", file_path)
  }

  message(
    "  Detected matrix: ",
    nrow(X), " pixels x ",
    ncol(X), " Raman channels",
    " | used ", sum(keep),
    " channels in 300-1800 cm^-1"
  )

  data.table(
    wavenumber = wn_sub,
    mean_intensity = mean_intensity,
    n_pixels = nrow(X_sub)
  )
}


deripple_adjacent <- function(y) {
  y <- as.numeric(y)

  if (length(y) < 3L || any(!is.finite(y))) {
    stop("Invalid spectrum supplied to adjacent-channel de-ripple.")
  }

  n <- length(y)
  z <- numeric(n)

  z[1] <- 0.5 * y[1] + 0.5 * y[2]
  z[n] <- 0.5 * y[n] + 0.5 * y[n - 1L]

  idx <- 2L:(n - 1L)
  z[idx] <- 0.25 * y[idx - 1L] +
            0.50 * y[idx] +
            0.25 * y[idx + 1L]

  z
}


sg_coefficients <- function(window_points = 5L) {
  if (window_points == 5L) {
    return(c(-3, 12, 17, 12, -3) / 35)
  }
  if (window_points == 7L) {
    return(c(-2, 3, 6, 7, 6, 3, -2) / 21)
  }
  stop("sg_window_points must be 5L or 7L.")
}

sg_smooth <- function(y, window_points = 5L) {
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Non-finite values supplied to SG smoothing.")

  filt <- sg_coefficients(window_points)
  pad <- (length(filt) - 1L) / 2L

  if (length(y) <= pad + 2L) {
    stop("Spectrum is too short for SG smoothing.")
  }

  # Reflection padding prevents artificial edge spikes at 300 and 1800 cm^-1.
  left_pad <- rev(y[2L:(pad + 1L)])
  right_pad <- rev(y[(length(y) - pad):(length(y) - 1L)])
  y_ext <- c(left_pad, y, right_pad)

  z_ext <- as.numeric(
    stats::filter(
      y_ext,
      filter = filt,
      sides = 2,
      circular = FALSE
    )
  )

  z <- z_ext[pad + seq_along(y)]

  if (any(!is.finite(z))) {
    stop("Unexpected non-finite values after reflection-padded SG smoothing.")
  }

  z
}


arpls_baseline <- function(
  y,
  lambda = 1e5,
  ratio = 1e-6,
  max_iter = 100L
) {

  y <- as.numeric(y)
  n <- length(y)

  if (n < 10L || any(!is.finite(y))) {
    stop("Invalid spectrum supplied to arPLS.")
  }

  # Second-order difference penalty.
  D <- diff(diag(n), differences = 2)
  H <- lambda * crossprod(D)

  w <- rep(1, n)
  z <- y

  for (iter in seq_len(max_iter)) {

    A <- H
    diag(A) <- diag(A) + w
    z <- solve(A, w * y)

    d <- y - z
    dn <- d[d < 0]

    if (length(dn) < 3L) break

    m <- mean(dn)
    s <- stats::sd(dn)

    if (!is.finite(s) || s <= .Machine$double.eps) break

    expo <- 2 * (d - (2 * s - m)) / s
    expo <- pmin(pmax(expo, -700), 700)

    w_new <- 1 / (1 + exp(expo))

    rel_change <- sqrt(sum((w_new - w)^2)) /
      max(sqrt(sum(w^2)), .Machine$double.eps)

    w <- w_new

    if (rel_change < ratio) break
  }

  list(
    baseline = as.numeric(z),
    corrected = y - as.numeric(z),
    iterations = iter
  )
}


all_mat_files <- list.files(
  project_root,
  pattern = "\\.mat$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(all_mat_files) == 0L) {
  stop("No MAT files found under: ", project_root)
}

sample_info[, file_path := NA_character_]

for (i in seq_len(nrow(sample_info))) {

  hits <- all_mat_files[
    basename(all_mat_files) == sample_info$filename[i]
  ]

  if (length(hits) == 0L) {
    stop(
      "Could not find file: ",
      sample_info$filename[i],
      "\nSearch root: ", project_root
    )
  }

  if (length(hits) > 1L) {
    stop(
      "More than one file matched: ",
      sample_info$filename[i],
      "\n",
      paste(hits, collapse = "\n")
    )
  }

  sample_info$file_path[i] <- hits[1]
}

fwrite(
  sample_info,
  file.path(out_dir, "sample_file_manifest.csv"),
  bom = TRUE
)


result_list <- vector("list", nrow(sample_info))

for (i in seq_len(nrow(sample_info))) {

  sid <- sample_info$sample_id[i]
  grp <- sample_info$group[i]
  f <- sample_info$file_path[i]

  message(
    "\n[", i, "/", nrow(sample_info),
    "] ", grp, " | ", sid
  )

  dt <- read_cord_mat_mean(
    file_path = f,
    analysis_min = analysis_min,
    analysis_max = analysis_max
  )

  # Step 1: suppress the common one-channel periodic ripple.
  mean_deripple <- deripple_adjacent(dt$mean_intensity)

  # Step 2: retain the same light SG(5,2) smoothing used in v6.
  mean_sg <- sg_smooth(
    mean_deripple,
    window_points = sg_window_points
  )

  # Step 3: estimate and subtract the slowly varying baseline.
  bc <- arpls_baseline(
    y = mean_sg,
    lambda = baseline_lambda,
    ratio = baseline_ratio,
    max_iter = baseline_max_iter
  )

  dt[, `:=`(
    sample_id = sid,
    group = grp,
    mean_intensity_deripple = mean_deripple,
    mean_intensity_SG = mean_sg,
    estimated_baseline = bc$baseline,
    corrected_intensity = bc$corrected,
    sg_window_points = sg_window_points,
    arpls_iterations = bc$iterations,
    arpls_lambda = baseline_lambda
  )]

  result_list[[i]] <- dt

  rm(dt, bc)
  gc()
}

all_spectra <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

all_spectra[, group := factor(
  group,
  levels = group_order
)]

setorder(all_spectra, group, wavenumber)

fwrite(
  all_spectra,
  file.path(
    out_dir,
    "cord_mean_spectra_original_deripple_SG_baseline_corrected_300_1800.csv"
  ),
  bom = TRUE
)


p_final <- ggplot(
  all_spectra,
  aes(
    x = wavenumber,
    y = corrected_intensity,
    color = group,
    group = group
  )
) +
  geom_line(
    linewidth = line_width,
    alpha = 0.98,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = group_colors[group_order],
    breaks = group_order,
    limits = group_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(analysis_min, analysis_max),
    breaks = seq(300, 1800, by = 300),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = expression("Raman shift (cm"^{-1}*")"),
    y = "Intensity (a.u.)",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.y = if (hide_y_tick_labels) {
      element_blank()
    } else {
      element_text(size = 10, color = "black")
    },
    axis.ticks.y = if (hide_y_tick_labels) {
      element_blank()
    } else {
      element_line()
    },
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.width = grid::unit(1.10, "cm"),
    legend.text = element_text(
      face = "bold",
      size = 9.5
    ),
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.15,
        alpha = 1
      )
    )
  )

ggsave(
  file.path(
    out_dir,
    "Cord_mean_Raman_spectra_DERIPPLE_SG5_arPLS_300_1800.png"
  ),
  p_final,
  width = 10.5,
  height = 6.2,
  dpi = export_png_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    out_dir,
    "Cord_mean_Raman_spectra_DERIPPLE_SG5_arPLS_300_1800.pdf"
  ),
  p_final,
  width = 10.5,
  height = 6.2,
  bg = "white",
  limitsize = FALSE
)

message("\nDONE")
message("Output folder: ", out_dir)
message(
  "Final figure: ",
  file.path(
    out_dir,
    "Cord_mean_Raman_spectra_DERIPPLE_SG5_arPLS_300_1800.png"
  )
)
message("Periodic ripple suppression = adjacent-channel [1,2,1]/4")
message("SG smoothing = SG(", sg_window_points, ",2)")
message("arPLS lambda = ", baseline_lambda)
