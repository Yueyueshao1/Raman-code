rm(list = ls())
gc()

# 0. R packages
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

# 1. USER SETTINGS
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
    "A胎盘2D_layer1.mat",
    "B胎盘2D_layer1.mat",
    "C胎盘2D_layer1.mat",
    "D胎盘2D_layer1.mat",
    "E胎盘2D_layer1.mat",
    "F胎盘2D_layer1.mat",
    "G胎盘2D_layer1.mat"
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
baseline_lambda <- 1e4
baseline_ratio <- 1e-6
baseline_max_iter <- 100L

# Figure settings
line_width <- 0.72
export_png_dpi <- 1200
hide_y_tick_labels <- TRUE

out_dir <- file.path(
  project_root,
  "results_placenta_2D_mean_spectra_arPLS_300_1800_v4_SCIPY"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. PYTHON / SCIPY MAT READER

if ("py_require" %in% getNamespaceExports("reticulate")) {
  reticulate::py_require(c("numpy", "scipy"))
}

# Trigger Python and check scipy availability.
if (!reticulate::py_module_available("numpy") ||
    !reticulate::py_module_available("scipy")) {
  stop(
    "\nPython numpy/scipy is not available.\n\n",
    "Run ONCE in R:\n",
    "  reticulate::py_install(c('numpy','scipy'), pip = TRUE)\n\n",
    "Then restart R and rerun this script.\n"
  )
}

reticulate::py_run_string("
import numpy as np
import scipy.io as sio

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

def read_placenta_mat_numeric(path):
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

    # Find a plausible spectral axis whose length matches one X dimension.
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

    # Prefer the longest matching vector.
    wn = sorted(candidates, key=lambda a: a.size, reverse=True)[0]

    # Orient X as pixels x Raman channels.
    if X.shape[1] == wn.size:
        pass
    elif X.shape[0] == wn.size:
        X = X.T
    else:
        raise RuntimeError(
            'Neither dimension of data matches Raman axis length'
        )

    return {
        'X': X,
        'wn': wn
    }
")

read_placenta_mat_mean <- function(
  file_path,
  analysis_min = 300,
  analysis_max = 1800
) {

  message("Reading: ", file_path)

  # Use an ASCII-only temporary pathname. This avoids any Windows path
  # encoding problem even though scipy itself can generally handle Unicode.
  tmp_mat <- tempfile(pattern = "placenta2D_", fileext = ".mat")
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
    reticulate::py$read_placenta_mat_numeric(py_path),
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

# ------------------------------------------------------------
# 3. arPLS BASELINE
# ------------------------------------------------------------
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

# 4. FIND THE 7 MAT FILES
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

# 5. 7 MEAN SPECTRA + arPLS
result_list <- vector("list", nrow(sample_info))

for (i in seq_len(nrow(sample_info))) {

  sid <- sample_info$sample_id[i]
  grp <- sample_info$group[i]
  f <- sample_info$file_path[i]

  message(
    "\n[", i, "/", nrow(sample_info),
    "] ", grp, " | ", sid
  )

  dt <- read_placenta_mat_mean(
    file_path = f,
    analysis_min = analysis_min,
    analysis_max = analysis_max
  )

  bc <- arpls_baseline(
    y = dt$mean_intensity,
    lambda = baseline_lambda,
    ratio = baseline_ratio,
    max_iter = baseline_max_iter
  )

  dt[, `:=`(
    sample_id = sid,
    group = grp,
    estimated_baseline = bc$baseline,
    corrected_intensity = bc$corrected,
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
    "placenta_2D_mean_spectra_original_baseline_corrected_300_1800.csv"
  ),
  bom = TRUE
)

# 6. FINAL FIGURE
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
    "Placenta_2D_mean_Raman_spectra_arPLS_300_1800.png"
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
    "Placenta_2D_mean_Raman_spectra_arPLS_300_1800.pdf"
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
    "Placenta_2D_mean_Raman_spectra_arPLS_300_1800.png"
  )
)
message("arPLS lambda = ", baseline_lambda)
