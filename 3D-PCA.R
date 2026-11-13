rm(list = ls())
gc()

# 0. Packages
required_pkgs <- c("data.table", "plotly", "htmlwidgets")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install missing packages first:\n",
    paste(missing_pkgs, collapse = ", "),
    "\n\nRun:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
  library(htmlwidgets)
})


# 1. Paths

project_root <- "G:/ML"
results_root <- file.path(
  project_root,
  "results_pairwise_PCA_2D_placenta_300_1800_DISEASE_ALIGNED"
)
pairwise_root <- file.path(results_root, "05_pairwise_PCA")

fig_out_dir <- file.path(
  results_root,
  "07_pairwise_3D_PCA_Fig2B_no_legend_zoom_rotation"
)
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(pairwise_root)) {
  stop("Cannot find pairwise PCA results folder: ", pairwise_root)
}

# 2. Settings copied from the uploaded Fig2B-like code

group_order <- c("Control", "SLE", "APS", "sPE", "FGR", "ICP", "GDM")
disease_groups <- setdiff(group_order, "Control")

# Original colors from Fig2B-like code
group_colors <- c(
  Control = "#3E5C8A",   # muted blue
  SLE     = "#D74F8B",   # teal
  APS     = "#6BC48F",   # green
  sPE     = "#B48AD7",   # lavender
  FGR     = "#FFB347",   # orange-yellow
  ICP     = "#FF8C5A",   # orange
  GDM     = "#1FB0B0"    # coral-red
)

# Same display density as the original Fig2B-like code
max_pixels_per_group <- 2000L

# Original 3D appearance
marker_size_pairwise <- 2.4
marker_opacity <- 0.95
camera_eye <- list(x = 1.65, y = 1.40, z = 0.90)
scene_bg <- "rgb(245,245,245)"
scene_grid <- "white"
scene_zero <- "rgb(210,210,210)"
paper_bg <- "white"
axis_line_col <- "rgb(150,150,150)"
html_width <- 880
html_height <- 520

show_camera_panel <- TRUE
zoom_min_percent <- 50L
zoom_max_percent <- 250L
zoom_step_percent <- 5L
camera_panel_top_px <- 10L
camera_panel_right_px <- 10L

# Original code had optional PNG export disabled
try_export_png <- FALSE

set.seed(20260610)

# 3. Helper functions

safe_group_name <- function(x) {
  x <- as.character(x)
  x_low <- tolower(x)

  ifelse(
    x_low %in% c("control", "ctrl", "con"), "Control",
    ifelse(
      x_low == "sle", "SLE",
      ifelse(
        x_low == "aps", "APS",
        ifelse(
          x_low %in% c("spe", "s_pe", "s-pe", "severe_pe", "severepe"), "sPE",
          ifelse(
            x_low == "fgr", "FGR",
            ifelse(
              x_low == "icp", "ICP",
              ifelse(x_low == "gdm", "GDM", x)
            )
          )
        )
      )
    )
  )
}

sample_rows <- function(dt, max_n) {
  if (nrow(dt) <= max_n) return(copy(dt))
  dt[sample.int(nrow(dt), max_n, replace = FALSE)]
}

read_pair_existing_scores <- function(pair_id, disease_group) {
  pair_dir <- file.path(pairwise_root, pair_id)
  score_dir <- file.path(pair_dir, "04_pixel_scores")
  variance_file <- file.path(
    pair_dir,
    "02_tables",
    "PCA_PC1_PC3_variance_explained.csv"
  )

  if (!dir.exists(score_dir)) {
    stop("Cannot find score folder: ", score_dir)
  }
  if (!file.exists(variance_file)) {
    stop("Cannot find variance table: ", variance_file)
  }

  score_files <- list.files(score_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(score_files) == 0L) {
    stop("No score RDS files found in: ", score_dir)
  }

  score_list <- lapply(score_files, function(f) {
    dt <- as.data.table(readRDS(f))
    required_cols <- c(
      "group",
      "PC1_disease_aligned",
      "PC2_disease_aligned",
      "PC3_disease_aligned"
    )
    missing_cols <- setdiff(required_cols, names(dt))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing required columns in ", f, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }

    dt[, group := safe_group_name(group)]
    dt[, .(
      group,
      PC1 = PC1_disease_aligned,
      PC2 = PC2_disease_aligned,
      PC3 = PC3_disease_aligned
    )]
  })

  scores <- rbindlist(score_list, fill = TRUE)
  scores <- scores[group %in% c("Control", disease_group)]

  if (!all(c("Control", disease_group) %in% unique(scores$group))) {
    stop("Both Control and ", disease_group, " were not found in ", pair_id)
  }

  plot_scores <- rbindlist(list(
    sample_rows(scores[group == "Control"], max_pixels_per_group),
    sample_rows(scores[group == disease_group], max_pixels_per_group)
  ))

  plot_scores[, group := factor(group, levels = c("Control", disease_group))]

  variance_dt <- fread(variance_file)
  if (!all(c("PC", "variance_explained_percent") %in% names(variance_dt))) {
    stop("Variance table lacks PC or variance_explained_percent: ", variance_file)
  }

  var_exp <- variance_dt[match(c("PC1", "PC2", "PC3"), PC), variance_explained_percent]
  if (length(var_exp) != 3L || any(!is.finite(var_exp))) {
    stop("Could not read valid PC1-PC3 variance percentages from: ", variance_file)
  }

  list(scores = plot_scores, var_exp = var_exp)
}

add_visible_camera_panel <- function(
  p,
  initial_eye,
  min_percent = 50L,
  max_percent = 250L,
  step_percent = 5L,
  top_px = 10L,
  right_px = 10L
) {
  js_code <- sprintf(
"function(el, x) {
  var initialEye = {x: %.12f, y: %.12f, z: %.12f};
  var initialDistance = Math.sqrt(
    initialEye.x * initialEye.x +
    initialEye.y * initialEye.y +
    initialEye.z * initialEye.z
  );

  el.style.position = 'relative';

  var panel = document.createElement('div');
  panel.className = 'pca-camera-panel';
  panel.style.position = 'absolute';
  panel.style.top = '%dpx';
  panel.style.right = '%dpx';
  panel.style.zIndex = '1000';
  panel.style.padding = '8px 10px';
  panel.style.background = 'rgba(255,255,255,0.90)';
  panel.style.border = '1px solid rgba(120,120,120,0.45)';
  panel.style.borderRadius = '4px';
  panel.style.fontFamily = 'Arial, sans-serif';
  panel.style.fontSize = '12px';
  panel.style.lineHeight = '1.35';
  panel.style.color = '#222';
  panel.style.boxShadow = '0 1px 4px rgba(0,0,0,0.12)';
  panel.style.userSelect = 'none';
  panel.style.minWidth = '190px';

  var title = document.createElement('div');
  title.textContent = '3D camera';
  title.style.fontWeight = 'bold';
  title.style.marginBottom = '5px';
  panel.appendChild(title);

  function addReadout(labelText, initialText) {
    var row = document.createElement('div');
    row.style.display = 'flex';
    row.style.justifyContent = 'space-between';
    row.style.gap = '14px';

    var label = document.createElement('span');
    label.textContent = labelText;

    var value = document.createElement('span');
    value.textContent = initialText;
    value.style.fontWeight = 'bold';
    value.style.fontVariantNumeric = 'tabular-nums';

    row.appendChild(label);
    row.appendChild(value);
    panel.appendChild(row);
    return value;
  }

  var zoomValue = addReadout('Zoom', '100%%');
  var azimuthValue = addReadout('Azimuth', '0.0 deg');
  var elevationValue = addReadout('Elevation', '0.0 deg');

  var controls = document.createElement('div');
  controls.style.display = 'flex';
  controls.style.alignItems = 'center';
  controls.style.gap = '6px';
  controls.style.marginTop = '6px';

  var slider = document.createElement('input');
  slider.type = 'range';
  slider.min = '%d';
  slider.max = '%d';
  slider.step = '%d';
  slider.value = '100';
  slider.style.width = '125px';

  var resetButton = document.createElement('button');
  resetButton.type = 'button';
  resetButton.textContent = 'Reset';
  resetButton.title = 'Reset zoom and rotation to the initial camera';
  resetButton.style.padding = '2px 6px';
  resetButton.style.fontSize = '11px';
  resetButton.style.cursor = 'pointer';
  resetButton.style.border = '1px solid rgba(110,110,110,0.55)';
  resetButton.style.borderRadius = '3px';
  resetButton.style.background = 'white';

  controls.appendChild(slider);
  controls.appendChild(resetButton);
  panel.appendChild(controls);
  el.appendChild(panel);

  function getEye() {
    var scene = el._fullLayout && el._fullLayout.scene;
    var camera = scene && scene.camera;
    var eye = camera && camera.eye;
    if (!eye) return {x: initialEye.x, y: initialEye.y, z: initialEye.z};
    return {x: Number(eye.x), y: Number(eye.y), z: Number(eye.z)};
  }

  function eyeDistance(eye) {
    return Math.sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z);
  }

  function cameraMetrics(eye) {
    var distance = eyeDistance(eye);
    var zoom = (!isFinite(distance) || distance <= 0)
      ? 100
      : initialDistance / distance * 100;

    var azimuth = Math.atan2(eye.y, eye.x) * 180 / Math.PI;
    if (azimuth < 0) azimuth += 360;

    var horizontalDistance = Math.sqrt(eye.x * eye.x + eye.y * eye.y);
    var elevation = Math.atan2(eye.z, horizontalDistance) * 180 / Math.PI;

    return {zoom: zoom, azimuth: azimuth, elevation: elevation};
  }

  function updateDisplayFromCamera() {
    var metrics = cameraMetrics(getEye());
    var boundedZoom = Math.max(%d, Math.min(%d, metrics.zoom));

    zoomValue.textContent = Math.round(metrics.zoom) + '%%';
    azimuthValue.textContent = metrics.azimuth.toFixed(1) + ' deg';
    elevationValue.textContent = metrics.elevation.toFixed(1) + ' deg';
    slider.value = String(Math.round(boundedZoom / %d) * %d);
  }

  function setZoom(percent) {
    var eye = getEye();
    var distance = eyeDistance(eye);
    if (!isFinite(distance) || distance <= 0) {
      eye = {x: initialEye.x, y: initialEye.y, z: initialEye.z};
      distance = initialDistance;
    }

    var targetDistance = initialDistance * 100 / percent;
    var scale = targetDistance / distance;

    Plotly.relayout(el, {
      'scene.camera.eye.x': eye.x * scale,
      'scene.camera.eye.y': eye.y * scale,
      'scene.camera.eye.z': eye.z * scale
    });
  }

  slider.addEventListener('input', function() {
    setZoom(Number(slider.value));
  });

  resetButton.addEventListener('click', function() {
    Plotly.relayout(el, {
      'scene.camera.eye.x': initialEye.x,
      'scene.camera.eye.y': initialEye.y,
      'scene.camera.eye.z': initialEye.z
    });
  });

  el.on('plotly_relayout', function(eventData) {
    if (
      eventData['scene.camera'] ||
      eventData['scene.camera.eye.x'] !== undefined ||
      eventData['scene.camera.eye.y'] !== undefined ||
      eventData['scene.camera.eye.z'] !== undefined
    ) {
      window.requestAnimationFrame(updateDisplayFromCamera);
    }
  });

  window.requestAnimationFrame(updateDisplayFromCamera);
}",
    as.numeric(initial_eye$x),
    as.numeric(initial_eye$y),
    as.numeric(initial_eye$z),
    as.integer(top_px),
    as.integer(right_px),
    as.integer(min_percent),
    as.integer(max_percent),
    as.integer(step_percent),
    as.integer(min_percent),
    as.integer(max_percent),
    as.integer(step_percent),
    as.integer(step_percent)
  )

  htmlwidgets::onRender(p, js_code)
}

make_scene_axis <- function(title_text) {
  list(
    title = list(text = title_text, font = list(size = 16, color = "black")),
    showbackground = TRUE,
    backgroundcolor = scene_bg,
    gridcolor = scene_grid,
    zerolinecolor = scene_zero,
    showline = TRUE,
    linecolor = axis_line_col,
    tickfont = list(size = 16, color = "black"),
    showspikes = FALSE
  )
}

make_pca_plot_3d <- function(plot_dt, var_exp, plot_title, marker_size = 3.0) {
  p <- plot_ly()

  groups_here <- unique(as.character(plot_dt$group))
  groups_here <- group_order[group_order %in% groups_here]

  for (g in groups_here) {
    dtg <- plot_dt[group == g]

    p <- add_trace(
      p,
      data = dtg,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      type = "scatter3d",
      mode = "markers",
      name = g,
      showlegend = FALSE,
      marker = list(
        size = marker_size,
        color = unname(group_colors[g]),
        opacity = marker_opacity,
        line = list(color = "rgba(80,80,80,0.15)", width = 0.05)
      ),
      hovertemplate = paste0(
        "<b>", g, "</b><br>",
        "PC1: %{x:.2f}<br>",
        "PC2: %{y:.2f}<br>",
        "PC3: %{z:.2f}<extra></extra>"
      )
    )
  }

  p <- layout(
    p,
    title = list(
      text = plot_title,
      x = 0.03,
      xanchor = "left",
      font = list(size = 18, color = "black")
    ),
    showlegend = FALSE,
    margin = list(l = 0, r = 0, b = 0, t = 40),
    paper_bgcolor = paper_bg,
    scene = list(
      xaxis = make_scene_axis(""),
      yaxis = make_scene_axis(""),
      zaxis = make_scene_axis(""),
      aspectmode = "cube",
      camera = list(eye = camera_eye)
    )
  )

  # Enable mouse-wheel zoom and keep the original Fig2B appearance.
  p <- plotly::config(
    p,
    scrollZoom = TRUE,
    displaylogo = FALSE
  )

  if (isTRUE(show_camera_panel)) {
    p <- add_visible_camera_panel(
      p = p,
      initial_eye = camera_eye,
      min_percent = zoom_min_percent,
      max_percent = zoom_max_percent,
      step_percent = zoom_step_percent,
      top_px = camera_panel_top_px,
      right_px = camera_panel_right_px
    )
  }

  p
}

save_plotly_outputs <- function(p, html_file, png_file = NULL) {
  htmlwidgets::saveWidget(
    widget = p,
    file = html_file,
    selfcontained = TRUE,
    title = basename(html_file)
  )

  if (isTRUE(try_export_png) && !is.null(png_file)) {
    ok <- FALSE
    try({
      plotly::save_image(
        p,
        file = png_file,
        width = html_width,
        height = html_height
      )
      ok <- TRUE
    }, silent = TRUE)

    if (!ok) {
      message(
        "PNG export skipped or failed for: ", basename(html_file),
        "\nHTML was still saved successfully."
      )
    }
  }
}

# 4. Six pairwise 3D PCA plots from existing results

pairwise_summary_list <- list()

for (dg in disease_groups) {
  pair_id <- paste0("Control_vs_", dg)
  message("Reading existing PCA scores and drawing: ", pair_id)

  pair_res <- read_pair_existing_scores(pair_id, dg)
  plot_dt <- pair_res$scores
  var_exp <- pair_res$var_exp

  html_file <- file.path(fig_out_dir, paste0("PCA3D_Control_vs_", dg, ".html"))
  png_file <- file.path(fig_out_dir, paste0("PCA3D_Control_vs_", dg, ".png"))

  p <- make_pca_plot_3d(
    plot_dt = plot_dt,
    var_exp = var_exp,
    plot_title = paste0("3D PCA: Control vs ", dg),
    marker_size = marker_size_pairwise
  )

  save_plotly_outputs(p, html_file, png_file)

  fwrite(
    plot_dt,
    file.path(fig_out_dir, paste0("PCA3D_Control_vs_", dg, "_plotted_scores.csv"))
  )

  pairwise_summary_list[[dg]] <- data.table(
    comparison = pair_id,
    PC1_percent = var_exp[1],
    PC2_percent = var_exp[2],
    PC3_percent = var_exp[3],
    n_Control_plotted = sum(plot_dt$group == "Control"),
    n_Disease_plotted = sum(plot_dt$group == dg),
    score_source = "existing disease-aligned PCA scores",
    PCA_refitted = FALSE
  )
}

pairwise_summary <- rbindlist(pairwise_summary_list)
fwrite(
  pairwise_summary,
  file.path(fig_out_dir, "pairwise_PCA3D_existing_results_summary.csv")
)

message("Done.")
message("All outputs saved under:\n", fig_out_dir)
message("The user-adjusted Fig2B-like font and marker settings were retained; live zoom and rotation angles were added.")
