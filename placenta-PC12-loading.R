# 0. User settings

project_root <- "G:/ML"

result_root <- file.path(
  project_root,
  "results_pairwise_PCA_2D_placenta_300_1800_DISEASE_ALIGNED_PROMINENT_EXTREMA"
)

pairwise_root <- file.path(result_root, "05_pairwise_PCA")

out_dir <- file.path(
  result_root,
  "08_PC1_PC2_SIX_DISEASES_STACKED_TWO_PANELS_GREY_BANDS"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Order of diseases (top to bottom)
disease_order <- c("SLE", "APS", "sPE", "FGR", "ICP", "GDM")
pair_ids <- paste0("Control_vs_", disease_order)

# Colors for the six disease curves
disease_colors <- c(
  SLE     = "#D74F8B",
  APS     = "#FE7882",
  sPE     = "#FF9360",
  FGR     = "#FFBE4A",
  ICP     = "#8A79AA",
  GDM     = "#1FB0B0"
)
# Display range (recommended)
display_min <- 900
display_max <- 1700

# If you want full-range QC, change to 300 and 1800.
# display_min <- 300
# display_max <- 1800

# Display style
use_sg7_for_display <- TRUE   # TRUE = SG7 smoothing for display only
max_abs_scale_for_display <- TRUE
curve_amplitude <- 0.46
offset_spacing <- 1.08

# Export settings
png_dpi <- 900
fig_width <- 8.4
fig_height <- 4.8

# Grey band appearance
band_fill <- "grey70"
band_alpha <- 0.32

# 1. User-defined characteristic bands (PC-specific)

selected_bands <- data.frame(
  PC = c("PC1", "PC1", "PC1", "PC2", "PC2", "PC2", "PC2"),
  xmin = c(
    998.8, 1151.1, 1509.8,
    1154.9, 1245.0, 1440.7, 1657.0
  ),
  xmax = c(
    1002.6, 1154.9, 1517.0,
    1177.5, 1286.0, 1448.0, 1664.1
  ),
  stringsAsFactors = FALSE
)

# 2. Packages
required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, dependencies = TRUE)
}

library(data.table)
library(ggplot2)

# 3. Helper function: SG7 smoother for display only
sg7_smooth <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 7L) return(x)

  coef <- c(-2, 3, 6, 7, 6, 3, -2) / 21
  y <- as.numeric(stats::filter(x, filter = coef, sides = 2))
  y[!is.finite(y)] <- x[!is.finite(y)]
  y
}

# 4. Read existing disease-aligned loadings
if (!dir.exists(pairwise_root)) {
  stop(
    "Cannot find pairwise result folder:\n", pairwise_root,
    "\nPlease modify project_root / result_root at the top of this script."
  )
}

all_loading <- list()

for (i in seq_along(disease_order)) {
  disease <- disease_order[i]
  pair_id <- pair_ids[i]

  loading_file <- file.path(
    pairwise_root,
    pair_id,
    "02_tables",
    "PCA_PC1_PC3_loadings_contributions_and_direction_audit.csv"
  )

  if (!file.exists(loading_file)) {
    stop("Missing loading file:\n", loading_file)
  }

  dt <- data.table::fread(loading_file)

  required_cols <- c(
    "wavenumber",
    "PC1_loading_disease_aligned",
    "PC2_loading_disease_aligned",
    "PC3_loading_disease_aligned"
  )
  if (!all(required_cols %in% names(dt))) {
    stop(
      "Required columns missing in:\n", loading_file,
      "\nExpected: ", paste(required_cols, collapse = ", ")
    )
  }

  long <- data.table::melt(
    dt,
    id.vars = "wavenumber",
    measure.vars = c(
      "PC1_loading_disease_aligned",
      "PC2_loading_disease_aligned",
      "PC3_loading_disease_aligned"
    ),
    variable.name = "PC",
    value.name = "loading_raw"
  )

  long[, PC := sub("_loading_disease_aligned$", "", PC)]
  long[, disease_group := disease]
  long[, pair_id := pair_id]
  long[, loading_sg7 := sg7_smooth(loading_raw), by = PC]

  all_loading[[pair_id]] <- long
}

loading_all <- data.table::rbindlist(all_loading, fill = TRUE)
loading_all <- loading_all[
  PC %in% c("PC1", "PC2") &
    wavenumber >= display_min &
    wavenumber <= display_max
]
loading_all[, disease_group := factor(disease_group, levels = disease_order)]
loading_all[, PC := factor(PC, levels = c("PC1", "PC2"))]

# 5. Build display data
make_display_data <- function(source_dt) {
  dt <- copy(source_dt)

  dt[, loading_display_source := if (isTRUE(use_sg7_for_display)) loading_sg7 else loading_raw]

  if (isTRUE(max_abs_scale_for_display)) {
    dt[, scale_denominator := max(abs(loading_display_source), na.rm = TRUE),
       by = .(disease_group, PC)]
    dt[!is.finite(scale_denominator) | scale_denominator == 0, scale_denominator := 1]
    dt[, loading_scaled := loading_display_source / scale_denominator]
  } else {
    dt[, scale_denominator := 1]
    dt[, loading_scaled := loading_display_source]
  }

  # top = SLE, bottom = GDM
  offset_table <- data.table(
    disease_group = factor(disease_order, levels = disease_order),
    offset = rev(seq(from = 0, by = offset_spacing, length.out = length(disease_order)))
  )

  dt <- merge(dt, offset_table, by = "disease_group", all.x = TRUE, sort = FALSE)
  dt[, y_plot := offset + curve_amplitude * loading_scaled]

  list(data = dt, offsets = offset_table)
}

plot_obj <- make_display_data(loading_all)
plot_dt <- plot_obj$data
offset_table <- plot_obj$offsets

# y-axis disease labels
axis_dt <- offset_table[order(offset)]
y_breaks <- axis_dt$offset
y_labels <- as.character(axis_dt$disease_group)

# grey band rectangles in the correct PC panel only
bands_dt <- as.data.table(selected_bands)
bands_dt[, PC := factor(PC, levels = c("PC1", "PC2"))]

# 6. Main figure
p <- ggplot(
  plot_dt,
  aes(
    x = wavenumber,
    y = y_plot,
    group = disease_group,
    color = disease_group
  )
) +
  geom_rect(
    data = bands_dt,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = band_fill,
    alpha = band_alpha
  ) +
  geom_hline(
    data = data.table::CJ(
      disease_group = factor(disease_order, levels = disease_order),
      PC = factor(c("PC1", "PC2"), levels = c("PC1", "PC2")),
      unique = TRUE
    )[, offset := setNames(offset_table$offset, offset_table$disease_group)[as.character(disease_group)]],
    aes(yintercept = offset),
    inherit.aes = FALSE,
    color = "grey86",
    linewidth = 0.28
  ) +
  geom_line(linewidth = 0.72, lineend = "round") +
  facet_grid(. ~ PC, scales = "fixed") +
  scale_color_manual(values = disease_colors, drop = FALSE) +
  scale_x_continuous(
    limits = c(display_min, display_max),
    breaks = seq(
      ceiling(display_min / 100) * 100,
      floor(display_max / 100) * 100,
      by = 100
    ),
    expand = expansion(mult = c(0.003, 0.01))
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    limits = c(min(offset_table$offset) - 0.45, max(offset_table$offset) + 0.45),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression("Wavenumber (cm"^{-1}*")"),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10.5, face = "bold", color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 8)),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.y = element_blank(),
    panel.spacing.x = grid::unit(0.6, "lines"),
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 8)
  )

# 7. Export
png_file <- file.path(out_dir, "Fig_PC1_PC2_SIX_DISEASES_STACKED_GREY_BANDS.png")
pdf_file <- file.path(out_dir, "Fig_PC1_PC2_SIX_DISEASES_STACKED_GREY_BANDS.pdf")

ggsave(
  png_file,
  p,
  width = fig_width,
  height = fig_height,
  dpi = png_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  pdf_file,
  p,
  width = fig_width,
  height = fig_height,
  bg = "white",
  limitsize = FALSE
)

# Also export the exact display data used for plotting.
data.table::fwrite(
  plot_dt,
  file.path(out_dir, "PLOT_DATA_PC1_PC2_SIX_DISEASES_STACKED_GREY_BANDS.csv")
)

# 8. README
readme <- c(
  "Placenta PC1-PC2 six-disease stacked loading figure with grey characteristic bands",
  "",
  "What this script does:",
  "- Reads existing disease-aligned loading tables from the pairwise PCA workflow.",
  "- Does NOT refit PCA.",
  "- Shows only PC1 and PC2.",
  "- Within each PC panel, all six disease groups are stacked vertically.",
  "- Grey vertical rectangles mark the user-defined characteristic Raman regions.",
  "- No peak-position numbers are drawn.",
  "- No legend is drawn.",
  "",
  "PC1 bands:",
  "- 998.8-1002.6 cm^-1",
  "- 1151.1-1154.9 cm^-1",
  "- 1509.8-1517.0 cm^-1",
  "",
  "PC2 bands:",
  "- 1154.9-1177.5 cm^-1",
  "- 1245.0-1286.0 cm^-1",
  "- 1440.7-1448.0 cm^-1",
  "- 1657.0-1664.1 cm^-1",
  "",
  "Display-only operations:",
  paste0("- SG7 smoothing: ", use_sg7_for_display),
  paste0("- Max-|loading| scaling: ", max_abs_scale_for_display),
  paste0("- Vertical offset spacing: ", offset_spacing),
  paste0("- Curve amplitude: ", curve_amplitude)
)

writeLines(readme, file.path(out_dir, "README_PC1_PC2_SIX_DISEASES_STACKED_GREY_BANDS.txt"))

message("============================================================")
message("PC1-PC2 six-disease stacked figure completed.")
message("Output folder: ", out_dir)
message("Main figure:")
message("  Fig_PC1_PC2_SIX_DISEASES_STACKED_GREY_BANDS")
message("============================================================")
