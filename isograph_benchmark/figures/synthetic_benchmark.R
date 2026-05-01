#!/usr/bin/env Rscript
# Publication-quality figures for the IsoGraph synthetic benchmark.
# Requires: ggplot2, patchwork, dplyr, tidyr, arrow, scales

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---------------------------------------------------------------------------
# Paths - resolve project root via .here marker
# ---------------------------------------------------------------------------
find_root <- function() {
  d <- normalizePath(getwd(), mustWork = FALSE)
  while (nchar(d) > 1) {
    if (file.exists(file.path(d, ".here"))) return(d)
    d <- dirname(d)
  }
  getwd()
}
ROOT <- find_root()
rel <- function(...) file.path(ROOT, ...)

RAW_PATH     <- rel("benchmark", "01_synthetic", "_m", "synthetic_results.parquet")
SUMMARY_PATH <- rel("benchmark", "02_metrics", "_m", "synthetic_metric_summary.parquet")
FIG_DIR      <- rel("benchmark", "02_metrics", "figures")
TABLE_DIR    <- rel("benchmark", "02_metrics", "_m")
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------
METHOD_ORDER <- c(
  "isograph_baseline",
  "isograph_latent",
  "isograph_graph",
  "isograph_cpu_latent",
  "isograph_gpu_latent",
  "isograph_vae",
  "isograph_vae_gpu",
  "wgcna_gene"
)
SCALE_METHOD_ORDER <- c(
  "isograph_cpu_latent",
  "isograph_gpu_latent",
  "isograph_vae",
  "isograph_vae_gpu",
  "wgcna_gene"
)

METHOD_LABELS <- c(
  isograph_baseline   = "IsoGraph Baseline",
  isograph_latent     = "IsoGraph Latent",
  isograph_graph      = "IsoGraph Graph",
  isograph_cpu_latent = "IsoGraph CPU Latent",
  isograph_gpu_latent = "IsoGraph GPU Latent",
  isograph_vae        = "IsoGraph VAE",
  isograph_vae_gpu    = "IsoGraph VAE GPU",
  wgcna_gene          = "WGCNA"
)

COMPUTE_LABELS <- c(
  isograph_baseline   = "CPU",
  isograph_latent     = "CPU",
  isograph_graph      = "CPU",
  isograph_cpu_latent = "CPU",
  isograph_gpu_latent = "GPU",
  isograph_vae        = "CPU",
  isograph_vae_gpu    = "GPU",
  wgcna_gene          = "CPU"
)

# Okabe-Ito colorblind-safe palette
METHOD_COLORS <- c(
  isograph_baseline   = "#0072B2",
  isograph_latent     = "#E69F00",
  isograph_graph      = "#009E73",
  isograph_cpu_latent = "#56B4E9",
  isograph_gpu_latent = "#F0E442",
  isograph_vae        = "#D55E00",
  isograph_vae_gpu    = "#CC79A7",
  wgcna_gene          = "#666666"
)

SCENARIO_ORDER <- c(
  "idealized_switching",
  "noise_stress",
  "feature_space_interactions",
  "non_switching_background",
  "unequal_isoform_abundance"
)

SCENARIO_LABELS <- c(
  idealized_switching        = "Idealized\nSwitching",
  noise_stress               = "Noise\nStress",
  feature_space_interactions = "Feature\nInteractions",
  non_switching_background   = "Non-Switching\nBackground",
  unequal_isoform_abundance  = "Unequal\nAbundance"
)

METRIC_LABELS <- c(
  metrics_module_recovery           = "Module recovery [0-1]",
  metrics_switch_gene_detection_rate = "Switch gene detection rate [0-1]",
  metrics_nonswitch_gene_module_rate = "Non-switching gene\nmodule rate [0-1]"
)

# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(size = 7),
      axis.title       = element_text(size = 8),
      legend.text      = element_text(size = 7),
      legend.title     = element_blank(),
      legend.key.size  = unit(0.35, "cm"),
      strip.text       = element_text(size = 7.5, face = "plain"),
      strip.background = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin      = margin(4, 8, 4, 4, "pt")
    )
}

scale_color_method <- function()
  scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS[METHOD_ORDER],
                     breaks = METHOD_ORDER)
scale_fill_method <- function()
  scale_fill_manual(values = METHOD_COLORS, labels = METHOD_LABELS[METHOD_ORDER],
                    breaks = METHOD_ORDER)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
load_raw <- function(path, method_order = METHOD_ORDER) {
  df <- read_parquet(path) |>
    filter(status == "completed", run_method %in% method_order) |>
    mutate(run_method = factor(run_method, levels = method_order))
  df
}

load_summary <- function(path, method_order = METHOD_ORDER) {
  df <- read_parquet(path) |>
    filter(method %in% method_order) |>
    mutate(method = factor(method, levels = method_order))
  df
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
method_scale <- function(df_col = "run_method") {
  list(
    scale_color_manual(
      values = METHOD_COLORS,
      labels = METHOD_LABELS[METHOD_ORDER],
      breaks = METHOD_ORDER
    ),
    scale_fill_manual(
      values = METHOD_COLORS,
      labels = METHOD_LABELS[METHOD_ORDER],
      breaks = METHOD_ORDER
    )
  )
}

dot_ci_panel <- function(summary_df, metric_name, metric_label,
                          scenarios = SCENARIO_ORDER) {
  sub <- summary_df |>
    filter(
      .data$metric == metric_name,
      .data$scenario %in% scenarios
    ) |>
    mutate(
      scenario = factor(.data$scenario, levels = scenarios,
                        labels = SCENARIO_LABELS[scenarios]),
      method   = factor(.data$method,   levels = METHOD_ORDER)
    )

  ggplot(sub, aes(x = .data$method, y = .data$mean, color = .data$method)) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed",
               color = "grey60") +
    geom_errorbar(aes(ymin = .data$ci_low, ymax = .data$ci_high),
                  width = 0.25, linewidth = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~ scenario, nrow = 1) +
    scale_x_discrete(labels = METHOD_LABELS[METHOD_ORDER]) +
    scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS[METHOD_ORDER],
                       breaks = METHOD_ORDER) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = NULL, y = metric_label, color = NULL) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
      legend.position = "none"
    )
}

runtime_violin <- function(raw_df) {
  ggplot(raw_df, aes(x = .data$run_method, y = .data$measurement_elapsed_sec,
                     fill = .data$run_method)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.08, outlier.size = 0.5, outlier.alpha = 0.3,
                 linewidth = 0.5, fill = "white", color = "grey30") +
    scale_x_discrete(labels = METHOD_LABELS[METHOD_ORDER]) +
    scale_fill_manual(values = METHOD_COLORS,
                      labels = METHOD_LABELS[METHOD_ORDER]) +
    scale_y_log10(labels = label_comma()) +
    labs(x = NULL, y = "Runtime (s, log10 scale)", fill = NULL) +
    theme_pub() +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
}

# ---------------------------------------------------------------------------
# Main figure (Figure 1)
# ---------------------------------------------------------------------------
make_fig1 <- function(raw, summary) {
  scenarios <- SCENARIO_ORDER[SCENARIO_ORDER %in% unique(summary$scenario)]
  methods   <- METHOD_ORDER[METHOD_ORDER %in% unique(summary$method)]

  pA <- dot_ci_panel(summary, "metrics_module_recovery",
                     METRIC_LABELS["metrics_module_recovery"], scenarios)
  pB <- dot_ci_panel(summary, "metrics_switch_gene_detection_rate",
                     METRIC_LABELS["metrics_switch_gene_detection_rate"], scenarios)
  pC <- dot_ci_panel(summary, "metrics_nonswitch_gene_module_rate",
                     METRIC_LABELS["metrics_nonswitch_gene_module_rate"], scenarios)
  pD <- runtime_violin(raw)

  pA <- pA + labs(tag = "A")
  pB <- pB + labs(tag = "B")
  pC <- pC + labs(tag = "C")
  pD <- pD +
    labs(tag = "D") +
    theme(legend.position = "right")

  fig <- (pA / pB / pC / pD) +
    plot_layout(heights = c(1, 1, 1, 0.75), guides = "collect") &
    theme(
      plot.tag = element_text(size = 10, face = "bold"),
      legend.position = "right"
    )
  fig
}

# ---------------------------------------------------------------------------
# Heatmap helper for S1-S3
# ---------------------------------------------------------------------------
heatmap_fig <- function(raw, scenario, metric, row_col, col_col,
                        row_label, col_label,
                        method_order = METHOD_ORDER) {
  sub <- raw |>
    filter(run_scenario == scenario, !is.na(.data[[metric]])) |>
    group_by(run_method, .data[[row_col]], .data[[col_col]]) |>
    summarise(mean_val = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") |>
    filter(run_method %in% method_order) |>
    mutate(run_method = factor(run_method, levels = method_order,
                               labels = METHOD_LABELS[method_order]))

  if (nrow(sub) == 0) return(NULL)

  ggplot(sub, aes(x = factor(.data[[col_col]]),
                  y = factor(.data[[row_col]]),
                  fill = mean_val)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", mean_val)),
              size = 2.2, color = ifelse(sub$mean_val < 0.55, "white", "black")) +
    facet_wrap(~ run_method, nrow = 2) +
    scale_fill_viridis_c(name = "Module recovery\n(mean)", limits = c(0, 1),
                         option = "D") +
    labs(x = col_label, y = row_label) +
    theme_pub() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# ---------------------------------------------------------------------------
# Line plot helper for S4, S5, S6
# ---------------------------------------------------------------------------
line_ci_fig <- function(raw, scenario, param_col, param_label,
                        metrics, metric_labels,
                        method_order = METHOD_ORDER,
                        y_log = FALSE) {
  sub <- raw |>
    filter(run_scenario == scenario,
           run_method %in% method_order,
           !is.na(.data[[param_col]])) |>
    pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") |>
    filter(!is.na(value)) |>
    group_by(run_method, .data[[param_col]], metric) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      se_val   = if (n() > 1) sd(value, na.rm = TRUE) / sqrt(n()) else 0,
      .groups  = "drop"
    ) |>
    mutate(
      metric     = factor(metric, levels = metrics, labels = metric_labels),
      run_method = factor(run_method, levels = method_order)
    )

  if (nrow(sub) == 0) return(NULL)

  p <- ggplot(sub, aes(x = .data[[param_col]], y = mean_val,
                       color = run_method, fill = run_method)) +
    geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = METHOD_COLORS[method_order],
                       labels = METHOD_LABELS[method_order]) +
    scale_fill_manual(values  = METHOD_COLORS[method_order],
                      labels = METHOD_LABELS[method_order]) +
    labs(x = param_label, y = NULL, color = NULL, fill = NULL) +
    theme_pub() +
    theme(legend.position = "right")

  if (!y_log) {
    p <- p + scale_y_continuous(limits = c(0, NA),
                                 expand = expansion(mult = c(0.02, 0.05)))
  } else {
    p <- p + scale_y_log10(labels = label_comma())
  }
  p
}

# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------
save_fig <- function(p, name, width, height) {
  pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(pdf_path, plot = p, width = width, height = height, units = "in",
         device = cairo_pdf)
  ggsave(png_path, plot = p, width = width, height = height, units = "in",
         dpi = 300)
  cat("  ", name, " saved\n", sep = "")
}

# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
make_table1 <- function(summary) {
  scenarios <- SCENARIO_ORDER[SCENARIO_ORDER %in% unique(summary$scenario)]
  methods   <- METHOD_ORDER[METHOD_ORDER %in% unique(summary$method)]

  metric_map <- c(
    metrics_module_recovery            = "Module recovery",
    metrics_switch_gene_detection_rate = "Switch detection",
    metrics_nonswitch_gene_module_rate = "False positive rate"
  )

  rows <- list()
  for (scen in scenarios) {
    for (meth in methods) {
      row <- list(
        Scenario = gsub("_", " ", scen) |> tools::toTitleCase(),
        Method   = METHOD_LABELS[[meth]],
        Compute  = COMPUTE_LABELS[[meth]]
      )
      for (metric in names(metric_map)) {
        sub <- summary |>
          filter(scenario == scen, method == meth, .data$metric == !!metric)
        col_name <- paste0(metric_map[[metric]], ", mean (95% CI)")
        if (nrow(sub) == 0) {
          row[[col_name]] <- "NA"
          row[["N"]] <- NA_integer_
        } else {
          row[["N"]] <- sub$n
          row[[col_name]] <- sprintf("%.3f [%.3f, %.3f]",
                                     sub$mean, sub$ci_low, sub$ci_high)
        }
      }
      rt <- summary |>
        filter(scenario == scen, method == meth, .data$metric == "measurement_elapsed_sec")
      row[["Runtime, s median"]] <- if (nrow(rt) > 0) sprintf("%.1f", rt$median) else "NA"
      rows[[length(rows) + 1]] <- as.data.frame(row, check.names = FALSE)
    }
  }

  tbl <- do.call(rbind, rows)
  out_path <- file.path(TABLE_DIR, "table1_benchmark_summary.csv")
  write.csv(tbl, out_path, row.names = FALSE)
  cat("  table1_benchmark_summary.csv:", nrow(tbl), "rows\n")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
cat("Loading data...\n")
if (!file.exists(RAW_PATH))     stop("Missing: ", RAW_PATH)
if (!file.exists(SUMMARY_PATH)) stop("Missing: ", SUMMARY_PATH)

raw     <- load_raw(RAW_PATH)
summary <- load_summary(SUMMARY_PATH)
cat("  Completed runs:", nrow(raw), " | Summary rows:", nrow(summary), "\n")

cat("Generating figures...\n")

# Fig 1 - Main overview
tryCatch({
  fig1 <- make_fig1(raw, summary)
  save_fig(fig1, "fig1_benchmark_overview", width = 14, height = 12)
}, error = function(e) warning("fig1 error: ", conditionMessage(e)))

# FigS1 - Idealized switching heatmap
tryCatch({
  p <- heatmap_fig(raw, "idealized_switching", "metrics_module_recovery",
                   "run_switching_fraction", "run_noise_sd",
                   "Switching fraction", "Noise SD")
  if (!is.null(p)) save_fig(p, "figS1_idealized_switching", width = 13, height = 8)
  else cat("  figS1: no data, skipping\n")
}, error = function(e) warning("figS1 error: ", conditionMessage(e)))

# FigS2 - Noise stress heatmap
tryCatch({
  p <- heatmap_fig(raw, "noise_stress", "metrics_module_recovery",
                   "run_count_dispersion", "run_noise_sd",
                   "Count dispersion", "Noise SD")
  if (!is.null(p)) save_fig(p, "figS2_noise_stress", width = 13, height = 8)
  else cat("  figS2: no data, skipping\n")
}, error = function(e) warning("figS2 error: ", conditionMessage(e)))

# FigS3 - Feature interactions heatmap
tryCatch({
  p <- heatmap_fig(raw, "feature_space_interactions", "metrics_module_recovery",
                   "run_interaction_strength", "run_interaction_fraction",
                   "Interaction strength", "Interaction fraction")
  if (!is.null(p)) save_fig(p, "figS3_feature_interactions", width = 13, height = 8)
  else cat("  figS3: no data, skipping\n")
}, error = function(e) warning("figS3 error: ", conditionMessage(e)))

# FigS4 - Non-switching background
tryCatch({
  p <- line_ci_fig(
    raw, "non_switching_background", "run_switching_fraction",
    "Background switching fraction",
    metrics = c("metrics_module_recovery", "metrics_nonswitch_gene_module_rate"),
    metric_labels = c("Module recovery [0-1]", "Non-switching gene module rate [0-1]")
  )
  if (!is.null(p)) save_fig(p, "figS4_nonswitching_background", width = 9, height = 3.8)
  else cat("  figS4: no data, skipping\n")
}, error = function(e) warning("figS4 error: ", conditionMessage(e)))

# FigS5 - Unequal abundance
tryCatch({
  p <- line_ci_fig(
    raw, "unequal_isoform_abundance", "run_abundance_imbalance",
    "Abundance imbalance ratio",
    metrics = c("metrics_module_recovery", "metrics_switch_gene_detection_rate"),
    metric_labels = c("Module recovery [0-1]", "Switch gene detection rate [0-1]")
  )
  if (!is.null(p)) save_fig(p, "figS5_unequal_abundance", width = 9, height = 3.8)
  else cat("  figS5: no data, skipping\n")
}, error = function(e) warning("figS5 error: ", conditionMessage(e)))

# FigS6 - Scale (n_genes)
tryCatch({
  scale_raw <- raw |> filter(run_method %in% SCALE_METHOD_ORDER)
  p <- line_ci_fig(
    scale_raw, "scale", "run_n_genes",
    "Number of genes",
    metrics = c("metrics_module_recovery", "measurement_elapsed_sec"),
    metric_labels = c("Module recovery [0-1]", "Runtime (s)"),
    method_order = SCALE_METHOD_ORDER,
    y_log = FALSE
  )
  if (!is.null(p)) {
    # Runtime panel should be log scale - rebuild with mixed scales
    sub_rt <- scale_raw |>
      filter(run_scenario == "scale", !is.na(run_n_genes)) |>
      group_by(run_method, run_n_genes) |>
      summarise(
        mean_val = mean(measurement_elapsed_sec, na.rm = TRUE),
        se_val   = if (n() > 1) sd(measurement_elapsed_sec, na.rm = TRUE) / sqrt(n()) else 0,
        .groups  = "drop"
      ) |>
      mutate(run_method = factor(run_method, levels = SCALE_METHOD_ORDER))

    sub_mr <- scale_raw |>
      filter(run_scenario == "scale", !is.na(run_n_genes)) |>
      group_by(run_method, run_n_genes) |>
      summarise(
        mean_val = mean(metrics_module_recovery, na.rm = TRUE),
        se_val   = if (n() > 1) sd(metrics_module_recovery, na.rm = TRUE) / sqrt(n()) else 0,
        .groups  = "drop"
      ) |>
      mutate(run_method = factor(run_method, levels = SCALE_METHOD_ORDER))

    p_mr <- ggplot(sub_mr, aes(run_n_genes, mean_val,
                                color = run_method, fill = run_method)) +
      geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_color_manual(values = METHOD_COLORS[SCALE_METHOD_ORDER],
                         labels = METHOD_LABELS[SCALE_METHOD_ORDER]) +
      scale_fill_manual(values  = METHOD_COLORS[SCALE_METHOD_ORDER],
                        labels = METHOD_LABELS[SCALE_METHOD_ORDER]) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
      labs(x = "Number of genes", y = "Module recovery [0-1]",
           color = NULL, fill = NULL) +
      theme_pub() + theme(legend.position = "none")

    p_rt <- ggplot(sub_rt, aes(run_n_genes, mean_val,
                                color = run_method, fill = run_method)) +
      geom_ribbon(aes(ymin = pmax(mean_val - se_val, 1),
                      ymax = mean_val + se_val),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_color_manual(values = METHOD_COLORS[SCALE_METHOD_ORDER],
                         labels = METHOD_LABELS[SCALE_METHOD_ORDER]) +
      scale_fill_manual(values  = METHOD_COLORS[SCALE_METHOD_ORDER],
                        labels = METHOD_LABELS[SCALE_METHOD_ORDER]) +
      scale_y_log10(labels = label_comma()) +
      labs(x = "Number of genes", y = "Runtime (s, log10 scale)",
           color = NULL, fill = NULL) +
      theme_pub() + theme(legend.position = "right")

    p_s6 <- p_mr + p_rt + plot_layout(guides = "collect") &
      theme(legend.position = "right")
    save_fig(p_s6, "figS6_scale", width = 9, height = 3.8)
  } else {
    cat("  figS6: no scale data, skipping\n")
  }
}, error = function(e) warning("figS6 error: ", conditionMessage(e)))

# Table 1
tryCatch(make_table1(summary), error = function(e) warning("table1 error: ", conditionMessage(e)))

cat("Done. Outputs in", FIG_DIR, "\n")
