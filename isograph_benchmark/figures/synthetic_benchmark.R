#!/usr/bin/env Rscript
# Publication-quality figures for the IsoGraph synthetic benchmark.
# Requires: ggpubr, ggplot2, patchwork, dplyr, tidyr, arrow, scales

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
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
LONG_PATH    <- rel("benchmark", "02_metrics", "_m", "synthetic_metric_long.parquet")
SUMMARY_PATH <- rel("benchmark", "02_metrics", "_m", "synthetic_metric_summary.parquet")
FIG_DIR      <- rel("benchmark", "02_metrics", "figures")
TABLE_DIR    <- rel("benchmark", "02_metrics", "_m")
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------
MAIN_METHOD_ORDER <- c(
  "isograph_baseline",
  "isograph_latent",
  "isograph_graph",
  "isograph_vae",
  "wgcna_gene"
)

COMPUTE_METHOD_ORDER <- c(
  "isograph_vae",
  "isograph_vae_gpu",
  "wgcna_gene"
)

METHOD_ORDER <- unique(c(MAIN_METHOD_ORDER, COMPUTE_METHOD_ORDER))

METHOD_LABELS <- c(
  isograph_baseline = "IsoGraph Baseline",
  isograph_latent   = "IsoGraph Latent",
  isograph_graph    = "IsoGraph Graph",
  isograph_vae      = "IsoGraph VAE",
  isograph_vae_gpu  = "IsoGraph VAE GPU",
  wgcna_gene        = "WGCNA"
)

METHOD_LABELS_SHORT <- c(
  isograph_baseline = "Baseline",
  isograph_latent   = "Latent",
  isograph_graph    = "Graph",
  isograph_vae      = "VAE",
  isograph_vae_gpu  = "VAE GPU",
  wgcna_gene        = "WGCNA"
)

COMPUTE_LABELS <- c(
  isograph_baseline = "CPU",
  isograph_latent   = "CPU",
  isograph_graph    = "CPU",
  isograph_vae      = "CPU",
  isograph_vae_gpu  = "GPU",
  wgcna_gene        = "CPU"
)

# Okabe-Ito colorblind-safe palette, with WGCNA as neutral reference.
METHOD_COLORS <- c(
  isograph_baseline = "#0072B2",
  isograph_latent   = "#E69F00",
  isograph_graph    = "#009E73",
  isograph_vae      = "#D55E00",
  isograph_vae_gpu  = "#CC79A7",
  wgcna_gene        = "#666666"
)

SCENARIO_ORDER <- c(
  "idealized_switching",
  "noise_stress",
  "feature_space_interactions",
  "non_switching_background",
  "unequal_isoform_abundance"
)

SCENARIO_LABELS <- c(
  idealized_switching        = "Idealized switching",
  noise_stress               = "Noise stress",
  feature_space_interactions = "Feature interactions",
  non_switching_background   = "Non-switching background",
  unequal_isoform_abundance  = "Unequal abundance",
  scale                      = "Scale"
)

MAIN_METRICS <- c(
  "metrics_module_recovery",
  "metrics_switch_gene_detection_rate",
  "metrics_nonswitch_gene_module_rate"
)

PWC_LABEL <- "{p.adj.signif}"
PWC_SYMNUM_ARGS <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  symbols = c("****", "***", "**", "*", "ns")
)

# ---------------------------------------------------------------------------
# Shared theme and helpers
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.38, "cm"),
      strip.text         = element_text(size = 8.3, face = "plain"),
      strip.background   = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}

save_fig <- function(p, name, width, height) {
  pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(pdf_path, plot = p, width = width, height = height, units = "in",
         device = cairo_pdf)
  ggsave(png_path, plot = p, width = width, height = height, units = "in",
         dpi = 300)
  cat("  ", name, " saved\n", sep = "")
}

format_param <- function(x) {
  if (is.numeric(x)) {
    out <- sprintf("%.2f", x)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    return(out)
  }
  as.character(x)
}

ci_summary <- function(df) {
  df |>
    summarise(
      n = sum(is.finite(.data$value)),
      mean_val = mean(.data$value, na.rm = TRUE),
      sd_val = sd(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      se_val = ifelse(.data$n > 1 & is.finite(.data$sd_val),
                      .data$sd_val / sqrt(.data$n), 0),
      ci_delta = ifelse(.data$n > 1,
                        qt(0.975, df = pmax(.data$n - 1, 1)) * .data$se_val, 0),
      ci_low = pmax(0, .data$mean_val - .data$ci_delta),
      ci_high = pmin(1, .data$mean_val + .data$ci_delta)
    )
}

iqr_string <- function(values, digits = 1) {
  values <- values[is.finite(values)]
  if (length(values) == 0) return("NA")
  qs <- quantile(values, c(0.25, 0.5, 0.75), na.rm = TRUE)
  fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
  sprintf(fmt, qs[[2]], qs[[1]], qs[[3]])
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
load_raw <- function(path, method_order = METHOD_ORDER) {
  read_parquet(path) |>
    filter(.data$status == "completed", .data$run_method %in% method_order) |>
    mutate(run_method = factor(.data$run_method, levels = method_order))
}

load_long <- function(path, method_order = METHOD_ORDER) {
  read_parquet(path) |>
    filter(.data$method %in% method_order) |>
    mutate(method = factor(.data$method, levels = method_order))
}

load_summary <- function(path, method_order = METHOD_ORDER) {
  read_parquet(path) |>
    filter(.data$method %in% method_order) |>
    mutate(method = factor(.data$method, levels = method_order))
}

# ---------------------------------------------------------------------------
# Main figure: CPU accuracy and specificity benchmark
# ---------------------------------------------------------------------------
scenario_count_labels <- function(long_df, methods = MAIN_METHOD_ORDER,
                                  scenarios = SCENARIO_ORDER) {
  counts <- long_df |>
    filter(
      .data$scenario %in% scenarios,
      .data$method %in% methods,
      .data$metric == "metrics_module_recovery",
      is.finite(.data$value)
    ) |>
    group_by(.data$scenario, .data$method) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(.data$scenario) |>
    summarise(n_min = min(.data$n), n_max = max(.data$n), .groups = "drop")

  labels <- vapply(scenarios, function(scen) {
    row <- counts[counts$scenario == scen, ]
    n_text <- if (nrow(row) == 0) {
      "n=0 per method"
    } else if (row$n_min == row$n_max) {
      paste0("n=", row$n_min, " per method")
    } else {
      paste0("n=", row$n_min, "-", row$n_max, " per method")
    }
    paste0(SCENARIO_LABELS[[scen]], "\n", n_text)
  }, character(1))
  setNames(labels, scenarios)
}

metric_box_data <- function(long_df, metric_name, methods = MAIN_METHOD_ORDER,
                            scenarios = SCENARIO_ORDER) {
  scenario_labels <- scenario_count_labels(long_df, methods, scenarios)

  long_df |>
    filter(
      .data$scenario %in% scenarios,
      .data$method %in% methods,
      .data$metric == metric_name,
      is.finite(.data$value)
    ) |>
    mutate(
      value_plot = ifelse(.data$metric == "metrics_nonswitch_gene_module_rate",
                          1 - .data$value, .data$value),
      value_plot = pmin(1, pmax(0, .data$value_plot)),
      method = factor(.data$method, levels = methods),
      scenario_label = factor(
        scenario_labels[as.character(.data$scenario)],
        levels = scenario_labels[scenarios]
      )
    )
}

box_metric_panel <- function(long_df, metric_name, y_label, tag,
                             show_x_labels = FALSE) {
  sub <- metric_box_data(long_df, metric_name)
  stat_sub <- sub |>
    filter(.data$method %in% c("isograph_vae", "wgcna_gene"))

  p <- ggpubr::ggboxplot(
    sub,
    x = "method",
    y = "value_plot",
    fill = "method",
    palette = METHOD_COLORS[MAIN_METHOD_ORDER],
    width = 0.62,
    outliers = TRUE,
    outlier.size = 0.45,
    outlier.alpha = 0.25,
    color = "grey25",
    xlab = FALSE,
    ylab = y_label,
    ggtheme = theme_pub()
  ) +
    facet_wrap(~ scenario_label, nrow = 1) +
    scale_x_discrete(labels = METHOD_LABELS_SHORT[MAIN_METHOD_ORDER]) +
    scale_fill_manual(
      values = METHOD_COLORS[MAIN_METHOD_ORDER],
      labels = METHOD_LABELS[MAIN_METHOD_ORDER],
      breaks = MAIN_METHOD_ORDER
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    ggpubr::geom_pwc(
      data = stat_sub,
      mapping = aes(x = .data$method, y = .data$value_plot,
                    group = .data$method),
      method = "wilcox_test",
      label = PWC_LABEL,
      p.adjust.method = "BH",
      p.adjust.by = "panel",
      symnum.args = PWC_SYMNUM_ARGS,
      hide.ns = TRUE,
      y.position = 1.01,
      tip.length = 0.01,
      bracket.nudge.y = 0,
      bracket.shorten = 0.04,
      size = 0.25,
      label.size = 2.6,
      inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(0, 1.08), clip = "off") +
    labs(x = NULL, y = y_label, fill = NULL, tag = tag) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
      plot.tag = element_text(size = 10, face = "bold")
    )

  if (!show_x_labels) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  p
}

make_fig1 <- function(long_df) {
  pA <- box_metric_panel(
    long_df,
    "metrics_module_recovery",
    "Module recovery (AUC; 1 = perfect)",
    "A"
  )
  pB <- box_metric_panel(
    long_df,
    "metrics_switch_gene_detection_rate",
    "Switching gene detection (recall; 1 = all recovered)",
    "B"
  )
  pC <- box_metric_panel(
    long_df,
    "metrics_nonswitch_gene_module_rate",
    "Non-switching gene exclusion (1 - false-positive rate)",
    "C",
    show_x_labels = TRUE
  )

  (pA / pB / pC) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# ---------------------------------------------------------------------------
# Parameter-response supplements
# ---------------------------------------------------------------------------
response_dot_fig <- function(long_df, scenario, x_col, facet_col,
                             x_label, facet_label,
                             methods = MAIN_METHOD_ORDER) {
  sub <- long_df |>
    filter(
      .data$scenario == !!scenario,
      .data$method %in% methods,
      .data$metric == "metrics_module_recovery",
      !is.na(.data[[x_col]]),
      !is.na(.data[[facet_col]]),
      is.finite(.data$value)
    )

  if (nrow(sub) == 0) return(NULL)

  plot_df <- sub |>
    mutate(
      x_value = .data[[x_col]],
      facet_value = .data[[facet_col]]
    ) |>
    group_by(.data$method, .data$x_value, .data$facet_value) |>
    ci_summary() |>
    mutate(method = factor(.data$method, levels = methods))

  x_levels <- sort(unique(plot_df$x_value))
  facet_levels <- sort(unique(plot_df$facet_value))

  facet_counts <- plot_df |>
    group_by(.data$facet_value) |>
    summarise(n_min = min(.data$n), n_max = max(.data$n), .groups = "drop")

  facet_labels <- vapply(facet_levels, function(v) {
    row <- facet_counts[facet_counts$facet_value == v, ]
    n_text <- if (row$n_min == row$n_max) {
      paste0("n=", row$n_min, " per cell")
    } else {
      paste0("n=", row$n_min, "-", row$n_max, " per cell")
    }
    paste0(facet_label, " = ", format_param(v), "\n", n_text)
  }, character(1))
  names(facet_labels) <- as.character(facet_levels)

  plot_df <- plot_df |>
    mutate(
      x_plot = factor(.data$x_value, levels = x_levels,
                      labels = format_param(x_levels)),
      facet_label = factor(
        facet_labels[as.character(.data$facet_value)],
        levels = facet_labels[as.character(facet_levels)]
      )
    )

  stat_df <- sub |>
    filter(.data$method %in% c("isograph_vae", "wgcna_gene")) |>
    mutate(
      x_value = .data[[x_col]],
      facet_value = .data[[facet_col]],
      x_plot = factor(.data$x_value, levels = x_levels,
                      labels = format_param(x_levels)),
      facet_label = factor(
        facet_labels[as.character(.data$facet_value)],
        levels = facet_labels[as.character(facet_levels)]
      ),
      method = factor(.data$method, levels = methods)
    )

  ggplot(plot_df, aes(x = .data$x_plot, y = .data$mean_val,
                      color = .data$method, group = .data$method)) +
    geom_errorbar(
      aes(ymin = .data$ci_low, ymax = .data$ci_high),
      width = 0.18,
      linewidth = 0.45,
      position = position_dodge(width = 0.52)
    ) +
    geom_point(size = 2.1, position = position_dodge(width = 0.52)) +
    facet_wrap(~ facet_label, nrow = 1) +
    scale_color_manual(
      values = METHOD_COLORS[methods],
      labels = METHOD_LABELS[methods],
      breaks = methods
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    ggpubr::geom_pwc(
      data = stat_df,
      mapping = aes(x = .data$x_plot, y = .data$value,
                    group = .data$method),
      method = "wilcox_test",
      group.by = "x.var",
      label = PWC_LABEL,
      p.adjust.method = "BH",
      p.adjust.by = "panel",
      symnum.args = PWC_SYMNUM_ARGS,
      hide.ns = TRUE,
      y.position = 1.01,
      dodge = 0.52,
      tip.length = 0.01,
      bracket.nudge.y = 0,
      bracket.shorten = 0.04,
      size = 0.25,
      label.size = 2.4,
      inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(0, 1.08), clip = "off") +
    labs(
      x = x_label,
      y = "Module recovery (AUC; 1 = perfect)",
      color = NULL
    ) +
    theme_pub() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

parameter_line_fig <- function(long_df, scenario, param_col, param_label,
                               metric_names, metric_labels,
                               transform_nonswitch = TRUE,
                               methods = MAIN_METHOD_ORDER) {
  sub <- long_df |>
    filter(
      .data$scenario == !!scenario,
      .data$method %in% methods,
      .data$metric %in% metric_names,
      !is.na(.data[[param_col]]),
      is.finite(.data$value)
    ) |>
    mutate(
      value = ifelse(
        transform_nonswitch &
          .data$metric == "metrics_nonswitch_gene_module_rate",
        1 - .data$value,
        .data$value
      ),
      value = pmin(1, pmax(0, .data$value)),
      param_value = .data[[param_col]]
    )

  if (nrow(sub) == 0) return(NULL)

  plot_df <- sub |>
    group_by(.data$method, .data$param_value, .data$metric) |>
    ci_summary() |>
    mutate(
      method = factor(.data$method, levels = methods),
      metric_label = factor(metric_labels[as.character(.data$metric)],
                            levels = metric_labels[metric_names])
    )

  ggplot(plot_df, aes(x = .data$param_value, y = .data$mean_val,
                      color = .data$method, group = .data$method)) +
    geom_errorbar(aes(ymin = .data$ci_low, ymax = .data$ci_high),
                  width = 0.02, linewidth = 0.45) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ metric_label, nrow = 1) +
    scale_color_manual(
      values = METHOD_COLORS[methods],
      labels = METHOD_LABELS[methods],
      breaks = methods
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(x = param_label, y = NULL, color = NULL) +
    theme_pub() +
    theme(legend.position = "right")
}

# ---------------------------------------------------------------------------
# Scale compute supplement
# ---------------------------------------------------------------------------
compute_box_panel <- function(plot_df, y_label, tag,
                              log_y = FALSE, show_x_labels = TRUE) {
  sub <- plot_df |>
    filter(.data$method %in% COMPUTE_METHOD_ORDER,
           !is.na(.data$n_genes),
           is.finite(.data$value)) |>
    mutate(
      method = factor(.data$method, levels = COMPUTE_METHOD_ORDER),
      n_genes = factor(.data$n_genes, levels = sort(unique(.data$n_genes)))
    )

  if (nrow(sub) == 0) return(NULL)

  p <- ggpubr::ggboxplot(
    sub,
    x = "n_genes",
    y = "value",
    fill = "method",
    palette = METHOD_COLORS[COMPUTE_METHOD_ORDER],
    width = 0.68,
    outliers = TRUE,
    outlier.size = 0.45,
    outlier.alpha = 0.25,
    color = "grey25",
    xlab = "Number of genes",
    ylab = y_label,
    ggtheme = theme_pub()
  ) +
    scale_fill_manual(
      values = METHOD_COLORS[COMPUTE_METHOD_ORDER],
      labels = METHOD_LABELS[COMPUTE_METHOD_ORDER],
      breaks = COMPUTE_METHOD_ORDER
    ) +
    ggpubr::geom_pwc(
      mapping = aes(x = .data$n_genes, y = .data$value,
                    group = .data$method),
      method = "wilcox_test",
      group.by = "x.var",
      ref.group = "isograph_vae",
      label = PWC_LABEL,
      p.adjust.method = "BH",
      p.adjust.by = "panel",
      symnum.args = PWC_SYMNUM_ARGS,
      hide.ns = TRUE,
      dodge = 0.78,
      tip.length = 0.01,
      bracket.nudge.y = 0.04,
      bracket.shorten = 0.03,
      size = 0.25,
      label.size = 2.4
    ) +
    labs(x = "Number of genes", y = y_label, fill = NULL, tag = tag) +
    theme(
      legend.position = "right",
      plot.tag = element_text(size = 10, face = "bold")
    )

  if (log_y) {
    p <- p + scale_y_log10(labels = label_comma())
  } else {
    p <- p + scale_y_continuous(labels = label_comma(),
                                expand = expansion(mult = c(0.02, 0.06)))
  }
  if (!show_x_labels) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  p
}

make_scale_compute_fig <- function(raw_df, long_df) {
  mem_col <- if ("measurement_max_rss_mb" %in% names(raw_df)) {
    "measurement_max_rss_mb"
  } else {
    "max_rss_mb"
  }

  runtime_df <- long_df |>
    filter(
      .data$scenario == "scale",
      .data$method %in% COMPUTE_METHOD_ORDER,
      .data$metric == "measurement_elapsed_sec",
      !is.na(.data$run_n_genes)
    ) |>
    transmute(
      method = as.character(.data$method),
      n_genes = .data$run_n_genes,
      value = .data$value
    )

  memory_df <- raw_df |>
    filter(.data$run_scenario == "scale",
           .data$run_method %in% COMPUTE_METHOD_ORDER,
           !is.na(.data$run_n_genes)) |>
    transmute(
      method = as.character(.data$run_method),
      n_genes = .data$run_n_genes,
      value = .data[[mem_col]] / 1024
    )

  pA <- compute_box_panel(
    runtime_df,
    "Runtime (s, log10 scale)",
    "A",
    log_y = TRUE,
    show_x_labels = FALSE
  )
  pB <- compute_box_panel(
    memory_df,
    "Peak host RAM (GB)",
    "B"
  )

  if (is.null(pA) || is.null(pB)) return(NULL)

  (pA / pB) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
make_table1 <- function(summary) {
  scenarios <- SCENARIO_ORDER[SCENARIO_ORDER %in% unique(summary$scenario)]
  methods   <- MAIN_METHOD_ORDER[MAIN_METHOD_ORDER %in% unique(summary$method)]

  metric_map <- c(
    metrics_module_recovery            = "Module recovery",
    metrics_switch_gene_detection_rate = "Switch detection",
    metrics_nonswitch_gene_module_rate =
      "Non-switching gene module rate (lower is better)"
  )

  rows <- list()
  for (scen in scenarios) {
    for (meth in methods) {
      row <- list(
        Scenario = SCENARIO_LABELS[[scen]],
        Method   = METHOD_LABELS[[meth]],
        Compute  = COMPUTE_LABELS[[meth]]
      )
      for (metric in names(metric_map)) {
        sub <- summary |>
          filter(.data$scenario == scen,
                 .data$method == meth,
                 .data$metric == !!metric)
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
        filter(.data$scenario == scen,
               .data$method == meth,
               .data$metric == "measurement_elapsed_sec")
      row[["Runtime, s median"]] <- if (nrow(rt) > 0) {
        sprintf("%.1f", rt$median)
      } else {
        "NA"
      }
      rows[[length(rows) + 1]] <- as.data.frame(row, check.names = FALSE)
    }
  }

  tbl <- do.call(rbind, rows)
  out_path <- file.path(TABLE_DIR, "table1_benchmark_summary.csv")
  write.csv(tbl, out_path, row.names = FALSE)
  cat("  table1_benchmark_summary.csv:", nrow(tbl), "rows\n")
}

make_compute_table <- function(raw_df, long_df) {
  mem_col <- if ("measurement_max_rss_mb" %in% names(raw_df)) {
    "measurement_max_rss_mb"
  } else {
    "max_rss_mb"
  }
  if (!"torch_gpu_peak_reserved_mb" %in% names(raw_df)) {
    raw_df$torch_gpu_peak_reserved_mb <- NA_real_
  }

  runtime_long <- long_df |>
    filter(.data$scenario == "scale",
           .data$method %in% COMPUTE_METHOD_ORDER,
           .data$metric == "measurement_elapsed_sec") |>
    transmute(
      method = as.character(.data$method),
      n_genes = .data$run_n_genes,
      runtime_sec = .data$value
    )

  resource_raw <- raw_df |>
    filter(.data$run_scenario == "scale",
           .data$run_method %in% COMPUTE_METHOD_ORDER) |>
    transmute(
      method = as.character(.data$run_method),
      n_genes = .data$run_n_genes,
      host_ram_peak_gb = .data[[mem_col]] / 1024,
      gpu_peak_reserved_gb = .data$torch_gpu_peak_reserved_mb / 1024
    )

  rows <- list()
  for (genes in sort(unique(resource_raw$n_genes))) {
    for (meth in COMPUTE_METHOD_ORDER) {
      rt <- runtime_long |>
        filter(.data$n_genes == genes, .data$method == meth)
      mem <- resource_raw |>
        filter(.data$n_genes == genes, .data$method == meth)
      if (nrow(rt) == 0 && nrow(mem) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        `Number of genes` = genes,
        Method = METHOD_LABELS[[meth]],
        Compute = COMPUTE_LABELS[[meth]],
        N = max(nrow(rt), nrow(mem)),
        `Runtime, s median [IQR]` =
          iqr_string(rt$runtime_sec, digits = 1),
        `Peak host RAM, GB median [IQR]` =
          iqr_string(mem$host_ram_peak_gb, digits = 1),
        `Peak GPU VRAM reserved, GB median [IQR]` =
          iqr_string(mem$gpu_peak_reserved_gb, digits = 1),
        check.names = FALSE
      )
    }
  }

  tbl <- do.call(rbind, rows)
  out_path <- file.path(TABLE_DIR, "tableS_scale_compute_summary.csv")
  write.csv(tbl, out_path, row.names = FALSE)
  cat("  tableS_scale_compute_summary.csv:", nrow(tbl), "rows\n")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
cat("Loading data...\n")
if (!file.exists(RAW_PATH))     stop("Missing: ", RAW_PATH)
if (!file.exists(LONG_PATH))    stop("Missing: ", LONG_PATH)
if (!file.exists(SUMMARY_PATH)) stop("Missing: ", SUMMARY_PATH)

raw     <- load_raw(RAW_PATH)
long    <- load_long(LONG_PATH)
summary <- load_summary(SUMMARY_PATH)
cat("  Completed runs:", nrow(raw),
    " | Long metric rows:", nrow(long),
    " | Summary rows:", nrow(summary), "\n")

cat("Generating figures...\n")

tryCatch({
  fig1 <- make_fig1(long)
  save_fig(fig1, "fig1_benchmark_overview", width = 12, height = 8.2)
}, error = function(e) warning("fig1 error: ", conditionMessage(e)))

tryCatch({
  p <- response_dot_fig(
    long, "idealized_switching",
    "run_switching_fraction", "run_noise_sd",
    "Switching gene fraction", "Noise SD"
  )
  if (!is.null(p)) save_fig(p, "figS1_idealized_switching", width = 10.5, height = 4.6)
  else cat("  figS1: no data, skipping\n")
}, error = function(e) warning("figS1 error: ", conditionMessage(e)))

tryCatch({
  p <- response_dot_fig(
    long, "noise_stress",
    "run_count_dispersion", "run_noise_sd",
    "Count dispersion", "Noise SD"
  )
  if (!is.null(p)) save_fig(p, "figS2_noise_stress", width = 12, height = 4.6)
  else cat("  figS2: no data, skipping\n")
}, error = function(e) warning("figS2 error: ", conditionMessage(e)))

tryCatch({
  p <- response_dot_fig(
    long, "feature_space_interactions",
    "run_interaction_strength", "run_interaction_fraction",
    "Interaction strength", "Interaction fraction"
  )
  if (!is.null(p)) save_fig(p, "figS3_feature_interactions", width = 10.5, height = 4.6)
  else cat("  figS3: no data, skipping\n")
}, error = function(e) warning("figS3 error: ", conditionMessage(e)))

tryCatch({
  p <- parameter_line_fig(
    long, "non_switching_background", "run_switching_fraction",
    "Switching gene fraction",
    metric_names = c("metrics_module_recovery",
                     "metrics_nonswitch_gene_module_rate"),
    metric_labels = c(
      metrics_module_recovery = "Module recovery (AUC; 1 = perfect)",
      metrics_nonswitch_gene_module_rate =
        "Non-switching gene exclusion (1 - false-positive rate)"
    )
  )
  if (!is.null(p)) save_fig(p, "figS4_nonswitching_background", width = 9, height = 3.9)
  else cat("  figS4: no data, skipping\n")
}, error = function(e) warning("figS4 error: ", conditionMessage(e)))

tryCatch({
  p <- parameter_line_fig(
    long, "unequal_isoform_abundance", "run_abundance_imbalance",
    "Isoform abundance imbalance ratio",
    metric_names = c("metrics_module_recovery",
                     "metrics_switch_gene_detection_rate"),
    metric_labels = c(
      metrics_module_recovery = "Module recovery (AUC; 1 = perfect)",
      metrics_switch_gene_detection_rate =
        "Switching gene detection (recall; 1 = all recovered)"
    )
  )
  if (!is.null(p)) save_fig(p, "figS5_unequal_abundance", width = 9, height = 3.9)
  else cat("  figS5: no data, skipping\n")
}, error = function(e) warning("figS5 error: ", conditionMessage(e)))

tryCatch({
  p <- make_scale_compute_fig(raw, long)
  if (!is.null(p)) save_fig(p, "figS6_scale_compute", width = 8.6, height = 6.2)
  else cat("  figS6: no data, skipping\n")
}, error = function(e) warning("figS6 error: ", conditionMessage(e)))

tryCatch(make_table1(summary), error = function(e) warning("table1 error: ", conditionMessage(e)))
tryCatch(make_compute_table(raw, long), error = function(e) warning("compute table error: ", conditionMessage(e)))

cat("Done. Outputs in", FIG_DIR, "\n")
