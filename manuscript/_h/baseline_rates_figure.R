# Supplementary real-data figure: honest three-baseline per-module rate comparison.
# IsoGraph is NOT globally superior on module-level metrics; phenotype signal lives in the
# switch features (both switch-fed methods beat both abundance-fed ones), GO-enrichment is
# low for IsoGraph by construction. Reads
# 04_module_characterization/_m/baseline_comparison/baseline_comparison_pooled.parquet, writes
# figBaselineRates.{pdf,png} to manuscript/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/baseline_rates_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
BC_FILE <- rel("04_module_characterization", "_m", "baseline_comparison", "baseline_comparison_pooled.parquet")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito; fill encodes the FEATURE class so the switch-vs-abundance story reads off colour.
FEATURE_COLORS <- c(`switch+abundance` = "#D55E00", `switch-only` = "#E69F00",
                    abundance = "#999999")
METHOD_LABELS  <- c(isograph = "IsoGraph", wgcna_switch_only = "WGCNA\nswitch-only",
                    wgcna_multiplex = "WGCNA\nmultiplex", wgcna_gene = "WGCNA\nabundance")
METRIC_LABELS  <- c(mean_frac_pheno_sig = "Phenotype-sig rate",
                    mean_frac_both = "Phenotype-sig AND GO rate",
                    mean_frac_go_enriched = "GO-enriched rate")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.position    = "bottom",
      strip.text         = element_text(size = 8.3),
      strip.background   = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}

save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300)
  cat("  ", name, " saved\n", sep = "")
}

bc <- as.data.frame(read_parquet(BC_FILE))

dd <- bc |>
  select(method, features, names(METRIC_LABELS)) |>
  pivot_longer(names(METRIC_LABELS), names_to = "metric", values_to = "rate") |>
  mutate(metric = factor(METRIC_LABELS[metric], levels = unname(METRIC_LABELS)),
         method = factor(method, names(METHOD_LABELS)),
         features = factor(features, names(FEATURE_COLORS)))

fig <- ggplot(dd, aes(method, rate, fill = features)) +
  geom_col(width = 0.74, alpha = 0.92) +
  geom_text(aes(label = sprintf("%.2f", rate)), vjust = -0.4, size = 2.4) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_manual(values = FEATURE_COLORS) +
  scale_x_discrete(labels = METHOD_LABELS) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Mean per-module rate (pooled, 17 regions)") +
  theme_pub() + theme(axis.text.x = element_text(size = 7))

save_fig(fig, "figBaselineRates", width = 7.2, height = 3.2)
cat("Done. Output in", FIG_DIR, "\n")
