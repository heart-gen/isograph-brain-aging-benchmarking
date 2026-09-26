# Supplementary real-data figure: how much of the switch-unique picture is the FDR line?
#
# A gene is switch-unique when its switch test clears FDR <= alpha conditional on abundance
# AND its abundance test conditional on the switch coordinate does not. The second half
# accepts a null at the same alpha, so the class is sensitive to where the line sits. The
# subsection's result is a regional pattern -- some analyses keep their switch-unique genes
# under composition adjustment, the two GTEx cortical ones do not -- and this figure asks
# whether that pattern is a property of the data or of alpha = 0.10.
#
# It is built so it can fail: panel C labels analyses whose call flips between alphas as
# "threshold-dependent" rather than rounding them to the convenient side, and several do.
#
# (A) switch-unique counts against alpha, one line per analysis. Counts must move with
#     alpha; the question is whether the ordering does.
# (B) gene-level persistence (n_overlap / switch-unique before adjustment) against alpha,
#     for the 12 analyses with a composition-adjusted arm.
# (C) the range of persistence across alphas per analysis, coloured by whether the analysis
#     collapses, retains, or changes its call. Open symbols mark analyses whose unadjusted
#     set falls below 10 genes at some alpha, where the ratio is single digits over single
#     digits and the call should not be leaned on.
#
# Reads 03_module_characterization/_m/switch_unique_threshold{,_calls}.csv (CSV, not
# parquet, so the figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figSwitchUniqueThreshold.{pdf,png}.
# Run: bash 03_module_characterization/_h/04c.switch_unique_threshold.sh
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
COMP    <- rel("03_module_characterization", "_m")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, matched to the other real-data supplements.
CLASS_COLORS <- c(`Limbic / striatal` = "#0072B2",
                  `Cortical`          = "#D55E00",
                  `Disease (SCZD)`    = "#CC79A7")
CALL_COLORS <- c(`collapses at every alpha`  = "#D55E00",
                 `retains at every alpha`    = "#0072B2",
                 `threshold-dependent`       = "#999999")
REFERENCE_ALPHA <- 0.10

CLASS_OF <- c(
  "SCZD (caudate)"                       = "Disease (SCZD)",
  "aging caudate"                        = "Limbic / striatal",
  "aging hippocampus"                    = "Limbic / striatal",
  "aging DLPFC"                          = "Cortical",
  "GTEx amygdala"                        = "Limbic / striatal",
  "GTEx hippocampus"                     = "Limbic / striatal",
  "GTEx caudate_basal_ganglia"           = "Limbic / striatal",
  "GTEx putamen_basal_ganglia"           = "Limbic / striatal",
  "GTEx nucleus_accumbens_basal_ganglia" = "Limbic / striatal",
  "GTEx anterior_cingulate_cortex_ba24"  = "Cortical",
  "GTEx frontal_cortex_ba9"              = "Cortical",
  "GTEx cortex"                          = "Cortical",
  "GTEx cerebellum"                      = "Cerebellar / other",
  "GTEx cerebellar_hemisphere"           = "Cerebellar / other",
  "GTEx hypothalamus"                    = "Cerebellar / other",
  "GTEx spinal_cord_cervical_c_1"        = "Cerebellar / other",
  "GTEx substantia_nigra"                = "Cerebellar / other")
CLASS_COLORS <- c(CLASS_COLORS, `Cerebellar / other` = "#009E73")

PRETTY <- c(
  "SCZD (caudate)"                       = "Caudate, SCZD (BrainSEQ)",
  "aging caudate"                        = "Caudate (BrainSEQ)",
  "aging hippocampus"                    = "Hippocampus (BrainSEQ)",
  "aging DLPFC"                          = "DLPFC (BrainSEQ)",
  "GTEx amygdala"                        = "Amygdala (GTEx)",
  "GTEx hippocampus"                     = "Hippocampus (GTEx)",
  "GTEx caudate_basal_ganglia"           = "Caudate (GTEx)",
  "GTEx putamen_basal_ganglia"           = "Putamen (GTEx)",
  "GTEx nucleus_accumbens_basal_ganglia" = "N. accumbens (GTEx)",
  "GTEx anterior_cingulate_cortex_ba24"  = "ACC BA24 (GTEx)",
  "GTEx frontal_cortex_ba9"              = "Frontal ctx BA9 (GTEx)",
  "GTEx cortex"                          = "Cortex (GTEx)",
  "GTEx cerebellum"                      = "Cerebellum (GTEx)",
  "GTEx cerebellar_hemisphere"           = "Cerebellar hem. (GTEx)",
  "GTEx hypothalamus"                    = "Hypothalamus (GTEx)",
  "GTEx spinal_cord_cervical_c_1"        = "Spinal cord C1 (GTEx)",
  "GTEx substantia_nigra"                = "Substantia nigra (GTEx)")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.2),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.34, "cm"),
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

thr   <- read.csv(file.path(COMP, "switch_unique_threshold.csv"))
calls <- read.csv(file.path(COMP, "switch_unique_threshold_calls.csv"))

thr <- thr |>
  mutate(class = factor(unname(CLASS_OF[label]), names(CLASS_COLORS)),
         pretty = unname(PRETTY[label]))

# ---------------------------------------------------------------------------
# Panel A - switch-unique counts against alpha
# ---------------------------------------------------------------------------
pA <- ggplot(thr, aes(alpha, switch_unique_base, group = pretty, colour = class)) +
  geom_vline(xintercept = REFERENCE_ALPHA, linewidth = 0.3, linetype = "dashed",
             colour = "grey55") +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = CLASS_COLORS, name = NULL) +
  # pseudo-log: the counts span 1 to >1,600 and the small analyses matter to the argument.
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 1, 3, 10, 30, 100, 300, 1000)) +
  scale_x_continuous(breaks = sort(unique(thr$alpha))) +
  labs(x = "FDR threshold", y = "Switch-unique genes") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_pub() +
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Panel B - gene-level persistence against alpha
# ---------------------------------------------------------------------------
datB <- thr |> filter(!is.na(persistence))
pB <- ggplot(datB, aes(alpha, persistence, group = pretty, colour = class)) +
  geom_vline(xintercept = REFERENCE_ALPHA, linewidth = 0.3, linetype = "dashed",
             colour = "grey55") +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_continuous(breaks = sort(unique(datB$alpha))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "FDR threshold",
       y = "Unadjusted switch-unique genes still\nswitch-unique after adjustment") +
  theme_pub()

# ---------------------------------------------------------------------------
# Panel C - the call, and whether it survives every alpha
# ---------------------------------------------------------------------------
# The "(small set)" suffix is dropped from the colour key and carried by the symbol fill
# instead, so the three calls stay legible while the caveat stays attached to its rows.
datC <- calls |>
  mutate(small = grepl("small set", call),
         call = sub(" \\(small set\\)$", "", call),
         call = factor(call, names(CALL_COLORS)),
         pretty = unname(PRETTY[label])) |>
  arrange(min_persistence) |>
  mutate(pretty = factor(pretty, levels = pretty))

pC <- ggplot(datC, aes(y = pretty, colour = call)) +
  geom_segment(aes(x = min_persistence, xend = max_persistence, yend = pretty),
               linewidth = 0.6) +
  geom_point(aes(x = min_persistence, fill = call, shape = small), size = 1.9,
             stroke = 0.6) +
  geom_point(aes(x = max_persistence, fill = call, shape = small), size = 1.9,
             stroke = 0.6) +
  scale_colour_manual(values = CALL_COLORS, name = NULL, drop = FALSE) +
  scale_fill_manual(values = CALL_COLORS, guide = "none", drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Persistence across FDR thresholds (range)", y = NULL,
       caption = "Open symbols: fewer than 10 unadjusted switch-unique genes at some threshold") +
  guides(colour = guide_legend(nrow = 3, override.aes = list(shape = NA))) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
        plot.caption = element_text(size = 6.3, colour = "grey40", hjust = 0),
        legend.position = "bottom")

fig <- ((pA | pB) + plot_layout(widths = c(1, 1))) / pC +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figSwitchUniqueThreshold", width = 7.4, height = 7.4)
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
