# Main real-data figure: how much of the DTU-without-DGE layer survives cell-type
# composition adjustment.
#
# "Is this just shifting cell-type proportions?" is the first objection put to any bulk
# brain DTU result, and the repository answers it -- but until now only in prose
# (03_module_characterization/_m/COMPOSITION_ADJUSTMENT_SUMMARY.md). The answer also moves
# a headline number: the SCZD composition-unique count falls 34 -> 2 under adjustment,
# while the aging layer largely survives and GTEx replicates it only partially.
#
# The figure is deliberately built so the losses are as visible as the survivals. A test
# that always survives is not evidence; this one plainly does not always survive.
#
# (A) base -> composition-adjusted composition-unique gene counts, paired within
#     cohort x region (log1p count axis: the counts span 0-545);
# (B) retained fraction per analysis, with the five GTEx regions that have no defensibly
#     matched snRNA reference shown as EXCLUDED rather than silently absent;
# (C) the marker cut -- module member genes are not over-represented for cell-type
#     markers, so the modules are not simply bags of marker genes.
#
# Reads 03_module_characterization/_m/composition_adjustment.parquet,
#       02_module_discovery/gtex/_m/composition/composition_adjustment_gtex.parquet,
#       02_module_discovery/*/*/_m/isograph_vae/celltype_composition/marker_enrichment.parquet
# Writes manuscript/_m/figures/figCompositionRobustness.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/composition_robustness_figure.R
suppressPackageStartupMessages({
  library(arrow)
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
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, shared with the other real-data figures. The split that matters here is
# anatomical, not cohort: the two GTEx regions that collapse entirely are both cortical.
CLASS_COLORS <- c(`Limbic / striatal` = "#0072B2",
                  `Cortical`          = "#D55E00",
                  `Disease (SCZD)`    = "#CC79A7")
COHORT_SHAPE <- c(BrainSEQ = 16, GTEx = 17)

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.36, "cm"),
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

# ---------------------------------------------------------------------------
# Load: BrainSEQ (disease + aging) and the GTEx aging replication arm
# ---------------------------------------------------------------------------
bs <- as.data.frame(read_parquet(
  rel("03_module_characterization", "_m", "composition_adjustment.parquet"))) |>
  mutate(cohort = "BrainSEQ")
gt <- as.data.frame(read_parquet(
  rel("02_module_discovery", "gtex", "_m", "composition",
      "composition_adjustment_gtex.parquet"))) |>
  mutate(cohort = "GTEx")

# Human labels; anatomical class drives the colour because that is where the pattern is.
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
  "GTEx cortex"                          = "Cortical")

PRETTY <- c(
  "SCZD (caudate)"                       = "Caudate, SCZD",
  "aging caudate"                        = "Caudate",
  "aging hippocampus"                    = "Hippocampus",
  "aging DLPFC"                          = "DLPFC",
  "GTEx amygdala"                        = "Amygdala",
  "GTEx hippocampus"                     = "Hippocampus",
  "GTEx caudate_basal_ganglia"           = "Caudate",
  "GTEx putamen_basal_ganglia"           = "Putamen",
  "GTEx nucleus_accumbens_basal_ganglia" = "N. accumbens",
  "GTEx anterior_cingulate_cortex_ba24"  = "ACC BA24",
  "GTEx frontal_cortex_ba9"              = "Frontal ctx BA9",
  "GTEx cortex"                          = "Cortex")

dat <- bind_rows(bs, gt) |>
  mutate(class  = factor(unname(CLASS_OF[region]), names(CLASS_COLORS)),
         label  = paste0(unname(PRETTY[region]), " (", cohort, ")"),
         cohort = factor(cohort, names(COHORT_SHAPE)))

# Hippocampus in BrainSEQ has zero composition-unique genes before adjustment, so there
# is nothing for adjustment to remove. Keep it in panel A (an honest zero) but drop it
# from the retained-fraction panel, where 0/0 is undefined rather than a loss.
stopifnot(all(!is.na(dat$class)))

# ---------------------------------------------------------------------------
# Panel A - paired base -> adjusted counts
# ---------------------------------------------------------------------------
pa <- dat |>
  select(label, class, cohort, base = comp_unique_base, adj = comp_unique_adj) |>
  pivot_longer(c(base, adj), names_to = "stage", values_to = "n") |>
  mutate(stage = factor(stage, c("base", "adj"),
                        labels = c("Unadjusted", "Composition-\nadjusted")))

pA <- ggplot(pa, aes(stage, n, group = label, colour = class)) +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  geom_point(aes(shape = cohort), size = 1.5) +
  scale_colour_manual(values = CLASS_COLORS, name = NULL) +
  scale_shape_manual(values = COHORT_SHAPE, name = NULL) +
  # log1p keeps the two collapses to zero visible alongside the 545-gene ACC arm.
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  labs(x = NULL, y = "Composition-unique genes") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 2, nrow = 2)) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 6, 0, 0, "pt"),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey92"))

# ---------------------------------------------------------------------------
# Panel B - retained fraction, with the non-deconvolved regions named as excluded
# ---------------------------------------------------------------------------
# The five GTEx regions without a defensibly matched Tran/LIBD snRNA reference
# (cerebellum, cerebellar hemisphere, hypothalamus, spinal cord, substantia nigra) were
# NOT deconvolved rather than forced against a mismatched panel. They are drawn as
# explicit "not deconvolved" rows: a reader must be able to see that this analysis
# covers 12 of 17 analyses, not all of them.
NOT_DECONVOLVED <- c("Cerebellum", "Cerebellar hem.", "Hypothalamus",
                     "Spinal cord C1", "Substantia nigra")

pb <- dat |>
  filter(is.finite(retained_frac)) |>
  select(label, class, retained_frac) |>
  arrange(retained_frac)
pb_excl <- data.frame(label = paste0(NOT_DECONVOLVED, " (GTEx)"),
                      class = NA_character_, retained_frac = NA_real_)
pb_all <- bind_rows(pb, pb_excl)
pb_all$label <- factor(pb_all$label, levels = pb_all$label)

pB <- ggplot(pb_all, aes(retained_frac, label)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = 0, xend = retained_frac, yend = label, colour = class),
               linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(colour = class), size = 1.8, na.rm = TRUE) +
  geom_text(data = subset(pb_all, is.na(retained_frac)),
            aes(x = 0.06, label = "not deconvolved"),
            hjust = 0, size = 2.2, colour = "grey45") +
  scale_colour_manual(values = CLASS_COLORS, na.value = "grey70", guide = "none") +
  scale_x_continuous(limits = c(0, 4.15), breaks = c(0, 1, 2, 3, 4)) +
  labs(x = "Retained fraction after adjustment", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel C - marker cut: modules are not bags of cell-type marker genes
# ---------------------------------------------------------------------------
# Rate, not count: the number of module x cell-type tests differs ~10x between regions,
# so raw counts would compare granularity rather than marker content.
mk_files <- Sys.glob(rel("02_module_discovery", "*", "*", "_m", "isograph_vae",
                         "celltype_composition", "marker_enrichment.parquet"))
mk <- bind_rows(lapply(mk_files, function(f) {
  d <- as.data.frame(read_parquet(f))
  parts <- strsplit(f, .Platform$file.sep)[[1]]
  n <- length(parts)
  data.frame(
    cohort   = ifelse(parts[n - 5] == "gtex", "GTEx", "BrainSEQ"),
    n_tests  = nrow(d),
    enriched = sum(d$fdr_enrich < 0.05, na.rm = TRUE),
    depleted = sum(d$fdr_deplete < 0.05, na.rm = TRUE))
})) |>
  summarise(across(c(n_tests, enriched, depleted), sum), .by = cohort) |>
  pivot_longer(c(enriched, depleted), names_to = "dir", values_to = "n") |>
  mutate(rate = n / n_tests,
         dir  = factor(dir, c("enriched", "depleted"),
                       labels = c("Marker-enriched", "Marker-depleted")),
         cohort = factor(cohort, names(COHORT_SHAPE)))

pC <- ggplot(mk, aes(dir, rate, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = NA) +
  geom_text(aes(label = sprintf("%d/%d", n, n_tests)),
            position = position_dodge(width = 0.72), vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = c(BrainSEQ = "#E69F00", GTEx = "#0072B2"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.019),
                     labels = scales::percent_format(accuracy = 0.5),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Module x cell-type tests\nsignificant (FDR < 0.05)") +
  theme_pub() + theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Assemble: A | B on top (the result), C beneath at half width (the control)
# ---------------------------------------------------------------------------
design <- "AABB
           AABB
           CCBB"
fig <- wrap_plots(pA, pB, pC, design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figCompositionRobustness", width = 7.2, height = 5.2)
cat("Done. Output in", FIG_DIR, "\n")
