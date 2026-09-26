# Supplementary real-data figure: abundance and isoform structure carry partially
# NON-REDUNDANT information — the two channels can be separated computationally and the
# separation adds information. (A) per-gene abundance-vs-switch axis correlation piles up
# near 0 in EVERY one of the 17 analyses -> the inferred axes are largely orthogonal, and
# that is a property of the representation, not of one dataset. (B) the de-confounded
# incremental test finds, in every cohort/region, a small specific set of switch-unique
# genes whose switch channel carries phenotype signal that total abundance cannot ->
# separation adds information. (C) one such gene: total abundance flat across diagnosis
# while the isoform-switch score shifts. Reads
# 03_module_characterization/_m/axis_orthogonality_{all,summary}.parquet (the _h/03b +
# _h/03c rollup over all stores) for A and
# 02_module_discovery/brainseq/caudate_sczd/_m/isograph_vae/abundance_structure/* (_h/03a)
# for B-C; writes figSeparation.{pdf,png} to manuscript/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/abundance_structure_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
AS_DIR  <- rel("02_module_discovery", "brainseq", "caudate_sczd", "_m", "isograph_vae", "abundance_structure")
ROLLUP  <- rel("03_module_characterization", "_m")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito. IsoGraph vermillion is reused for the switch channel throughout the manuscript.
ISO_COL    <- "#D55E00"   # switch / switch-unique
AB_COL     <- "#999999"   # abundance channel
COHORT_COL <- c(`BrainSEQ SCZD` = "#D55E00", `BrainSEQ aging` = "#E69F00",
                `GTEx aging` = "#0072B2")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.34, "cm"),
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

pfmt <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 1e-4) return(sprintf("p = %.0e", p))
  sprintf("p = %.3g", p)
}

ortho <- as.data.frame(read_parquet(file.path(ROLLUP, "axis_orthogonality_all.parquet")))
osum  <- as.data.frame(read_parquet(file.path(ROLLUP, "axis_orthogonality_summary.parquet")))
incr  <- as.data.frame(read_parquet(file.path(AS_DIR, "incremental_summary.parquet")))
eg    <- as.data.frame(read_parquet(file.path(AS_DIR, "example_gene.parquet")))
egs   <- fromJSON(file.path(AS_DIR, "example_gene_stats.json"))

# ---------------------------------------------------------------------------
# Panel A - abundance vs switch axes are largely orthogonal per gene, in every analysis
# ---------------------------------------------------------------------------
clean_region <- function(x) {
  x <- gsub("_basal_ganglia", "", x)
  x <- gsub("_c_1", "", x)
  gsub("_", " ", x)
}
COHORT_TAG <- c(`BrainSEQ SCZD` = "SCZD", `BrainSEQ aging` = "BSeq",
                `GTEx aging` = "GTEx")
# Short two-line strip labels: 6 facets share 7.2 in, so ~20 characters per line.
short_region <- function(x) {
  m <- c(anterior_cingulate_cortex_ba24 = "ant. cingulate BA24",
         frontal_cortex_ba9 = "frontal cortex BA9",
         cerebellar_hemisphere = "cerebellar hemi.",
         spinal_cord_cervical_c_1 = "spinal cord", dlpfc = "DLPFC")
  ifelse(x %in% names(m), m[x], clean_region(x))
}
# facet order = rollup order (SCZD, BrainSEQ aging regions, GTEx regions)
osum <- osum |>
  mutate(facet = sprintf("%s\n%s", COHORT_TAG[cohort], short_region(region)),
         facet = factor(facet, levels = unique(facet)),
         cohort = factor(cohort, names(COHORT_COL)),
         lab = sprintf("med |r| %.2f\n%.0f%% < 0.1", median_abs_r, 100 * frac_abs_r_lt_0.1))
datA <- ortho |>
  mutate(facet = factor(sprintf("%s\n%s", COHORT_TAG[cohort], short_region(region)),
                        levels = levels(osum$facet)),
         cohort = factor(cohort, names(COHORT_COL)))
stopifnot(!anyNA(datA$facet), nlevels(osum$facet) == nrow(osum))
cat(sprintf("  panel A: %d analyses, median |r| %.3f-%.3f, frac |r|<0.1 %.2f-%.2f\n",
            nrow(osum), min(osum$median_abs_r), max(osum$median_abs_r),
            min(osum$frac_abs_r_lt_0.1), max(osum$frac_abs_r_lt_0.1)))
# y = fraction of the panel's genes (after_stat sums within each facet), so one shared
# axis serves all 17 panels and the BrainSEQ/GTEx gene-count difference does not show.
pA <- ggplot(datA, aes(pearson_r, fill = cohort)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), bins = 40, colour = NA,
                 alpha = 0.85, boundary = 0) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey35") +
  geom_text(data = osum, aes(x = -0.97, y = Inf, label = lab), hjust = 0, vjust = 1.25,
            size = 2.0, lineheight = 0.9, colour = "grey25", inherit.aes = FALSE) +
  facet_wrap(~ facet, ncol = 6) +
  scale_fill_manual(values = COHORT_COL, guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.45)), breaks = c(0, 0.005, 0.01, 0.015),
                     labels = function(v) sprintf("%.1f", 100 * v)) +
  labs(x = "Per-gene correlation of abundance vs switch axis (Pearson r)",
       y = "Genes (% of analysis)") +
  theme_pub() +
  theme(strip.text = element_text(size = 6.5, lineheight = 0.9,
                                  margin = margin(1, 0, 1, 0, "pt")),
        axis.text = element_text(size = 6.2),
        panel.spacing.x = unit(9, "pt"),
        panel.spacing.y = unit(4, "pt"))

# ---------------------------------------------------------------------------
# Panel B - switch-unique genes: switch adds phenotype signal beyond abundance
# ---------------------------------------------------------------------------
datB <- incr |>
  mutate(region_lab = sprintf("%s: %s", COHORT_TAG[cohort], clean_region(region)),
         cohort = factor(cohort, names(COHORT_COL))) |>
  arrange(composition_unique) |>
  mutate(row_lab = factor(region_lab, levels = region_lab))
pB <- ggplot(datB, aes(composition_unique, row_lab, fill = cohort)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_text(aes(label = composition_unique), hjust = -0.25, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = COHORT_COL) +
  scale_x_sqrt(expand = expansion(mult = c(0, 0.18)),
               breaks = c(0, 10, 50, 150, 300, 545)) +
  # "switch-unique", not "composition-unique": the latter reads as cell composition,
  # which is the confounder this layer is defended against, not what the class means.
  labs(x = "Switch-unique genes (switch-significant, abundance not)", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
        axis.text.y = element_text(size = 6.8),
        legend.position = c(0.82, 0.30))

# ---------------------------------------------------------------------------
# Panel C - example gene: stable total abundance, real isoform switch
# ---------------------------------------------------------------------------
gene_lab <- sub("\\..*$", "", egs$gene_id)
if (egs$group_col == "Dx") {
  lv <- intersect(c("Control", "SCZD"), unique(eg$group))
  eg$group <- factor(eg$group, levels = c(lv, setdiff(unique(eg$group), lv)))
}
CH_LAB <- c(abundance_z = "Abundance (z)", switch_score = "Switch score")
egl <- eg |>
  pivot_longer(c(abundance_z, switch_score), names_to = "channel", values_to = "value") |>
  mutate(channel = factor(CH_LAB[channel], levels = unname(CH_LAB)))
labC <- data.frame(
  channel = factor(unname(CH_LAB), levels = unname(CH_LAB)),
  lab = c(pfmt(egs$p_abund_given_switch), pfmt(egs$p_switch_given_abund)))
pC <- ggplot(egl, aes(group, value)) +
  geom_boxplot(aes(fill = channel), width = 0.6, outlier.shape = NA, alpha = 0.85,
               linewidth = 0.3) +
  geom_jitter(width = 0.14, height = 0, size = 0.5, alpha = 0.35, colour = "grey20") +
  facet_wrap(~ channel, scales = "free_y") +
  geom_text(data = labC, aes(x = 0.62, y = Inf, label = lab), hjust = 0, vjust = 1.5,
            size = 2.6, colour = "grey25", inherit.aes = FALSE) +
  scale_fill_manual(values = c(`Abundance (z)` = AB_COL,
                               `Switch score` = ISO_COL), guide = "none") +
  labs(x = NULL, y = gene_lab) +
  theme_pub() + theme(panel.spacing = unit(8, "pt"), strip.text = element_text(size = 7.5))

# free() releases A from B's long axis labels (A is otherwise squeezed to B's panel).
pA_cell <- if (exists("free", where = asNamespace("patchwork"), inherits = FALSE)) patchwork::free(pA) else pA
fig <- pA_cell / ((pB | pC) + plot_layout(widths = c(1.25, 1))) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figSeparation", width = 7.2, height = 7.0)
cat("Done. Output in", FIG_DIR, "\n")
