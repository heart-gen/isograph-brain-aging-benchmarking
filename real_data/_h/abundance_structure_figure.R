# Supplementary real-data figure: abundance and isoform structure carry partially
# NON-REDUNDANT information — the two channels can be separated computationally and the
# separation adds information. (A) per-gene abundance-vs-switch axis correlation piles up
# near 0 -> the inferred axes are largely orthogonal. (B) the de-confounded incremental
# test finds, in every cohort/region, a small specific set of composition-unique genes
# whose switch channel carries phenotype signal that total abundance cannot -> separation
# adds information. (C) one such gene: total abundance flat across diagnosis while the
# isoform-switch score shifts. Reads
# 02_module_discovery/brainseq/caudate_sczd/_m/isograph_vae/abundance_structure/*, writes
# figSeparation.{pdf,png} to real_data/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        real_data/_h/abundance_structure_figure.R
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
AS_DIR  <- rel("real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae", "abundance_structure")
FIG_DIR <- rel("real_data", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito. IsoGraph vermillion is reused for the switch channel throughout the manuscript.
ISO_COL    <- "#D55E00"   # switch / composition-unique
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

ortho <- as.data.frame(read_parquet(file.path(AS_DIR, "axis_orthogonality.parquet")))
incr  <- as.data.frame(read_parquet(file.path(AS_DIR, "incremental_summary.parquet")))
eg    <- as.data.frame(read_parquet(file.path(AS_DIR, "example_gene.parquet")))
egs   <- fromJSON(file.path(AS_DIR, "example_gene_stats.json"))

# ---------------------------------------------------------------------------
# Panel A - abundance vs switch axes are largely orthogonal per gene
# ---------------------------------------------------------------------------
med_abs <- median(ortho$abs_r)
frac_lt <- mean(ortho$abs_r < 0.1)
pA <- ggplot(ortho, aes(pearson_r)) +
  geom_histogram(bins = 60, fill = ISO_COL, colour = NA, alpha = 0.85,
                 boundary = 0) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey35") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  annotate("text", x = -0.98, y = Inf, hjust = 0, vjust = 1.4, size = 2.5,
           colour = "grey25",
           label = sprintf("median |r| = %.2f\n%.0f%% of genes |r| < 0.1",
                           med_abs, 100 * frac_lt)) +
  labs(x = "Per-gene correlation of abundance vs switch axis (Pearson r)",
       y = "Genes") +
  theme_pub()

# ---------------------------------------------------------------------------
# Panel B - composition-unique genes: switch adds phenotype signal beyond abundance
# ---------------------------------------------------------------------------
clean_region <- function(x) {
  x <- gsub("_basal_ganglia", "", x)
  x <- gsub("_c_1", "", x)
  gsub("_", " ", x)
}
COHORT_TAG <- c(`BrainSEQ SCZD` = "SCZD", `BrainSEQ aging` = "BSeq",
                `GTEx aging` = "GTEx")
datB <- incr |>
  mutate(region_lab = sprintf("%s: %s", COHORT_TAG[cohort], clean_region(region)),
         cohort = factor(cohort, names(COHORT_COL))) |>
  arrange(composition_unique) |>
  mutate(row_lab = factor(region_lab, levels = region_lab))
pB <- ggplot(datB, aes(composition_unique, row_lab, fill = cohort)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_text(aes(label = composition_unique), hjust = -0.25, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = COHORT_COL) +
  scale_x_sqrt(expand = expansion(mult = c(0, 0.10)),
               breaks = c(0, 10, 50, 150, 300, 545)) +
  labs(x = "Composition-unique genes (switch-sig, abundance-not)", y = NULL) +
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
  theme_pub() + theme(panel.spacing = unit(8, "pt"))

fig <- ((pA | pC) + plot_layout(widths = c(1, 1.05))) / pB +
  plot_layout(heights = c(0.95, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figSeparation", width = 7.2, height = 5.8)
cat("Done. Output in", FIG_DIR, "\n")
