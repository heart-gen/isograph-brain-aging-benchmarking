# Supplementary figure: module-level colocalization convergence, all five traits.
#
# Fig 4E showed per-module SCZ coloc counts and nothing else. Counts alone cannot
# support a convergence claim: a big module collects more hits than a small one for no
# biological reason, and a gene can only colocalize if it was TESTED. This figure adds
# the two comparisons the count panel was missing, and extends it from SCZ to the four
# aging traits that are the paper's actual subject.
#
# (A) Global concentration against a size-matched null. The statistic is the sum of
#     squared per-module hit counts; the null redraws the same number of hit genes from
#     the same tested pool, so module sizes are held fixed by construction. NO trait, in
#     either partition, is more concentrated than chance.
# (B) Per-module counts for all five traits — the extended Fig 4E content, kept so the
#     reader can see what the underlying data look like, with modules too small in the
#     tested pool to be testable drawn as open symbols.
# (C) Why the published SCZ enrichment does not survive. Coloc genes look enriched in
#     MAGMA-anchored modules only when the denominator is ALL module genes. Anchored
#     modules are MAGMA-enriched for the same GWAS that decides which genes enter the
#     coloc test, so the tested pool is already anchored-rich; against it the enrichment
#     is null. Both denominators are plotted so the reader can see the ascertainment.
#
# Reads 05_genetic_anchoring/_m/module_coloc_convergence/{convergence,global}.parquet.
# Writes manuscript/_m/figures/figColocConvergence.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/coloc_convergence_figure.R
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
CONV    <- rel("05_genetic_anchoring", "_m", "module_coloc_convergence")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

OBS_COL  <- "#D55E00"   # observed
NULL_COL <- "grey45"    # size-matched null
COHORT_COL <- c(GTEx = "#0072B2", BrainSEQ = "#E69F00")
COHORT_LAB <- c(gtex_caudate_bg = "GTEx", brainseq_caudate = "BrainSEQ")
TRAIT_ORDER <- c("AD", "PD", "LBD", "ALS", "SCZ")

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

g <- as.data.frame(read_parquet(file.path(CONV, "global.parquet"))) |>
  mutate(cohort = unname(COHORT_LAB[source]),
         trait  = factor(trait, TRAIT_ORDER),
         lab    = paste0(trait, " / ", cohort))
c_ <- as.data.frame(read_parquet(file.path(CONV, "convergence.parquet"))) |>
  mutate(cohort = unname(COHORT_LAB[source]),
         trait  = factor(trait, TRAIT_ORDER))

cat(sprintf("  %d (trait, source) cells; concentration p range %.3f-%.3f; none < 0.05\n",
            nrow(g), min(g$concentration_p, na.rm = TRUE),
            max(g$concentration_p, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Panel A - observed concentration vs the size-matched null
# ---------------------------------------------------------------------------
gA <- g |> filter(n_coloc_genes > 0) |> arrange(trait, cohort)
gA$lab <- factor(gA$lab, levels = rev(gA$lab))
gA$plab <- sprintf("P = %.2f", gA$concentration_p)
# Headroom on the right for the P labels, computed from the data so a larger future
# value cannot push the label off-panel.
xmax <- max(c(gA$concentration_obs, gA$concentration_null_mean), na.rm = TRUE) * 1.42

pA <- ggplot(gA, aes(y = lab)) +
  geom_segment(aes(x = concentration_null_mean, xend = concentration_obs, yend = lab),
               linewidth = 0.4, colour = "grey70") +
  geom_point(aes(x = concentration_null_mean, shape = "Size-matched null"),
             colour = NULL_COL, size = 1.9) +
  geom_point(aes(x = concentration_obs, shape = "Observed"),
             colour = OBS_COL, size = 2.2) +
  geom_text(aes(x = xmax, label = plab), hjust = 1, size = 2.2, colour = "grey30") +
  scale_shape_manual(values = c(Observed = 16, `Size-matched null` = 18), name = NULL) +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0.01, 0))) +
  labs(x = "Concentration of coloc genes across modules\n(sum of squared per-module counts)",
       y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel B - per-module counts, all five traits (extended Fig 4E content)
# ---------------------------------------------------------------------------
# Untestable = too few genes in the tested pool for the hypergeometric to mean anything
# (a 1-of-1 module returns a "significant" p purely because it is tiny). Drawn open
# rather than dropped, so the reader sees them and does not read them as evidence.
cB <- c_ |>
  mutate(ylab = paste0(cohort, " ", module_id)) |>
  arrange(trait, n_coloc_genes)
cB$row <- seq_len(nrow(cB))

pB <- ggplot(cB, aes(n_coloc_genes, reorder(paste(row, ylab), row), colour = cohort)) +
  geom_segment(aes(x = 0, xend = n_coloc_genes,
                   yend = reorder(paste(row, ylab), row)), linewidth = 0.45) +
  geom_point(aes(shape = testable), size = 1.7) +
  facet_grid(trait ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = COHORT_COL, name = "Aging cohort") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "testable", `FALSE` = "too few genes in pool"),
                     name = NULL) +
  scale_y_discrete(labels = function(x) sub("^[0-9]+ ", "", x)) +
  # limit derived from the data: a hardcoded bound silently DROPPED the two
  # 6-gene SCZ modules on the previous build ("Removed 2 rows ... outside the
  # scale range"), which is the one failure mode a count panel must not have.
  scale_x_continuous(breaks = seq(0, max(cB$n_coloc_genes), 1),
                     limits = c(0, max(cB$n_coloc_genes) * 1.07),
                     expand = expansion(mult = c(0.01, 0))) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2)),
         shape  = guide_legend(order = 2)) +
  labs(x = "Colocalizing switch genes in module", y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 5.6),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel C - the denominator, which is the whole result
# ---------------------------------------------------------------------------
cC <- g |>
  filter(n_coloc_genes > 0) |>
  select(lab, trait, cohort, frac_coloc_anchored,
         frac_pool_anchored, frac_allgenes_anchored) |>
  pivot_longer(c(frac_pool_anchored, frac_allgenes_anchored),
               names_to = "denom", values_to = "background") |>
  mutate(denom = factor(denom,
                        levels = c("frac_allgenes_anchored", "frac_pool_anchored"),
                        labels = c("all module genes (as published)",
                                   "CLPP-tested pool (correct)")))
cC$lab <- factor(cC$lab, levels = rev(unique(g$lab[g$n_coloc_genes > 0])))

pC <- ggplot(cC, aes(y = lab)) +
  geom_segment(aes(x = background, xend = frac_coloc_anchored, yend = lab),
               linewidth = 0.4, colour = "grey78") +
  geom_point(aes(x = background, colour = denom), size = 1.9) +
  geom_point(aes(x = frac_coloc_anchored, shape = "Observed coloc genes"),
             colour = OBS_COL, size = 2.1) +
  scale_colour_manual(values = c(`all module genes (as published)` = "#999999",
                                 `CLPP-tested pool (correct)` = "#009E73"),
                      name = "Background (denominator)") +
  scale_shape_manual(values = c(`Observed coloc genes` = 16), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.72), expand = expansion(mult = c(0.01, 0.02))) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 2)) +
  labs(x = "Fraction in MAGMA-anchored modules", y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

design <- "AABBB
           CCBBB"
fig <- wrap_plots(pA, pB, pC, design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figColocConvergence", width = 7.2, height = 7.4)
cat("Done. Output in", FIG_DIR, "\n")
