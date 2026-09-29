# Module context and held-out DTU evidence (supplementary figure).
#
# The question is whether a gene's IsoGraph module neighbours, scored on one donor half,
# carry information about that gene's satuRn DTU evidence in the other half beyond the gene's
# own discovery evidence -- and beyond its 50 most co-expressed genes. It is a complementary
# test, not a comparison with gene-wise DTU, and the figure is built to show the
# heterogeneity rather than a pooled effect: three of the six analyses carry the signal.
#
# (A) Partial correlation per split direction (5 seeds x 2 directions), before and after the
#     co-expression comparator; the large mark is the mean, which is the permutation
#     statistic. q is the BH-adjusted stratified permutation q of the adjusted arm.
# (B) Its practical reading: among genes that did not meet the discovery DTU criterion, the
#     held-out replication rate (P < 0.05) for the bottom vs top tertile of module context,
#     with the median adjusted OR over directions. Descriptive; (A) is the inference.
#
# Reads 04_module_trust/_m/dtu_added_value/{region_summary,replicates}.csv (CSV, not parquet,
# so the figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figDtuAddedValue.{pdf,png}.
# Run: Rscript manuscript/_h/dtu_added_value_figure.R
suppressPackageStartupMessages({
  library(dplyr)
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
DTU     <- rel("04_module_trust", "_m", "dtu_added_value")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ISOGRAPH_COL <- "#D55E00"   # Okabe-Ito vermillion, as everywhere else in the project
PLAIN_COL    <- "grey55"

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.y = element_blank(),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(9, "pt"),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}

need <- function(f) {
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun: python -m isograph_benchmark.real_data.module_dtu_added_value summarize",
         call. = FALSE)
  }
  read.csv(f)
}
s   <- need(file.path(DTU, "region_summary.csv"))
rep <- need(file.path(DTU, "replicates.csv"))

# Anatomical pairs next to each other, BrainSEQ above GTEx within a pair.
LABELS <- c(
  "brainseq:caudate"              = "Caudate (BrainSEQ)",
  "gtex:caudate_basal_ganglia"    = "Caudate (GTEx)",
  "brainseq:hippocampus"          = "Hippocampus (BrainSEQ)",
  "gtex:hippocampus"              = "Hippocampus (GTEx)",
  "brainseq:dlpfc"                = "DLPFC (BrainSEQ)",
  "gtex:frontal_cortex_ba9"       = "Frontal ctx BA9 (GTEx)"
)
lab <- function(cohort, region) {
  factor(unname(LABELS[paste(cohort, region, sep = ":")]), levels = rev(unname(LABELS)))
}
s   <- s   |> mutate(analysis = lab(cohort, region))
rep <- rep |> mutate(analysis = lab(cohort, region))
stopifnot(!anyNA(s$analysis), !anyNA(rep$analysis))
n_dir <- max(s$n_replicates)

ARMS <- c("Module context", "Module context | co-expression")
arm_cols <- setNames(c(PLAIN_COL, ISOGRAPH_COL), ARMS)
nudge <- setNames(c(0.17, -0.17), ARMS)

pts <- bind_rows(
  rep |> transmute(analysis, arm = ARMS[1], r = r_module),
  rep |> transmute(analysis, arm = ARMS[2], r = r_module_given_abund)
) |>
  mutate(arm = factor(arm, ARMS), y = as.numeric(analysis) + nudge[as.character(arm)])
means <- bind_rows(
  s |> transmute(analysis, arm = ARMS[1], r = r_module),
  s |> transmute(analysis, arm = ARMS[2], r = r_module_given_abund)
) |>
  mutate(arm = factor(arm, ARMS), y = as.numeric(analysis) + nudge[as.character(arm)])

fmt_q <- function(q) ifelse(q < 0.01, sprintf("q = %.3f", q), sprintf("q = %.2f", q))
qlab <- s |>
  transmute(analysis, y = as.numeric(analysis),
            label = sprintf("%s  %d/%d", fmt_q(q_given_abund_strat),
                            r_module_given_abund_n_positive, n_dir),
            strong = q_given_abund_strat < 0.05)

x_hi <- max(pts$r) + 0.02
pA <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey40") +
  geom_point(data = pts, aes(r, y, colour = arm), size = 1.1, alpha = 0.55,
             shape = 16, position = position_jitter(height = 0.05, width = 0, seed = 7)) +
  geom_point(data = means, aes(r, y, fill = arm), shape = 23, size = 2.4,
             colour = "black", stroke = 0.35) +
  geom_text(data = qlab, aes(x_hi, y, label = label, fontface = ifelse(strong, "bold", "plain")),
            hjust = 0, size = 2.5, colour = "grey15") +
  scale_colour_manual(values = arm_cols, breaks = ARMS) +
  scale_fill_manual(values = arm_cols, breaks = ARMS) +
  scale_y_continuous(breaks = seq_along(levels(s$analysis)), labels = levels(s$analysis),
                     expand = expansion(add = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.02))) +
  coord_cartesian(clip = "off") +
  guides(colour = "none",
         fill = guide_legend(override.aes = list(size = 2.4), nrow = 1)) +
  labs(x = sprintf("Partial correlation with held-out DTU evidence\n(%d split directions per analysis)",
                   n_dir),
       y = NULL) +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 78, 4, 4, "pt"))

# (B) sub-threshold replication, bottom -> top tertile of module context.
sub <- s |>
  transmute(analysis, y = as.numeric(analysis),
            bottom = subthr_rate_bottom, top = subthr_rate_top,
            label = sprintf("OR %.2f  %d/%d", subthr_or_median, subthr_or_n_above_1, n_dir))
TERT <- c("Bottom tertile", "Top tertile")
sub_long <- bind_rows(
  sub |> transmute(analysis, y, rate = bottom, tertile = TERT[1]),
  sub |> transmute(analysis, y, rate = top, tertile = TERT[2])
) |> mutate(tertile = factor(tertile, TERT))
x_hi_b <- max(sub_long$rate) + 0.02
pB <- ggplot(sub) +
  geom_segment(aes(x = bottom, xend = top, y = y, yend = y), linewidth = 0.6, colour = "grey70") +
  geom_point(data = sub_long, aes(rate, y, fill = tertile), shape = 21, size = 2.3,
             colour = "black", stroke = 0.3) +
  geom_text(aes(x_hi_b, y, label = label), hjust = 0, size = 2.5, colour = "grey15") +
  scale_fill_manual(values = setNames(c("white", ISOGRAPH_COL), TERT)) +
  scale_y_continuous(breaks = seq_along(levels(s$analysis)), labels = NULL,
                     expand = expansion(add = 0.5)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Held-out replication (P < 0.05)\namong sub-threshold genes", y = NULL) +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        axis.ticks.y = element_blank(),
        plot.margin = margin(4, 62, 4, 4, "pt"))

fig <- (pA | pB) + plot_layout(widths = c(1.45, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

W <- 7.1; H <- 3.2
ggsave(file.path(FIG_DIR, "figDtuAddedValue.pdf"), fig, width = W, height = H,
       units = "in", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "figDtuAddedValue.png"), fig, width = W, height = H,
       units = "in", dpi = 300)
cat("  figDtuAddedValue saved\n")
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
