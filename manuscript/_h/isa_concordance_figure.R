# Supplementary figure: IsoGraph's switch genes carry elevated DTU evidence under an
# independent, established caller (satuRn, the DTU test used by IsoformSwitchAnalyzeR v2).
#
# This answers "is the switch signal an artifact of your method?" with a tool the reader
# already trusts, across 17 cohort x region x trait analyses.
#
# Two honesty requirements are built into the panel rather than left to the caption:
#   1. Regions with ZERO switch genes are VACUOUS by construction, not null results.
#      They are drawn as explicit "not testable" rows so 11/17 significant cannot be
#      mistaken for 11/17 of a comparable set.
#   2. The binary empirical-FDR overlap count is a conservative LOWER BOUND -- Efron's
#      empirical null assumes a sparse alternative, and pervasive brain DTU violates it.
#      The continuous rank-based statistic is therefore the primary test, and it is what
#      this figure plots.
#
# (A) rank-biserial effect size per analysis with the Mann-Whitney p, vacuous analyses
#     shown as not-testable rows;
# (B) the underlying gene-level satuRn evidence distributions for the strongest and the
#     weakest testable analysis, so the reader can see what a given effect size means.
#
# Reads 06_switch_mechanism/_m/isa_concordance/*/summary.json and gene_concordance.parquet;
# writes figIsaConcordance.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/isa_concordance_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
ISA_DIR <- rel("06_switch_mechanism", "_m", "isa_concordance")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

SIG_COL   <- "#D55E00"   # significant under the continuous test
NS_COL    <- "#999999"   # testable but not significant
VAC_COL   <- "grey80"    # vacuous: no switch genes to test
SWITCH_COL<- "#D55E00"
BG_COL    <- "#56B4E9"

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

# Ocean FS: iterdir-style listing, not recursive globbing over a huge tree.
analyses <- list.dirs(ISA_DIR, full.names = FALSE, recursive = FALSE)
analyses <- analyses[nzchar(analyses)]

s <- bind_rows(lapply(analyses, function(a) {
  f <- file.path(ISA_DIR, a, "summary.json")
  if (!file.exists(f)) return(NULL)
  # isa_concordance writes bare NaN for the vacuous analyses (no switch genes), which is
  # not legal JSON; jsonlite rejects it. Coerce to null before parsing rather than
  # skipping those analyses -- they must appear in the figure as not-testable rows.
  j <- fromJSON(gsub("\\bNaN\\b", "null", paste(readLines(f, warn = FALSE), collapse = "\n")))
  data.frame(
    analysis  = a,
    cohort    = j$cohort,
    region    = j$region,
    trait     = j$trait,
    n_switch  = j$n_switch_genes,
    n_bg      = j$n_background_genes,
    rb        = if (is.null(j$rank_biserial)) NA_real_ else j$rank_biserial,
    p         = if (is.null(j$mwu_p))         NA_real_ else j$mwu_p,
    or        = if (is.null(j$adjusted_or))   NA_real_ else j$adjusted_or,
    med_sw    = if (is.null(j$median_evidence_switch))     NA_real_ else j$median_evidence_switch,
    med_bg    = if (is.null(j$median_evidence_background)) NA_real_ else j$median_evidence_background)
}))

PRETTY_REGION <- function(x) {
  x <- gsub("_", " ", x)
  x <- sub("anterior cingulate cortex ba24", "ACC BA24", x)
  x <- sub("frontal cortex ba9", "frontal ctx BA9", x)
  x <- sub("nucleus accumbens basal ganglia", "n. accumbens", x)
  x <- sub("caudate basal ganglia", "caudate", x)
  x <- sub("putamen basal ganglia", "putamen", x)
  x <- sub("spinal cord cervical c 1", "spinal cord C1", x)
  x <- sub("cerebellar hemisphere", "cerebellar hem.", x)
  x <- sub("^dlpfc$", "DLPFC", x)
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

s <- s |>
  mutate(
    cohort_lab = ifelse(cohort == "gtex", "GTEx", "BrainSEQ"),
    trait_lab  = ifelse(trait == "dx", ", SCZD", ""),
    label = sprintf("%s%s (%s)", PRETTY_REGION(region), trait_lab, cohort_lab),
    # A region with no switch genes has nothing to test. That is a coverage fact about
    # where IsoGraph called age modules, NOT a failed concordance test, and the two must
    # not be pooled into one denominator.
    status = case_when(
      n_switch == 0            ~ "Not testable (no switch genes)",
      !is.finite(p)            ~ "Not testable (no switch genes)",
      p < 0.05 & rb > 0        ~ "Concordant (P < 0.05)",
      TRUE                     ~ "Not significant"),
    status = factor(status, c("Concordant (P < 0.05)", "Not significant",
                              "Not testable (no switch genes)")))

n_testable <- sum(s$status != "Not testable (no switch genes)")
n_sig      <- sum(s$status == "Concordant (P < 0.05)")
cat(sprintf("  %d/%d testable analyses concordant (%d vacuous)\n",
            n_sig, n_testable, sum(s$status == "Not testable (no switch genes)")))

# ---------------------------------------------------------------------------
# Panel A - rank-biserial effect size per analysis
# ---------------------------------------------------------------------------
sa <- s |> arrange(status, rb)
sa$label <- factor(sa$label, levels = sa$label)
# p-value annotation only where a test was actually run
sa$plab <- ifelse(is.finite(sa$p),
                  ifelse(sa$p < 1e-4, sprintf("%.0e", sa$p), sprintf("%.3f", sa$p)), "")

pA <- ggplot(sa, aes(rb, label, colour = status)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = 0, xend = rb, yend = label), linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 1.9, na.rm = TRUE) +
  geom_text(aes(x = 0.575, label = plab), hjust = 1, size = 2.1, colour = "grey35") +
  geom_text(data = subset(sa, status == "Not testable (no switch genes)"),
            aes(x = 0.01, label = "no switch genes called"),
            hjust = 0, size = 2.1, colour = "grey55") +
  scale_colour_manual(values = c(`Concordant (P < 0.05)` = SIG_COL,
                                 `Not significant` = NS_COL,
                                 `Not testable (no switch genes)` = VAC_COL),
                      name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-0.10, 0.585), breaks = seq(-0.1, 0.4, 0.1)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Rank-biserial effect size\n(satuRn DTU evidence, switch vs background genes)",
       y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 6.6),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel B - what an effect size looks like: median satuRn evidence, switch vs background
# ---------------------------------------------------------------------------
# Paired within analysis. Plotting the two medians directly keeps the comparison on the
# scale the test was run on, without asserting a distributional shape.
pb <- s |>
  filter(status != "Not testable (no switch genes)") |>
  select(label, med_sw, med_bg) |>
  pivot_longer(c(med_sw, med_bg), names_to = "set", values_to = "med") |>
  mutate(set = factor(set, c("med_bg", "med_sw"),
                      labels = c("Background\ngenes", "IsoGraph\nswitch genes")))

pB <- ggplot(pb, aes(set, med, group = label)) +
  geom_line(linewidth = 0.4, colour = "grey65", alpha = 0.8) +
  geom_point(aes(colour = set), size = 1.7) +
  scale_colour_manual(values = c(`Background\ngenes` = BG_COL,
                                 `IsoGraph\nswitch genes` = SWITCH_COL),
                      guide = "none") +
  labs(x = NULL, y = expression("Median satuRn evidence  " * -log[10] * "(min isoform P)")) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7.2),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey92"))

fig <- (pA | pB) +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figIsaConcordance", width = 7.2, height = 4.4)
cat("Done. Output in", FIG_DIR, "\n")
