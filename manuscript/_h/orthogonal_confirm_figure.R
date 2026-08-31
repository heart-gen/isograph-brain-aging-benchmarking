# Supplementary figure: the genetically anchored switch pairs behave like switches on an
# independent long-read platform.
#
# Every switch call in this paper otherwise rests on short-read quantification, so
# "your switches are quantification artifacts" is an objection no other panel answers.
# This one does, on an independent platform, lab, cohort and quantifier (Aguzzoli-Heberle
# 2024 NBT; ONT DLPFC BA9/46; Bambu; n = 12), against an abundance-matched null drawn
# from 67,167 non-splicing-led IsoGraph switch pairs.
#
# (A) observed switch-like rate vs the matched null, for all detected anchored pairs and
#     for the abundance-qualified subset; the unmatched background rate is drawn as a
#     reference so the reader can see what the matching is correcting for.
# (B) per-gene confirmation for the twelve splicing-led genes, with the genes that FAIL
#     the abundance qualification marked rather than dropped -- their anchored isoform is
#     too lowly expressed in long-read for any usage correlation to be interpretable.
#
# Matching is on the abundance decile of the better-expressed pair member, because at
# n = 12 that is what governs whether a usage correlation is estimable at all.
#
# Reads 06_switch_mechanism/_m/switch_orthogonal_confirm/{anchored_summary.json,
#       anchored_gene_confirmation.parquet}; writes figOrthogonalConfirm.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/orthogonal_confirm_figure.R
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
OC_DIR  <- rel("06_switch_mechanism", "_m", "switch_orthogonal_confirm")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, shared with the other real-data figures.
OBS_COL  <- "#D55E00"   # observed anchored set
NULL_COL <- "grey45"    # abundance-matched null
FAIL_COL <- "#999999"   # fails the abundance qualification

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

s <- fromJSON(file.path(OC_DIR, "anchored_summary.json"))

# ---------------------------------------------------------------------------
# Panel A - observed switch-like rate against the abundance-matched null
# ---------------------------------------------------------------------------
# Every number is read from the summary JSON rather than transcribed, so the annotation
# cannot drift from switch_orthogonal_confirm.py's output.
mk_row <- function(key, label) {
  m <- s$matched_null[[key]]
  data.frame(set = label, n = m$n_focal, observed = m$observed,
             null_mean = m$null_mean, lo = m$null_q025, hi = m$null_q975,
             p = m$p_empirical_two_sided)
}
a <- bind_rows(
  mk_row("switch_like_rate",
         sprintf("All detected pairs\n(n = %d)", s$matched_null$switch_like_rate$n_focal)),
  mk_row("switch_like_rate_usable_only",
         sprintf("Anchored isoform usably\nexpressed (n = %d)",
                 s$matched_null$switch_like_rate_usable_only$n_focal)))
a$set <- factor(a$set, levels = a$set)
a$plab <- sprintf("P = %s", formatC(a$p, format = "g", digits = 2))

pA <- ggplot(a, aes(set)) +
  # unmatched background: what the rate would be compared against WITHOUT matching
  geom_hline(yintercept = s$background$switch_like_rate_unmatched,
             linewidth = 0.35, linetype = "dotted", colour = "grey45") +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = NULL_COL,
                 linewidth = 2.6, alpha = 0.35) +
  geom_point(aes(y = null_mean), colour = NULL_COL, size = 1.9, shape = 18) +
  geom_point(aes(y = observed), colour = OBS_COL, size = 2.6) +
  geom_text(aes(y = observed, label = sprintf("%.2f", observed)),
            colour = OBS_COL, size = 2.5, hjust = -0.45) +
  geom_text(aes(y = hi, label = plab), size = 2.4, colour = "grey25",
            vjust = -0.9) +
  annotate("text", x = 0.55, y = s$background$switch_like_rate_unmatched,
           label = "unmatched background", hjust = 0, vjust = -0.5,
           size = 2.2, colour = "grey45") +
  scale_y_continuous(limits = c(0, 0.72), breaks = seq(0, 0.7, 0.1)) +
  labs(x = NULL, y = "Switch-like rate in long-read") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

# ---------------------------------------------------------------------------
# Panel B - per-gene confirmation, with the abundance-disqualified genes named
# ---------------------------------------------------------------------------
g <- as.data.frame(read_parquet(file.path(OC_DIR, "anchored_gene_confirmation.parquet"))) |>
  mutate(status = case_when(
           confirmed_at_usable_abundance ~ "Confirmed",
           orthogonally_confirmed        ~ "Anchored isoform too lowly expressed",
           TRUE                          ~ "Not confirmed"),
         status = factor(status, c("Confirmed",
                                   "Anchored isoform too lowly expressed",
                                   "Not confirmed")),
         gene_name = reorder(gene_name, n_switch_like))

# The anchored isoform fraction is the qualification that decides whether a gene's
# correlation is interpretable at all; carry it as an explicit right-hand annotation
# rather than leaving it to the caption.
g$if_lab <- sprintf("%.3f", g$max_anchored_if)

pB <- ggplot(g, aes(n_switch_like, gene_name)) +
  geom_segment(aes(x = 0, xend = n_switch_like, yend = gene_name, colour = status),
               linewidth = 0.55) +
  geom_point(aes(colour = status), size = 2) +
  geom_text(aes(x = 5.35, label = if_lab), hjust = 1, size = 2.2, colour = "grey35") +
  annotate("text", x = 5.35, y = 12.9, label = "anchored\nisoform IF",
           hjust = 1, vjust = 0.5, size = 2.1, colour = "grey35", lineheight = 0.9) +
  scale_colour_manual(values = c(Confirmed = OBS_COL,
                                 `Anchored isoform too lowly expressed` = "#E69F00",
                                 `Not confirmed` = FAIL_COL),
                      name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(0, 5.45), breaks = 0:4,
                     expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(ylim = c(0.5, 13.6), clip = "off") +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Switch-like anchored pairs in long-read", y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
fig <- (pA | pB) +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figOrthogonalConfirm", width = 7.2, height = 3.5)
cat("Done. Output in", FIG_DIR, "\n")
