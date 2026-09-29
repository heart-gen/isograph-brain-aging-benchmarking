# Supplementary figure: the genetically anchored switch pairs behave like switches on an
# independent long-read platform.
#
# Every switch call in this paper otherwise rests on short-read quantification, so
# "your switches are quantification artifacts" is an objection no other panel answers.
# This one does, on an independent platform, lab, cohort and quantifier (Aguzzoli-Heberle
# 2024 NBT; ONT DLPFC BA9/46; Bambu; n = 12), against an abundance-matched null drawn
# from 67,167 non-splicing-led IsoGraph switch pairs.
#
# (A) observed switch-like rate vs the matched null, for three sets: every IsoGraph
#     switch pair scored genome-wide (the global arm, where the margin over its own null
#     is thin), all detected anchored pairs, and the abundance-qualified anchored subset.
#     The genome-wide row moved here from figConceptOverview on 2026-09-22, so the two
#     long-read rows that panel carried now sit beside the analysis that produced them.
#     The unmatched background rate is drawn under the anchored sets so the reader can
#     see what the matching is correcting for.
# (B) per-gene confirmation for the splicing-led genes, with the genes that FAIL the
#     abundance qualification marked rather than dropped -- their anchored isoform is too
#     lowly expressed in long-read for any usage correlation to be interpretable.
#
# Matching is on the abundance decile of the better-expressed pair member, because at
# n = 12 that is what governs whether a usage correlation is estimable at all.
#
# Event layer (2026-09-29): the anchored sets now default to the signal-level
# colocalization layer (`signal_coloc/`: 55 genes whose colocalizing junction is carried by
# a tissue-matched switch-pair isoform), matching the coloc-led genetics section. The CLPP
# layer (30 genes) stays available with `--events clpp` and writes a separate file.
# Reads 06_switch_mechanism/_m/switch_orthogonal_confirm/{global_null_summary.json,
#       [signal_coloc/]{anchored_summary.json, anchored_gene_confirmation.parquet}};
#       writes figOrthogonalConfirm_{a,b}.{pdf,png} (coloc) or figOrthogonalConfirm_clpp_{a,b}
#       (CLPP), one file per panel at final print size for the hand-composed Fig. 4.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/orthogonal_confirm_figure.R [--events signal|clpp]
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
args    <- commandArgs(trailingOnly = TRUE)
EVENTS  <- if (length(args) >= 2 && args[1] == "--events") args[2] else "signal"
stopifnot(EVENTS %in% c("signal", "clpp"))
EV_DIR   <- if (EVENTS == "signal") file.path(OC_DIR, "signal_coloc") else OC_DIR
FIG_NAME <- if (EVENTS == "signal") "figOrthogonalConfirm" else "figOrthogonalConfirm_clpp"
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

s <- fromJSON(file.path(EV_DIR, "anchored_summary.json"))
g0 <- fromJSON(file.path(OC_DIR, "global_null_summary.json"))

# ---------------------------------------------------------------------------
# Panel A - observed switch-like rate against the abundance-matched null
# ---------------------------------------------------------------------------
# Every number is read from the summary JSONs rather than transcribed, so the annotation
# cannot drift from switch_orthogonal_confirm.py's output. The genome-wide row comes from
# the --mode global-null run, whose null is drawn from a different background than the
# anchored arm's; the two arms are shown on one axis because the reader needs to see that
# the anchored margin (observed minus null) is an order of magnitude larger than the
# genome-wide one, not that either beats the other's null.
mk_row <- function(src, key, label, arm) {
  m <- src$matched_null[[key]]
  data.frame(set = label, arm = arm, n = m$n_focal, observed = m$observed,
             null_mean = m$null_mean, lo = m$null_q025, hi = m$null_q975,
             p = m$p_empirical_two_sided)
}
a <- bind_rows(
  mk_row(g0, "switch_like_rate",
         sprintf("Genome-wide,\nall switch pairs\n(n = %s)",
                 format(g0$matched_null$switch_like_rate$n_focal, big.mark = ",")),
         "global"),
  mk_row(s, "switch_like_rate",
         sprintf("Anchored,\nall detected\n(n = %d)", s$matched_null$switch_like_rate$n_focal),
         "anchored"),
  mk_row(s, "switch_like_rate_usable_only",
         sprintf("Anchored, isoform\nusably expressed\n(n = %d)",
                 s$matched_null$switch_like_rate_usable_only$n_focal),
         "anchored"))
a$set <- factor(a$set, levels = a$set)
a$plab <- sprintf("P = %s", formatC(a$p, format = "g", digits = 2))
a$margin <- a$observed - a$null_mean
bg_unmatched <- s$background$switch_like_rate_unmatched
anchored_x <- range(as.integer(a$set[a$arm == "anchored"])) + c(-0.4, 0.4)

pA <- ggplot(a, aes(set)) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = NULL_COL,
                 linewidth = 2.6, alpha = 0.35) +
  # unmatched background of the ANCHORED arm: what its rate would be compared against
  # WITHOUT matching; drawn under the anchored sets only, since the global arm has its
  # own background. These numeric-x annotations must follow a discrete-x layer, or
  # ggplot trains the x scale as continuous and the discrete `set` column then fails.
  annotate("segment", x = anchored_x[1], xend = anchored_x[2],
           y = bg_unmatched, yend = bg_unmatched,
           linewidth = 0.35, linetype = "dotted", colour = "grey45") +
  annotate("text", x = anchored_x[1] + 0.05, y = bg_unmatched,
           label = "unmatched background", hjust = 0, vjust = -0.5,
           size = 2.1, colour = "grey45") +
  geom_point(aes(y = null_mean), colour = NULL_COL, size = 1.9, shape = 18) +
  geom_point(aes(y = observed), colour = OBS_COL, size = 2.6) +
  geom_text(aes(y = observed, label = sprintf("%.2f", observed)),
            colour = OBS_COL, size = 2.5, hjust = -0.45) +
  geom_text(aes(y = null_mean, label = sprintf("null %.2f", null_mean)),
            colour = NULL_COL, size = 2.1, hjust = 1.25) +
  geom_text(aes(y = pmax(hi, observed), label = plab), size = 2.3, colour = "grey25",
            vjust = -0.9) +
  scale_y_continuous(limits = c(0, 0.86), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Switch-like rate in long-read") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 6.5, lineheight = 0.9))

cat("  long-read rates (observed vs matched null):\n")
for (i in seq_len(nrow(a))) cat(sprintf("    %-40s %.3f vs %.3f (margin %+.3f, %s)\n",
    gsub("\n", " ", a$set[i]), a$observed[i], a$null_mean[i], a$margin[i], a$plab[i]))

# ---------------------------------------------------------------------------
# Panel B - per-gene confirmation, with the abundance-disqualified genes named
# ---------------------------------------------------------------------------
g <- as.data.frame(read_parquet(file.path(EV_DIR, "anchored_gene_confirmation.parquet"))) |>
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
# The gene count is read from the table, not assumed: the legacy expression-filter run
# had 12 genes and a fixed 13.6 y-limit clipped the 30-gene switching-filter result.
n_genes <- nrow(g)
max_n   <- max(4, max(g$n_switch_like, na.rm = TRUE))
x_max   <- max_n * 1.22 + 1          # room for the right-hand IF column past the longest bar
x_step  <- if (max_n > 8) 5 else 1
cat(sprintf("  per-gene panel: %d genes, %d confirmed (%d at usable abundance)\n",
            n_genes, sum(g$orthogonally_confirmed), sum(g$confirmed_at_usable_abundance)))

# Only the long-read-confirmed genes get a row, so the panel keeps a fixed per-row pitch
# and a height close to the 30-gene layout the composed Fig. 4 was built around; the
# unconfirmed genes are named in one wrapped line under the axis rather than dropped.
gb <- g |> filter(status != "Not confirmed") |> droplevels() |>
  mutate(gene_name = reorder(as.character(gene_name), n_switch_like))
not_conf <- g |> filter(status == "Not confirmed") |> arrange(as.character(gene_name)) |>
  pull(gene_name) |> as.character()
n_rows <- nrow(gb)
stopifnot(n_rows + length(not_conf) == n_genes)
nc_lab <- paste(strwrap(sprintf("Not confirmed (%d): %s", length(not_conf),
                                paste(not_conf, collapse = ", ")), width = 70),
                collapse = "\n")

pB <- ggplot(gb, aes(n_switch_like, gene_name)) +
  geom_segment(aes(x = 0, xend = n_switch_like, yend = gene_name, colour = status),
               linewidth = 0.55) +
  geom_point(aes(colour = status), size = 2) +
  geom_text(aes(x = x_max - 0.1, label = if_lab), hjust = 1, size = 2.2, colour = "grey35") +
  annotate("text", x = x_max - 0.1, y = n_rows + 0.9, label = "anchored\nisoform IF",
           hjust = 1, vjust = 0.5, size = 2.1, colour = "grey35", lineheight = 0.9) +
  scale_colour_manual(values = c(Confirmed = OBS_COL,
                                 `Anchored isoform too lowly expressed` = "#E69F00"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, x_max), breaks = seq(0, max_n, by = x_step),
                     expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(ylim = c(0.5, n_rows + 1.6), clip = "off") +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Switch-like anchored pairs in long-read", y = NULL, caption = nc_lab) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30",
                                    lineheight = 1.05),
        plot.caption.position = "plot",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Save each panel at its final printed size
# ---------------------------------------------------------------------------
# Fig. 4 is composed by hand, so a and b are written as separate PDFs at the size they
# occupy in the 180-mm figure (fonts at print size); place them at 100% with no rescaling.
# Panel a matches the column width of panel c below it; panel b's height follows a fixed
# 3.4-mm row pitch plus room for the axis, legend and not-confirmed line. No tags: the
# composed figure carries the panel letters.
MM <- 1 / 25.4
save_fig(pA, paste0(FIG_NAME, "_a"), width = 82 * MM, height = 78 * MM)
save_fig(pB, paste0(FIG_NAME, "_b"), width = 100 * MM,
         height = (34 + 3.4 * n_rows) * MM)
cat("Done. Output in", FIG_DIR, "\n")
