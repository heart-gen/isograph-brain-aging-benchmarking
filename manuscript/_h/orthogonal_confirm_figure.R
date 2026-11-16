# Supplementary figure: the genetically anchored switch pairs behave like switches on an
# independent long-read platform.
#
# Every switch call in this paper otherwise rests on short-read quantification, so
# "your switches are quantification artifacts" is an objection no other panel answers.
# This one does, on an independent platform, lab, cohort and quantifier (Aguzzoli-Heberle
# 2024 NBT; ONT DLPFC BA9/46; Bambu; n = 12), against matched nulls drawn from IsoGraph
# switch pairs of genes outside the prioritized set.
#
# (A) observed switch-like rate vs two matched nulls, for three sets: every IsoGraph
#     switch pair scored genome-wide (the global arm, where the margin over its own null
#     is thin), all detected anchored pairs, and the abundance-qualified anchored subset.
#     The joint-matched null (abundance quintile x detection-rate, isoform-count and
#     usage-SD tertiles) is the primary comparison, drawn heavier; the abundance-only null
#     is kept beside it as the secondary comparison (review, 2026-11-16). The P values
#     printed are the joint null's. The genome-wide row moved here from
#     figConceptOverview on 2026-09-22.
# (B) per-gene confirmation for the splicing-led genes, with the genes that FAIL the
#     abundance qualification marked rather than dropped -- their anchored isoform is too
#     lowly expressed in long-read for any usage correlation to be interpretable.
#
# Abundance (of the better-expressed pair member) governs whether a usage correlation is
# estimable at all at n = 12; detection, isoform count and usage variability also move the
# chance rate of a negative correlation, which is what the joint null adds.
#
# Event layer (2026-09-29): the anchored sets now default to the signal-level
# colocalization layer (`signal_coloc/`: 55 genes whose colocalizing junction is carried by
# a tissue-matched switch-pair isoform), matching the coloc-led genetics section. The CLPP
# layer (30 genes) stays available with `--events clpp` and writes a separate file.
# Reads 06_switch_mechanism/_m/switch_orthogonal_confirm/{global_null_summary.json,
#       [signal_coloc/]{anchored_summary.json, anchored_gene_confirmation.parquet}};
#       writes figOrthogonalConfirm_{a,b}.{pdf,png} (coloc) or figOrthogonalConfirm_clpp_{a,b}
#       (CLPP), one file per panel at final print size for Fig. 6, which the manuscript's
#       content/images/assemble_main_figures.py composes.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/orthogonal_confirm_figure.R [--events signal|clpp] [--panels a|b|ab]
# Panel a reads only the summary JSONs, so `--panels a` runs without the arrow package.
suppressPackageStartupMessages({
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
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
EVENTS  <- get_arg("--events", "signal")
PANELS  <- strsplit(get_arg("--panels", "ab"), "")[[1]]
stopifnot(EVENTS %in% c("signal", "clpp"), length(PANELS) > 0, all(PANELS %in% c("a", "b")))
if ("b" %in% PANELS) suppressPackageStartupMessages(library(arrow))
EV_DIR   <- if (EVENTS == "signal") file.path(OC_DIR, "signal_coloc") else OC_DIR
FIG_NAME <- if (EVENTS == "signal") "figOrthogonalConfirm" else "figOrthogonalConfirm_clpp"
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, shared with the other real-data figures.
OBS_COL  <- "#D55E00"   # observed anchored set
NULL_COL <- "grey45"    # matched null (panel a: joint primary, abundance secondary below)
NULL_COLS <- c(`Joint-matched null` = "grey25", `Abundance-matched null` = "grey62")
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
# Anchored sets are drawn on distinct transcript pairs, the primary unit: a pair anchored from
# both sides is two orientations with one long-read outcome. Run
# `switch_orthogonal_confirm --mode distinct --events <layer>` if the summary lacks them.
if (is.null(s$distinct_pairs)) {
  stop("anchored_summary.json has no distinct_pairs; run switch_orthogonal_confirm ",
       "--mode distinct --events ", EVENTS, call. = FALSE)
}
sd <- s$distinct_pairs
g0 <- fromJSON(file.path(OC_DIR, "global_null_summary.json"))

# ---------------------------------------------------------------------------
# Panel A - observed switch-like rate against the joint- and abundance-matched nulls
# ---------------------------------------------------------------------------
# Every number is read from the summary JSONs rather than transcribed, so the annotation
# cannot drift from switch_orthogonal_confirm.py's output. The genome-wide row comes from
# the --mode global-null run, whose nulls are drawn from a different background than the
# anchored arm's; the two arms are shown on one axis because the reader needs to see that
# the anchored margin (observed minus null) is much larger than the genome-wide one, not
# that either beats the other's null. Each row carries both nulls side by side: the joint
# null (heavier, left) is the primary comparison and supplies the printed P.
NULLS <- c(matched_null_joint = "Joint-matched null", matched_null = "Abundance-matched null")
null_rows <- function(src, key, x) {
  bind_rows(lapply(names(NULLS), function(nm) {
    m <- src[[nm]][[key]]
    if (is.null(m)) stop("summary lacks ", nm, "$", key, "; rerun switch_orthogonal_confirm",
                         call. = FALSE)
    data.frame(x = x, type = NULLS[[nm]], n = m$n_focal, observed = m$observed,
               null_mean = m$null_mean, lo = m$null_q025, hi = m$null_q975,
               p = m$p_empirical_two_sided)
  }))
}
rows <- list(
  list(src = g0, key = "switch_like_rate", fmt = "Genome-wide,\nall switch pairs\n(n = %s)"),
  list(src = sd, key = "switch_like_rate", fmt = "Anchored,\nall detected\n(n = %s)"),
  list(src = sd, key = "switch_like_rate_usable_only",
       fmt = "Anchored, isoform\nwell expressed\n(n = %s)"))
nl <- bind_rows(lapply(seq_along(rows), function(i)
  null_rows(rows[[i]]$src, rows[[i]]$key, i)))
nl$type <- factor(nl$type, unname(NULLS))
# Each row's observed rate and n are the joint null's. The joint null drops focal pairs
# whose joint stratum is absent from the background (6 genome-wide pairs, none anchored),
# so its observed rate can differ from the abundance null's in the fourth decimal.
chk <- nl |> group_by(x) |> summarise(gap = diff(range(observed)))
stopifnot(all(chk$gap < 0.001))
a <- nl |> filter(type == NULLS[["matched_null_joint"]]) |>
  mutate(set = vapply(x, function(i) sprintf(rows[[i]]$fmt,
                                             format(n[x == i], big.mark = ",")), ""),
         plab = sprintf("P = %s", formatC(p, format = "g", digits = 2)),
         ytop = pmax(observed, hi) + 0.045)
OFFSET <- 0.19
nl$xpos <- nl$x + ifelse(nl$type == NULLS[["matched_null_joint"]], -OFFSET, OFFSET)
nl$hj   <- ifelse(nl$type == NULLS[["matched_null_joint"]], 1.3, -0.3)

pA <- ggplot() +
  geom_linerange(data = nl, aes(xpos, ymin = lo, ymax = hi, colour = type, linewidth = type),
                 alpha = 0.45) +
  geom_point(data = nl, aes(xpos, null_mean, colour = type, shape = type), size = 1.9) +
  geom_text(data = nl, aes(xpos, null_mean, label = sprintf("%.2f", null_mean),
                           colour = type, hjust = hj), size = 2.1, show.legend = FALSE) +
  geom_point(data = a, aes(x, observed), colour = OBS_COL, size = 2.6) +
  geom_text(data = a, aes(x, observed, label = sprintf("%.2f", observed)),
            colour = OBS_COL, size = 2.5, vjust = -0.9) +
  geom_text(data = a, aes(x, ytop, label = plab), size = 2.3, colour = "grey25",
            vjust = -0.7) +
  scale_colour_manual(values = NULL_COLS, name = NULL) +
  scale_shape_manual(values = c(`Joint-matched null` = 18, `Abundance-matched null` = 5),
                     name = NULL) +
  scale_linewidth_manual(values = c(`Joint-matched null` = 2.6,
                                    `Abundance-matched null` = 1.3), name = NULL) +
  scale_x_continuous(breaks = a$x, labels = a$set, limits = c(0.45, 3.55),
                     expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 0.86), breaks = seq(0, 0.8, 0.2)) +
  labs(x = NULL, y = "Switch-like rate in long-read") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 6.5, lineheight = 0.9),
        legend.position = "top",
        legend.text = element_text(size = 6.5),
        legend.key.width = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "pt"))

cat("  long-read rates (observed vs matched nulls):\n")
for (i in seq_len(nrow(nl))) cat(sprintf("    %-40s %-24s %.3f vs %.3f [%.3f-%.3f] (P = %.2g)\n",
    gsub("\n", " ", a$set[nl$x[i]]), nl$type[i], nl$observed[i], nl$null_mean[i],
    nl$lo[i], nl$hi[i], nl$p[i]))

# ---------------------------------------------------------------------------
# Panel B - per-gene confirmation, with the abundance-disqualified genes named
# ---------------------------------------------------------------------------
if ("b" %in% PANELS) {
# Bars count distinct switch-like pairs; gene status is unit-invariant (a gene has a
# switch-like orientation exactly when it has a switch-like distinct pair).
dp_gene <- as.data.frame(read_parquet(file.path(EV_DIR, "distinct_pair_confirmation.parquet"))) |>
  group_by(gene) |>
  summarise(n_distinct_switch_like = sum(pair_detected & switch_like), .groups = "drop")
g <- as.data.frame(read_parquet(file.path(EV_DIR, "anchored_gene_confirmation.parquet"))) |>
  left_join(dp_gene, by = "gene") |>
  mutate(status = case_when(
           confirmed_at_usable_abundance ~ "Confirmed",
           orthogonally_confirmed        ~ "Anchored isoform too lowly expressed",
           TRUE                          ~ "Not confirmed"),
         status = factor(status, c("Confirmed",
                                   "Anchored isoform too lowly expressed",
                                   "Not confirmed")),
         gene_name = reorder(gene_name, n_distinct_switch_like))

# The anchored isoform fraction is the qualification that decides whether a gene's
# correlation is interpretable at all; carry it as an explicit right-hand annotation
# rather than leaving it to the caption.
g$if_lab <- sprintf("%.3f", g$max_anchored_if)
# The gene count is read from the table, not assumed: the legacy expression-filter run
# had 12 genes and a fixed 13.6 y-limit clipped the 30-gene switching-filter result.
n_genes <- nrow(g)
max_n   <- max(4, max(g$n_distinct_switch_like, na.rm = TRUE))
x_max   <- max_n * 1.22 + 1          # room for the right-hand IF column past the longest bar
x_step  <- if (max_n > 8) 5 else 1
cat(sprintf("  per-gene panel: %d genes, %d confirmed (%d at usable abundance)\n",
            n_genes, sum(g$orthogonally_confirmed), sum(g$confirmed_at_usable_abundance)))

# Only the long-read-confirmed genes get a row, so the panel keeps a fixed per-row pitch
# and a height close to the 30-gene layout the composed Fig. 4 was built around; the
# unconfirmed genes are named in one wrapped line under the axis rather than dropped.
gb <- g |> filter(status != "Not confirmed") |> droplevels() |>
  mutate(gene_name = reorder(as.character(gene_name), n_distinct_switch_like))
not_conf <- g |> filter(status == "Not confirmed") |> arrange(as.character(gene_name)) |>
  pull(gene_name) |> as.character()
n_rows <- nrow(gb)
stopifnot(n_rows + length(not_conf) == n_genes)
nc_lab <- paste(strwrap(sprintf("Not confirmed (%d): %s", length(not_conf),
                                paste(not_conf, collapse = ", ")), width = 70),
                collapse = "\n")

pB <- ggplot(gb, aes(n_distinct_switch_like, gene_name)) +
  geom_segment(aes(x = 0, xend = n_distinct_switch_like, yend = gene_name, colour = status),
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
  labs(x = "Switch-like transcript pairs in long-read", y = NULL, caption = nc_lab) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30",
                                    lineheight = 1.05),
        plot.caption.position = "plot",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

}

# ---------------------------------------------------------------------------
# Save each panel at its final printed size
# ---------------------------------------------------------------------------
# Fig. 6 is composed by assemble_main_figures.py, so a and b are written as separate PDFs
# at the size they occupy in the 180-mm figure (fonts at print size).
# Panel a matches the column width of panel c below it; panel b's height follows a fixed
# 3.4-mm row pitch plus room for the axis, legend and not-confirmed line. No tags: the
# composed figure carries the panel letters.
MM <- 1 / 25.4
if ("a" %in% PANELS) save_fig(pA, paste0(FIG_NAME, "_a"), width = 82 * MM, height = 78 * MM)
if ("b" %in% PANELS) save_fig(pB, paste0(FIG_NAME, "_b"), width = 100 * MM,
                              height = (34 + 3.4 * n_rows) * MM)
cat("Done. Output in", FIG_DIR, "\n")
