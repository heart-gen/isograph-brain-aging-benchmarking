# Main figure for the module-reproducibility results subsection: are the modules
# reproducible, and does their aging structure carry to another cohort?
#
# Six panels, each attached to one sentence of that subsection, and every one of them
# drawn beside the matched WGCNA baseline it has to be read against:
#
# (A) the four reproducibility claims on one fraction scale, with both numerators and
#     denominators printed. The panel is deliberately unflattering where the data are:
#     IsoGraph sits BELOW WGCNA on the first two rows and above it on the transfer rows.
# (B) Q1 stability -- per-module split-half co-assignment density against the size-matched
#     permutation null (93/118 vs 62/72 trusted).
# (C) Q2 drivers -- shared-gene driver-loading rho across split halves, per region. This is
#     the panel that pins the regional medians the text quotes; read them off the panel
#     rather than from any hand-written summary.
# (D) Q3 within-cohort -- split-half age-effect scatter, conditioned on the pairs where the
#     effect is detectable in BOTH halves (22/23 and 14/14).
# (E) cross-cohort transfer -- per-module age z in the source cohort against the target,
#     for frozen eigengene weights. Null-standardised on purpose: raw projected age
#     correlations inherit each cohort's own age-correlated structure (see
#     EIGENGENE_PROJECTION.md), so the raw statistic measures the target cohort, not the
#     module.
# (F) what survives when the gene partition does not -- matched cross-cohort pairs against
#     a size-matched null on GO-term overlap and cell-type profile. IsoGraph clears both;
#     WGCNA's GO overlap sits on its null. This is the positive half of a paragraph whose
#     first half is a null result, and it is the reason the subsection does not end on one.
#
# The structural-driver panel that used to close this figure moved to
# driver_structure_figure.R: the structural classes belong with the transcript-evidence
# results, not with reproducibility.
#
# Reads 04_module_trust/_m/stability/module_trust_tables/*.csv (CSV, not parquet, so the
# figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figTrustFunnel.{pdf,png}.
# Run: bash 04_module_trust/_h/04e.module_trust_tables.sh
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
TBL     <- rel("04_module_trust", "_m", "stability", "module_trust_tables")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito: IsoGraph vermillion, WGCNA reddish purple (shared with synthetic_benchmark.R).
METHOD_COLORS <- c(isograph = "#D55E00", wgcna = "#CC79A7")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna = "WGCNA")
TRUST_COLORS  <- c(`TRUE` = "#D55E00", `FALSE` = "grey75")
# The published age model. The spline arm exists for both methods and is carried in
# crosscohort_permutation.csv; this figure shows the linear arm and says so.
PUBLISHED_MODEL <- "linear"

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.38, "cm"),
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

# Region labels carry the cohort: two of the six regions are hippocampus, and the paired
# BrainSEQ DLPFC / GTEx frontal cortex BA9 are the same anatomy under two names.
REGION_LABELS <- c(
  "brainseq/caudate"            = "Caudate (BrainSEQ)",
  "brainseq/dlpfc"              = "DLPFC (BrainSEQ)",
  "brainseq/hippocampus"        = "Hippocampus (BrainSEQ)",
  "gtex/caudate_basal_ganglia"  = "Caudate (GTEx)",
  "gtex/frontal_cortex_ba9"     = "Frontal ctx BA9 (GTEx)",
  "gtex/hippocampus"            = "Hippocampus (GTEx)")
region_label <- function(cohort, region) {
  key <- paste(cohort, region, sep = "/")
  ifelse(key %in% names(REGION_LABELS), unname(REGION_LABELS[key]), key)
}

read_tbl <- function(name) {
  f <- file.path(TBL, paste0(name, ".csv"))
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun: bash 04_module_trust/_h/04e.module_trust_tables.sh",
         call. = FALSE)
  }
  df <- read.csv(f)
  # pandas writes booleans as True/False, which read.csv leaves as character. Convert them
  # back, or every count in this figure silently becomes a sum over strings.
  for (col in names(df)) {
    v <- df[[col]]
    if (is.character(v) && all(v %in% c("True", "False", NA))) df[[col]] <- v == "True"
  }
  df
}

claims <- read_tbl("funnel_claims")
stab   <- read_tbl("split_half_modules")
pairs  <- read_tbl("split_half_pairs")
proj   <- read_tbl("projection_modules")
psum   <- read_tbl("projection_summary")
func   <- read_tbl("functional_preservation")

# ---------------------------------------------------------------------------
# Panel A - the four claims, observed vs the matched WGCNA baseline
# ---------------------------------------------------------------------------
# Every value comes from funnel_claims.csv, which the CLI computes from the same ledgers
# the Results text is written against; nothing is typed in here.
# Single-line labels: the panel is one row tall and two-line labels collide with their
# neighbours. "to" rather than an arrow glyph -- the export font drops U+2192 in axis text.
CLAIM_LABELS <- c(
  modules_chance_trusted                 = "Modules chance-trusted",
  split_half_age_sign_concordance        = "Split-half age sign concordance",
  aging_axis_transfers__brainseq_to_gtex = "Aging axis transfers, BrainSEQ to GTEx",
  aging_axis_transfers__gtex_to_brainseq = "Aging axis transfers, GTEx to BrainSEQ")

ev <- claims |>
  mutate(label = unname(CLAIM_LABELS[claim]),
         label = factor(label, rev(unname(CLAIM_LABELS))),
         lab = sprintf("%d/%d", num, den),
         # IsoGraph labels above the point, WGCNA below, so both counts stay legible
         vjust = ifelse(method == "isograph", -0.9, 1.9))
ev_seg <- ev |> select(label, method, value) |>
  pivot_wider(names_from = method, values_from = value)

pA <- ggplot(ev, aes(y = label)) +
  geom_segment(data = ev_seg, aes(x = wgcna, xend = isograph, yend = label),
               colour = "grey70", linewidth = 0.45) +
  geom_point(aes(x = value, colour = method), size = 2) +
  geom_text(aes(x = value, label = lab, colour = method, vjust = vjust),
            size = 2.15, show.legend = FALSE) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS, guide = "none") +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.7))) +
  labs(x = paste("Fraction of modules (split-half concordance: of both-significant pairs;",
                 "transfer: of age-testable modules)"),
       y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7, lineheight = 0.9),
        axis.title.x = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel B - Q1 stability: per-module co-assignment density vs permutation null
# ---------------------------------------------------------------------------
trust_counts <- stab |>
  group_by(method) |>
  summarise(trusted = sum(trusted), M = n(), .groups = "drop") |>
  # Report the rate, not just the count. The claim this panel supports is that IsoGraph
  # modules are trustworthy in absolute terms at a much finer granularity, NOT that
  # IsoGraph beats WGCNA on trust -- it does not, and the percentage keeps the annotation
  # honest against the paper's own not-globally-superior framing.
  mutate(lab = sprintf("%s\n%d/%d trusted\n(%.0f%%)",
                       METHOD_LABELS[method], trusted, M, 100 * trusted / M))
null_band <- stab |> group_by(method) |>
  summarise(null = median(null_mean, na.rm = TRUE), .groups = "drop")

pB <- ggplot(stab, aes(method, coassign_density)) +
  geom_violin(aes(fill = method), colour = NA, alpha = 0.18, scale = "width") +
  geom_jitter(aes(colour = trusted), width = 0.18, height = 0, size = 0.5, alpha = 0.7) +
  geom_crossbar(data = null_band, aes(x = method, y = null, ymin = null, ymax = null),
                width = 0.55, linewidth = 0.4, colour = "grey25", linetype = "dashed") +
  geom_text(data = trust_counts, aes(x = method, y = 1.03, label = lab),
            size = 2.35, vjust = 0, lineheight = 0.95) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  scale_colour_manual(values = TRUST_COLORS, labels = c(`TRUE` = "Trusted (FDR<0.05)",
                                                        `FALSE` = "Not trusted")) +
  scale_x_discrete(labels = METHOD_LABELS) +
  coord_cartesian(ylim = c(0, 1.24), clip = "off") +
  labs(x = NULL, y = "Co-assignment density") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.margin = margin(14, 6, 4, 4, "pt"))

# ---------------------------------------------------------------------------
# Panel C - Q2 driver-loading reproducibility (IsoGraph only; WGCNA has no tx drivers)
# ---------------------------------------------------------------------------
wc <- pairs |> filter(method == "isograph", is.finite(driver_load_rho)) |>
  mutate(region_lab = factor(region_label(cohort, region), unname(REGION_LABELS)))
wc_med <- wc |> group_by(region_lab) |>
  summarise(med = median(driver_load_rho), .groups = "drop")

pC <- ggplot(wc, aes(region_lab, driver_load_rho)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(fill = METHOD_COLORS[["isograph"]], colour = NA, alpha = 0.22,
              scale = "width") +
  geom_boxplot(width = 0.16, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_text(data = wc_med, aes(region_lab, 1.06, label = sprintf("%.2f", med)),
            size = 2.4) +
  coord_cartesian(ylim = c(-0.6, 1.12)) +
  labs(x = NULL, y = expression("Driver-loading "*rho*" (split-half)")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------------
# Panel D - Q3 WITHIN-cohort split-half aging concordance.
#
# This panel was the BrainSEQ -> GTEx cross-cohort scatter until the PI ruled on
# 2026-09-19 that the cross-cohort comparison is not a fair replication test: BrainSEQ
# transcripts are Salmon-quantified and GTEx RSEM, so a failure there confounds pipeline
# with biology. That comparison is a supplement (crosscohort_replication_figure.R) and the
# within-cohort split-half analogue took its place -- same question, one quantifier at a
# time.
#
# Gene-Jaccard would be the wrong statistic here: it is granularity-confounded (a method
# with a few giant modules wins by construction), which is why panel B uses gene-level
# co-assignment. This panel conditions on the pairs where the age effect is detectable in
# BOTH halves and asks only whether the halves agree on its direction.
# ---------------------------------------------------------------------------
wd <- pairs |> filter(is.finite(age_effect_a), is.finite(age_effect_b)) |>
  mutate(detected = as.logical(both_age_sig))

conc_counts <- wd |> group_by(method) |>
  summarise(both = sum(detected),
            conc = sum(detected & as.logical(sign_concordant)), .groups = "drop") |>
  mutate(lab = sprintf("%s  %d/%d", METHOD_LABELS[method], conc, both))
conc_lab <- paste(c("Same age direction in both halves", conc_counts$lab), collapse = "\n")

wlim <- max(abs(c(wd$age_effect_a, wd$age_effect_b)), na.rm = TRUE) * 1.04

pD <- ggplot(wd, aes(age_effect_a, age_effect_b)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted",
              colour = "grey60") +
  geom_point(aes(colour = method, alpha = detected, size = detected)) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.3), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.7), guide = "none") +
  annotate("text", x = -wlim, y = wlim, hjust = 0, vjust = 1, size = 2.3,
           lineheight = 0.95, label = conc_lab) +
  coord_equal(xlim = c(-wlim, wlim), ylim = c(-wlim, wlim)) +
  labs(x = "Age effect, split half A", y = "Age effect, split half B") +
  theme_pub() + theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Panel E - cross-cohort transfer of the aging axis, per module
#
# Until now this claim appeared only as two annotated points in panel A. Each point is one
# trusted module whose eigengene weights were frozen in the source cohort and applied
# unchanged to the target; the axes are that module's age correlation expressed as a z
# against same-weight random projections, so a point in the lower-left or upper-right
# quadrant agrees in direction across cohorts. The z spans roughly -59 to +52, so both
# axes are pseudo-log scaled; a linear scale would compress every point that matters into
# the origin.
# ---------------------------------------------------------------------------
DIR_LABELS <- c(brainseq_to_gtex = "Weights frozen in BrainSEQ, applied to GTEx",
                gtex_to_brainseq = "Weights frozen in GTEx, applied to BrainSEQ")
pe_dat <- proj |>
  filter(is.finite(age_z_source), is.finite(age_z_target)) |>
  mutate(dir_lab = factor(unname(DIR_LABELS[direction]), unname(DIR_LABELS)),
         agrees = as.logical(sign_match))
pe_lab <- psum |>
  mutate(dir_lab = factor(unname(DIR_LABELS[direction]), unname(DIR_LABELS)),
         lab = sprintf("%s  %d/%d", METHOD_LABELS[method], sign_match, n_age_testable)) |>
  group_by(dir_lab) |>
  summarise(lab = paste(c("Same age direction:", lab), collapse = "\n"),
            .groups = "drop")

pE <- ggplot(pe_dat, aes(age_z_source, age_z_target)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(aes(colour = method, shape = agrees), size = 1.3, alpha = 0.85,
             stroke = 0.5) +
  geom_label(data = pe_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.02,
             vjust = 1.05, size = 2.2, lineheight = 0.95, inherit.aes = FALSE,
             fill = "white", alpha = 0.75, label.size = 0, label.padding = unit(1.5, "pt")) +
  facet_wrap(~ dir_lab, ncol = 1) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "Direction agrees",
                                `FALSE` = "Direction disagrees")) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(-30, -10, -3, 0, 3, 10, 30)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(-30, -10, -3, 0, 3, 10, 30)) +
  labs(x = "Module age z, source cohort", y = "Module age z, target cohort",
       caption = paste("z against same-weight random projections;",
                       "axes pseudo-log scaled")) +
  theme_pub() +
  theme(axis.title = element_text(size = 7.5),
        plot.caption = element_text(size = 6.3, colour = "grey40", hjust = 0),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey92"))

# ---------------------------------------------------------------------------
# Panel F - what survives when the gene partition does not
#
# The cohorts do not reproducibly recover the same gene partition (median gene Jaccard
# 0.04, and the matched-pair age-concordance count does not clear its permutation null).
# These two measures ask the separate question of whether the matched pairs still share
# biology, against a null that re-pairs each module with a random target module FROM THE
# SAME GENE-COUNT DECILE -- necessary because all of these similarities grow with module
# size.
#
# `structure_r` had no computable pair and is omitted from the panel rather than plotted
# as zero: it is an absent upstream input, not a null result. Cell-type profiles exist
# only in the IsoGraph artifact tree, so WGCNA has no cell-type row for the same reason.
# Both absences are named in the caption, not hidden by the filter.
# ---------------------------------------------------------------------------
MEASURE_LABELS <- c(go_jaccard = "GO-term overlap (Jaccard)",
                    celltype_r = "Cell-type profile correlation")
pf_dat <- func |>
  filter(model == PUBLISHED_MODEL, n_finite > 0, measure %in% names(MEASURE_LABELS)) |>
  mutate(measure_lab = factor(unname(MEASURE_LABELS[measure]), unname(MEASURE_LABELS)),
         method_lab = factor(unname(METHOD_LABELS[method]), unname(METHOD_LABELS)),
         plab = sprintf("italic(p)==%s",
                        format(signif(p_emp, 2), scientific = FALSE)))

pf_missing <- expand.grid(measure_lab = levels(pf_dat$measure_lab),
                          method_lab = levels(pf_dat$method_lab),
                          stringsAsFactors = FALSE) |>
  anti_join(pf_dat, by = c("measure_lab", "method_lab")) |>
  mutate(measure_lab = factor(measure_lab, levels(pf_dat$measure_lab)),
         method_lab = factor(method_lab, levels(pf_dat$method_lab)))

# Method on the x axis and one facet per measure, rather than dodged methods within a
# measure: a measure can be missing for one method, and a dodge then silently offsets the
# null band from the point it belongs to.
pF <- ggplot(pf_dat, aes(method_lab, mean_matched, colour = method)) +
  geom_linerange(aes(ymin = null_mean - null_sd, ymax = null_mean + null_sd),
                 position = position_nudge(x = -0.17), linewidth = 3.4,
                 colour = "grey85") +
  geom_point(aes(y = null_mean), position = position_nudge(x = -0.17), size = 3,
             shape = 95, colour = "grey35") +
  geom_point(position = position_nudge(x = 0.08), size = 2.4) +
  geom_text(aes(label = plab), position = position_nudge(x = 0.08), parse = TRUE,
            vjust = -1.2, size = 2.3, show.legend = FALSE) +
  # Name the absent combinations in the panel. Leaving an empty x slot invites the reader
  # to read it as zero; it is an absent upstream input.
  geom_text(data = pf_missing, aes(x = method_lab, y = 0.05, label = "not\ncomputable"),
            inherit.aes = FALSE, size = 2.1, lineheight = 0.95, colour = "grey45") +
  facet_wrap(~ measure_lab, nrow = 1) +
  scale_colour_manual(values = METHOD_COLORS, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  coord_cartesian(ylim = c(-0.03, 0.31)) +
  labs(x = NULL, y = "Similarity of matched pairs",
       caption = paste("Grey: size-matched permutation null (mean +/- SD). Cell-type",
                       "profiles exist only in the IsoGraph artifact tree and\ntranscript",
                       "structure had no computable pair; both are absent inputs, not",
                       "null results.")) +
  theme_pub() +
  theme(axis.title.y = element_text(size = 7.5),
        plot.caption = element_text(size = 6.3, colour = "grey40", hjust = 0))

# ---------------------------------------------------------------------------
# Assemble: A full width, then (B | C), (D | F), E full width
# ---------------------------------------------------------------------------
design <- "AAAAAA
           BBBCCC
           BBBCCC
           DDDEEE
           DDDEEE
           FFFFFF"
# A's row labels are wide; if its panel were axis-aligned with B's, B would be squeezed
# until its two count annotations collide. free() (patchwork >= 1.2) releases A from that
# alignment; older patchwork falls back to the aligned layout.
pA_cell <- if (exists("free", where = asNamespace("patchwork"), inherits = FALSE)) {
  patchwork::free(pA)
} else pA
fig <- wrap_plots(pA_cell, pB, pC, pD, pE, pF, design = design) +
  plot_layout(guides = "collect", heights = c(1.25, 1, 1, 1, 1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 8, 0, 0, "pt"))

save_fig(fig, "figTrustFunnel", width = 7.2, height = 9.6)

cat("  funnel summary panel:\n")
for (cl in rev(levels(ev$label))) {
  r <- ev[ev$label == cl, ]
  cat(sprintf("    %-40s IsoGraph %s (%.3f)  WGCNA %s (%.3f)\n", gsub("\n", " ", cl),
              r$lab[r$method == "isograph"], r$value[r$method == "isograph"],
              r$lab[r$method == "wgcna"], r$value[r$method == "wgcna"]))
}
cat("  driver-loading rho, regional medians: ",
    paste(sprintf("%s %.2f", wc_med$region_lab, wc_med$med), collapse = "; "), "\n",
    sep = "")
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
