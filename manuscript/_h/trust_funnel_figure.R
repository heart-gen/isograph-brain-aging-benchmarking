# Publication figure for the IsoGraph module trust funnel (real data).
# (A) the funnel in one panel: every reproducibility claim beside the matched WGCNA
#     baseline it must be read against -- trusted-module rate, split-half age sign
#     concordance, and cross-cohort eigengene projection in both directions. This panel
#     closed figConceptOverview until 2026-09-22; it opens this figure now because the
#     split-half (D) and projection results it summarises are this figure's content.
# (B) Q1 stability -> (C) Q2 driver reproducibility -> (D) Q3 within-cohort split-half
# aging concordance -> (E) Q4 structural-switch drivers. (Q3 was cross-cohort until
# 2026-09-19; see the panel D header and crosscohort_replication_figure.R.)
# Reads 04_module_trust/_m/stability/module_trust/*.parquet and
# 04_module_trust/_m/stability/eigengene_projection/eigengene_projection_summary.parquet;
# writes figTrustFunnel.{pdf,png} to manuscript/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/trust_funnel_figure.R
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
MT_DIR  <- rel("04_module_trust", "_m", "stability", "module_trust")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito: IsoGraph vermillion, WGCNA reddish purple (shared with synthetic_benchmark.R).
METHOD_COLORS <- c(isograph = "#D55E00", wgcna = "#CC79A7")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna = "WGCNA")
TRUST_COLORS  <- c(`TRUE` = "#D55E00", `FALSE` = "grey75")

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

# Concatenate every per-region file for a funnel stage (method encoded in the filename).
load_stage <- function(prefix) {
  files <- list.files(MT_DIR, pattern = paste0("^", prefix, ".*\\.parquet$"), full.names = TRUE)
  files <- files[!grepl("_pooled__", files)]   # pooled tables are empty by the n=3 floor
  bind_rows(lapply(files, function(f) {
    df <- as.data.frame(read_parquet(f))
    method <- if (grepl("__wgcna\\.parquet$", f)) "wgcna" else "isograph"
    df$method <- method
    df
  }))
}

REGION_LABELS <- c(
  caudate = "Caudate", hippocampus = "Hippocampus", dlpfc = "DLPFC",
  caudate_basal_ganglia = "Caudate (GTEx)", frontal_cortex_ba9 = "Frontal ctx (GTEx)"
)
relabel_region <- function(x) ifelse(x %in% names(REGION_LABELS), REGION_LABELS[x], x)

# ---------------------------------------------------------------------------
# Panel A - the funnel as observed vs the matched WGCNA baseline
#
# Every real-data reproducibility claim on a common fraction scale beside the thing it
# must be read against. The panel is deliberately unflattering where the data are:
# IsoGraph sits BELOW WGCNA on the first two rows and above it on the transfer rows, and
# both counts are printed so the annotation cannot imply a difference the data do not
# contain. Every value is read from a result file; nothing is typed in by hand.
# ---------------------------------------------------------------------------
stab <- load_stage("module_stability__")
wc_all <- load_stage("within_cohort__")
proj <- as.data.frame(read_parquet(rel("04_module_trust", "_m", "stability",
                                       "eigengene_projection",
                                       "eigengene_projection_summary.parquet")))

frac <- function(num, den) list(num = as.integer(num), den = as.integer(den))
trusted_frac <- function(meth) {
  d <- stab[stab$method == meth, ]
  frac(sum(d$trusted), nrow(d))
}
within_frac <- function(meth) {
  d <- wc_all[wc_all$method == meth & wc_all$both_age_sig, ]
  frac(sum(d$sign_concordant), nrow(d))
}
proj_frac <- function(meth, dir) {
  r <- proj[proj$method == meth & proj$direction == dir, ]
  stopifnot(nrow(r) == 1)
  frac(r$sign_match, r$n_age_testable)
}
claim_row <- function(claim, iso, wg) {
  data.frame(claim = claim,
             method = c("isograph", "wgcna"),
             num = c(iso$num, wg$num), den = c(iso$den, wg$den))
}
ev <- bind_rows(
  claim_row("Modules chance-trusted",
            trusted_frac("isograph"), trusted_frac("wgcna")),
  # plain "to" rather than an arrow glyph: the export font drops U+2192 in axis text
  claim_row("Split-half age sign\nconcordance",
            within_frac("isograph"), within_frac("wgcna")),
  claim_row("Aging axis transfers,\nBrainSEQ to GTEx",
            proj_frac("isograph", "brainseq_to_gtex"), proj_frac("wgcna", "brainseq_to_gtex")),
  claim_row("Aging axis transfers,\nGTEx to BrainSEQ",
            proj_frac("isograph", "gtex_to_brainseq"), proj_frac("wgcna", "gtex_to_brainseq"))
) |>
  mutate(value = num / den,
         claim = factor(claim, rev(unique(claim))),
         lab = sprintf("%d/%d", num, den),
         # IsoGraph labels above the point, WGCNA below, so both counts are always legible
         vjust = ifelse(method == "isograph", -0.9, 1.9))
ev_seg <- ev |> select(claim, method, value) |>
  pivot_wider(names_from = method, values_from = value)

pA <- ggplot(ev, aes(y = claim)) +
  geom_segment(data = ev_seg, aes(x = wgcna, xend = isograph, yend = claim),
               colour = "grey70", linewidth = 0.45) +
  geom_point(aes(x = value, colour = method), size = 2) +
  geom_text(aes(x = value, label = lab, colour = method, vjust = vjust),
            size = 2.15, show.legend = FALSE) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS, guide = "none") +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  labs(x = "Fraction of modules (split-half concordance: of both-significant pairs; projection: of age-testable modules)",
       y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7, lineheight = 0.9),
        axis.title.x = element_text(size = 7.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

cat("  funnel summary panel:\n")
for (cl in levels(ev$claim)) {
  r <- ev[ev$claim == cl, ]
  cat(sprintf("    %-58s IsoGraph %s (%.3f)  WGCNA %s (%.3f)\n", gsub("\n", " ", cl),
              r$lab[r$method == "isograph"], r$value[r$method == "isograph"],
              r$lab[r$method == "wgcna"], r$value[r$method == "wgcna"]))
}

# ---------------------------------------------------------------------------
# Panel B - Q1 stability: per-module co-assignment density vs permutation null
# ---------------------------------------------------------------------------
trust_counts <- stab |>
  group_by(method) |>
  summarise(trusted = sum(trusted), M = n(), .groups = "drop") |>
  # Report the rate, not just the count. IsoGraph 236/266 and WGCNA 64/73 are the SAME
  # trusted fraction (~89% vs ~88%); the claim this panel supports is that IsoGraph
  # modules are trustworthy in absolute terms at ~3.6x finer granularity, NOT that
  # IsoGraph beats WGCNA on trust. Showing the percentage keeps the annotation honest
  # against the paper's own not-globally-superior framing.
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
wc <- wc_all |> filter(method == "isograph", is.finite(driver_load_rho))
wc$region_lab <- factor(relabel_region(wc$region),
                        levels = relabel_region(names(REGION_LABELS)))
wc_med <- wc |> group_by(region_lab) |>
  summarise(med = median(driver_load_rho), .groups = "drop")

pC <- ggplot(wc, aes(region_lab, driver_load_rho)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(fill = METHOD_COLORS[["isograph"]], colour = NA, alpha = 0.22, scale = "width") +
  geom_boxplot(width = 0.16, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_text(data = wc_med, aes(region_lab, 1.06, label = sprintf("%.2f", med)),
            size = 2.4) +
  coord_cartesian(ylim = c(-0.6, 1.12)) +
  labs(x = NULL, y = expression("Driver-loading "*rho*" (split-half)")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------------
# Panel D - Q3 WITHIN-cohort split-half aging concordance.
#
# This panel used to be the BrainSEQ -> GTEx cross-cohort scatter. The PI ruled on
# 2026-09-19 that the cross-cohort comparison is not a fair replication test: BrainSEQ
# transcripts are Salmon-quantified and GTEx RSEM, so a failure there confounds pipeline
# with biology. It moved to the supplement (crosscohort_replication_figure.R) and the
# within-cohort split-half analogue took its place -- same question, same visual grammar,
# one cohort and one quantifier at a time.
#
# Gene-Jaccard would be the wrong statistic here: it is granularity-confounded (a method
# with a few giant modules wins by construction), which is why Q1 above uses gene-level
# co-assignment. This panel conditions on the pairs where the age effect is detectable in
# BOTH halves and asks only whether the two halves agree on its direction.
# ---------------------------------------------------------------------------
wc <- wc_all |>
  filter(is.finite(age_effect_a), is.finite(age_effect_b))
wc$method_lab <- METHOD_LABELS[wc$method]
wc$detected <- wc$both_age_sig

conc_counts <- wc |> group_by(method) |>
  summarise(both = sum(both_age_sig),
            conc = sum(both_age_sig & sign_concordant),
            M = n(), .groups = "drop") |>
  mutate(lab = sprintf("%s  %d/%d", METHOD_LABELS[method], conc, both))
conc_lab <- paste(c("Same age direction in both halves", conc_counts$lab), collapse = "\n")

wlim <- max(abs(c(wc$age_effect_a, wc$age_effect_b)), na.rm = TRUE) * 1.04

pD <- ggplot(wc, aes(age_effect_a, age_effect_b)) +
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
# Panel E - Q4 structural class of driver isoform switches (age-significant modules)
# ---------------------------------------------------------------------------
comp <- load_stage("module_complementarity__") |> filter(method == "isograph", age_sig)
SWITCH_COLS <- c(drv_cds_changed = "CDS change", drv_utr_changed = "UTR change",
                 drv_biotype_switch = "Biotype switch",
                 drv_coding_status_change = "Coding-status change")
dd <- comp |>
  select(region, all_of(names(SWITCH_COLS))) |>
  pivot_longer(-region, names_to = "switch", values_to = "frac") |>
  filter(is.finite(frac)) |>
  group_by(switch) |>
  summarise(mean_frac = mean(frac), se = sd(frac) / sqrt(n()), .groups = "drop") |>
  mutate(switch = factor(SWITCH_COLS[switch], levels = unname(SWITCH_COLS)))

pE <- ggplot(dd, aes(switch, mean_frac)) +
  geom_col(fill = METHOD_COLORS[["isograph"]], width = 0.68, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, mean_frac - se), ymax = pmin(1, mean_frac + se)),
                width = 0.2, linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Mean driver fraction") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------------
# Assemble: summary row A full width, then (B | C) / (D | E)
# ---------------------------------------------------------------------------
design <- "AAAA
           BBCC
           BBCC
           DDEE
           DDEE"
# A's row labels are wide; if its panel were axis-aligned with B's, B would be squeezed
# until its two count annotations collide. free() (patchwork >= 1.2) releases A from
# that alignment; older patchwork falls back to the aligned layout.
pA_cell <- if (exists("free", where = asNamespace("patchwork"), inherits = FALSE)) {
  patchwork::free(pA)
} else pA
fig <- wrap_plots(pA_cell, pB, pC, pD, pE, design = design) +
  plot_layout(guides = "collect", heights = c(1.5, 1, 1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 8, 0, 0, "pt"))

save_fig(fig, "figTrustFunnel", width = 7.2, height = 9.0)
cat("Done. Output in", FIG_DIR, "\n")
