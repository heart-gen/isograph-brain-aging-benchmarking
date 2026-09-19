# Main real-data figure: genetic anchoring resolves disease variants to isoform switches.
# (A) SNCA cross-disease vignette - LBD and PD risk alleles both raise usage of the SAME
#     alternative-first-exon junction, mapping onto one IsoGraph switch pair (GO-invisible);
# (B) S-LDSC partitioned heritability of the aging switch layer across five traits
#     (single-annot cis/sQTL/eQTL enrichment) - splicing- vs expression-lean by trait;
# (C) the splicing-led coloc genes resolved to a switch pair (max CLPP, coloured by trait;
#     genes concordant with more than one trait share one 'multiple traits' colour);
# (D) verdict breakdown over all colocalized genes (n computed from the panel table);
# (E) Colocalizing genes do NOT concentrate in particular modules, in any trait: observed
#     concentration against a size-matched null that holds module sizes fixed.
#     RETRACTION NOTE: panel E previously showed per-module SCZ coloc counts (folded in
#     from the former standalone figSczConvergence) and was read as convergence. It was
#     retracted 2026-08-30. Counts alone cannot support that claim -- a gene can only
#     colocalize if it was coloc-TESTED, and MAGMA-anchored modules are enriched for the
#     same GWAS that decides which genes enter the test, so the tested pool is already
#     anchored-rich. Against ALL module genes SCZ/GTEx gives P = 0.0095; against the
#     tested pool the same counts give P = 0.30. This panel now shows the honest
#     all-trait test instead, which is null in all 10 (trait, cohort) cells.
# NOTE on panel A: SNCA is the WORKED EXAMPLE, not proof the method was needed -- it
# colocalizes identically in the all-genes background pool (PP4_sQTL 0.975). Its switch is
# confirmed on the exact contrast drawn here by BrainSEQ short-read junctions (minor-form
# usage 0.189 DLPFC / 0.234 caudate; 06_switch_mechanism/_m/junction_coloc_confirm). The
# ONT long-read 0.29% (max_anchored_if = 0.0029) measured a whole-transcript proxy and is an
# assay limitation; never cite it against SNCA. The set-level evidence is panels B-D plus
# figOrthogonalConfirm; the legend must say so.
# Reads 05_genetic_anchoring/_m/coloc/coloc_isoform_events_combined.parquet,
#       08_integration/_m/deep_dive/{snca_transcript_exons.tsv,deep_dive_panel.parquet},
#       05_genetic_anchoring/_m/ldsc/ldsc_partitioned.parquet,
#       05_genetic_anchoring/_m/module_coloc_convergence/global.parquet;
#       writes figGeneticAnchoring.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/genetic_anchoring_figure.R
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
DD_DIR  <- rel("08_integration", "_m", "deep_dive")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, shared with the other real-data figures.
TRAIT_COLORS <- c(AD = "#E69F00", ALS = "#0072B2", LBD = "#009E73",
                  PD = "#CC79A7", SCZ = "#D55E00")
ANNOT_COLORS <- c(cis = "#000000", sQTL = "#D55E00", eQTL = "#56B4E9")
ANCHOR_COLORS <- c(LBD = "#009E73", PD = "#CC79A7")

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

sig_star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                        ifelse(p < 0.05, "*", "")))

# ===========================================================================
# Panel A - SNCA cross-disease switch vignette
# ===========================================================================
ex <- read.delim(file.path(DD_DIR, "snca_transcript_exons.tsv"))
# transcript rows (canonical at bottom, the two switch isoforms above)
TX <- c(ENST00000336904 = "…336904\n(canonical)",
        ENST00000508895 = "…508895\n(switch)",
        ENST00000618500 = "…618500\n(switch)")
ex <- ex |> filter(transcript %in% names(TX))
ex$ty <- match(ex$transcript, names(TX))            # 1 canonical .. 3 switch
W_LO <- 89835.3; W_HI <- 89838.6                     # kb window on the switch region
ex$s_kb <- ex$start / 1000; ex$e_kb <- ex$end / 1000
exon_win <- ex |> filter(e_kb >= W_LO, s_kb <= W_HI)
intron <- exon_win |> group_by(transcript, ty) |>
  summarise(lo = min(s_kb), hi = max(e_kb), .groups = "drop")
# shared differential junction chr4:89835692-89836127 (switch transcripts only)
JUNC_LO <- 89835.692; JUNC_HI <- 89836.127
JUNC_MID <- (JUNC_LO + JUNC_HI) / 2
junc_df <- data.frame(ty = c(2, 3), x = JUNC_LO, xend = JUNC_HI)

pA <- ggplot() +
  geom_segment(data = intron, aes(x = lo, xend = hi, y = ty, yend = ty),
               colour = "grey70", linewidth = 0.4) +
  geom_segment(data = junc_df, aes(x = x, xend = xend, y = ty, yend = ty),
               colour = "#D55E00", linewidth = 1.2) +
  geom_rect(data = exon_win, aes(xmin = s_kb, xmax = e_kb, ymin = ty - 0.2, ymax = ty + 0.2),
            fill = "grey30", colour = NA) +
  # junction callout in the open intronic space
  annotate("text", x = 89837.05, y = 3.42, size = 2.5, fontface = "italic", colour = "#D55E00",
           label = "alt. first-exon junction\n(5' UTR switch, GO-invisible)") +
  annotate("segment", x = 89836.38, xend = JUNC_MID + 0.02, y = 3.36, yend = 3.04,
           arrow = arrow(length = unit(0.1, "cm")), colour = "#D55E00", linewidth = 0.4) +
  # disease anchors: both risk alleles raise junction usage
  annotate("text", x = 89838.05, y = 2.5, hjust = 0, size = 2.5, fontface = "bold",
           colour = "grey20", label = "risk alleles increase junction usage:") +
  annotate("text", x = 89838.05, y = 2.16, hjust = 0, size = 2.5,
           colour = ANCHOR_COLORS[["LBD"]], label = "LBD  rs7680557 (A)") +
  annotate("text", x = 89838.05, y = 1.86, hjust = 0, size = 2.5,
           colour = ANCHOR_COLORS[["PD"]],  label = "PD   rs1471483 (C)") +
  scale_x_reverse(limits = c(W_HI + 0.05, W_LO - 0.05),
                  breaks = c(89836, 89837, 89838)) +
  scale_y_continuous(breaks = 1:3, labels = unname(TX[names(TX)]), limits = c(0.6, 3.7)) +
  labs(x = "SNCA, chr4 position (kb)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 0.85),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

# ===========================================================================
# Panel B - S-LDSC single-annot enrichment of the aging switch layer, by trait
# ===========================================================================
ld <- as.data.frame(read_parquet(rel("05_genetic_anchoring", "_m", "ldsc", "ldsc_partitioned.parquet")))
lb <- ld |>
  filter(annotation == "aging", model %in% c("cis_only", "sqtl_only", "eqtl_only")) |>
  mutate(annot = recode(model, cis_only = "cis", sqtl_only = "sQTL", eqtl_only = "eQTL"),
         annot = factor(annot, c("cis", "sQTL", "eQTL")),
         trait = factor(toupper(trait), c("AD", "PD", "LBD", "ALS", "SCZ")),
         star = sig_star(enrichment_p))

pB <- ggplot(lb, aes(trait, enrichment, fill = annot)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = NA) +
  # stars run vertically: three dodged bars are narrower than "***" set horizontally
  geom_text(aes(label = star, y = enrichment + 0.06), angle = 90, hjust = 0, vjust = 0.75,
            position = position_dodge(width = 0.75), size = 2.6) +
  scale_fill_manual(values = ANNOT_COLORS, name = NULL) +
  # legend sits above the panel: an inset legend hid the tops of the tallest bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL,
       y = expression(atop("Partitioned " * italic(h)^2 * " enrichment",
                           "(aging switch layer)"))) +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal",
                      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))

# ===========================================================================
# Panel C - the resolved splicing-led coloc genes (max CLPP, by trait span)
# ===========================================================================
panel <- as.data.frame(read_parquet(file.path(DD_DIR, "deep_dive_panel.parquet")))
pc <- panel |>
  filter(resolved_to_switch_pair) |>
  mutate(trait_span = ifelse(grepl(",", concordant_traits), "multiple traits", concordant_traits),
         trait_span = factor(trait_span, intersect(c(names(TRAIT_COLORS), "multiple traits"),
                                                   trait_span)),
         gene = reorder(gene, max_clpp)) |>
  arrange(desc(max_clpp))

pC <- ggplot(pc, aes(max_clpp, gene)) +
  geom_segment(aes(x = 0, xend = max_clpp, yend = gene), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(colour = trait_span, shape = multi_locus), size = 2) +
  scale_colour_manual(values = c(TRAIT_COLORS, "multiple traits" = "grey35"),
                      name = "concordant trait(s)") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("single", "multi-locus"), name = NULL) +
  # square-root axis: a few high-CLPP genes would otherwise crush the rest against zero
  scale_x_sqrt(breaks = c(0, 0.05, 0.25, 0.5)) +
  labs(x = "max eCAVIAR CLPP (square-root scale)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "right", legend.box = "vertical",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Panel D - verdict breakdown over all colocalized genes
# ===========================================================================
vd <- panel |>
  mutate(v = case_when(
    grepl("splicing-led", verdict)        ~ "splicing-led\n(resolved switch)",
    grepl("not resolved", verdict)        ~ "splicing\n(unresolved)",
    grepl("expression-led", verdict)      ~ "expression-led\n(eQTL)",
    TRUE                                  ~ "other")) |>
  count(v) |>
  mutate(v = factor(v, c("splicing-led\n(resolved switch)", "splicing\n(unresolved)",
                         "expression-led\n(eQTL)")),
         frac = n / sum(n))
VD_COLORS <- c("splicing-led\n(resolved switch)" = "#D55E00",
               "splicing\n(unresolved)" = "#E69F00",
               "expression-led\n(eQTL)" = "#56B4E9")

pD <- ggplot(vd, aes(n, v, fill = v)) +
  geom_col(width = 0.7, colour = NA) +
  geom_text(aes(label = n), hjust = -0.35, size = 3) +
  scale_fill_manual(values = VD_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  labs(x = sprintf("colocalized genes (n = %d)", sum(vd$n)), y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Panel E - colocalizing genes do not concentrate in modules, in any trait
# Replaces the retracted SCZ count panel (see the retraction note in the header).
# The statistic is the sum of squared per-module hit counts; the null redraws the same
# number of hit genes from the same TESTED pool, so module sizes are held fixed by
# construction and a big module cannot look convergent just for being big. This is
# panel A of S-real-10 (figColocConvergence), which carries the full three-panel
# treatment including both denominators.
# ===========================================================================
COHORT_LAB <- c(gtex_caudate_bg = "GTEx", brainseq_caudate = "BrainSEQ")
TRAIT_ORDER <- c("AD", "PD", "LBD", "ALS", "SCZ")
OBS_COL  <- "#D55E00"
NULL_COL <- "grey45"

conv <- as.data.frame(read_parquet(rel("05_genetic_anchoring", "_m",
                                       "module_coloc_convergence", "global.parquet"))) |>
  filter(n_coloc_genes > 0) |>
  mutate(cohort = unname(COHORT_LAB[source]),
         trait  = factor(trait, TRAIT_ORDER),
         lab    = paste0(trait, " / ", cohort),
         plab   = sprintf("P = %.2f", concentration_p)) |>
  arrange(trait, cohort)
conv$lab <- factor(conv$lab, levels = rev(conv$lab))

stopifnot(nrow(conv) > 0, all(conv$concentration_p >= 0.05))

# Headroom for the P labels derived from the data, so a larger future value cannot be
# pushed off-panel (the same defect that silently dropped two points from S-real-10).
xmaxE <- max(c(conv$concentration_obs, conv$concentration_null_mean), na.rm = TRUE) * 1.42

pE <- ggplot(conv, aes(y = lab)) +
  geom_segment(aes(x = concentration_null_mean, xend = concentration_obs, yend = lab),
               linewidth = 0.4, colour = "grey70") +
  geom_point(aes(x = concentration_null_mean, shape = "Size-matched null"),
             colour = NULL_COL, size = 1.9) +
  geom_point(aes(x = concentration_obs, shape = "Observed"),
             colour = OBS_COL, size = 2.2) +
  geom_text(aes(x = xmaxE, label = plab), hjust = 1, size = 2.1, colour = "grey30") +
  scale_shape_manual(values = c(Observed = 16, `Size-matched null` = 18), name = NULL) +
  scale_x_continuous(limits = c(0, xmaxE), expand = expansion(mult = c(0.01, 0))) +
  labs(x = "Concentration of coloc genes across modules\n(sum of squared per-module counts)",
       y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.6),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Assemble: A on top (illustrative, deliberately the shortest row); B | C; D | E.
# Panel E fills the slot that once held an empty plot_spacer(); it now carries the
# convergence NULL rather than the retracted count panel.
# ===========================================================================
fig <- (pA) / (pB | pC) / (pD | pE) +
  plot_layout(heights = c(0.85, 1.15, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figGeneticAnchoring", width = 7.2, height = 8.6)
cat("Done. Output in", FIG_DIR, "\n")
