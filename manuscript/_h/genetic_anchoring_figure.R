# Main real-data figure: genetic anchoring resolves disease variants to isoform switches.
# (A) PRDM2 vignette - the ALS risk allele lowers usage of the junction that joins PRDM2's
#     distal terminal exon, shifting the gene toward an early-terminating form; the junction
#     maps onto the tissue-matched IsoGraph switch pair in cortex;
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
# NOTE on panel A: PRDM2 replaced SNCA on 2026-09-19. SNCA's LBD junction is NOT in the
# tissue-matched switch pair under the switching filter, so it is scored non-concordant and
# is kept in the manuscript as the falsification example (PI, 2026-09-18) -- never redraw it
# here as resolved.
# PRDM2 is the WORKED EXAMPLE, not proof the method was needed. It is the only one of the 30
# concordant genes corroborated by every genetic layer at once: coloc.abf PP4_sQTL 0.964 vs
# PP4_eQTL 0.494, colocalizing in 6 tissues on splicing and 0 on expression; SMR 7/7 probes
# supported with no HEIDI rejection and no eQTL instrument at all; and the only concordant
# event whose switch replicates in an independent BrainSEQ cohort (DLPFC). Its honest costs
# belong in the legend: CLPP is only 0.056, the module is GO-VISIBLE (unlike the legacy SNCA
# panel), and the junction/switch polarity r is 0.11 -- the junction rides the switch axis
# weakly.
# Orthogonal support for this exact pair, from the two assays that do cover it:
#   ONT long-read (switch_orthogonal_confirm / longread_switch_confirm, DLPFC BA9 n = 12):
#     ENST00000311066 and ENST00000413440 are BOTH detected and co-expressed, usage
#     Spearman -0.448, `switch_like = True`, `coding_consequence = True`; PRDM2 is
#     `confirmed = True` at gene level in frontal_cortex_ba9 and cortex.
#   BrainSEQ allele-aware junction recount (ase_junction_switch/dlpfc): the junction drawn
#     here is present as an isoform-1-specific junction of this very pair, testable over
#     292 donors / 4,095 fragments, with 495 of 498 donors carrying BOTH forms. Minor-form
#     usage is 0.045 -- the early-3'-end form is a real but minority form. Against the other
#     early-terminating transcript (ENST00000343137) minor-form usage is 0.108, but ONT does
#     not detect that one, so the pair drawn here is the better-supported switch.
# The LIBD PSI arm (`junction_coloc_confirm`) returns `junction_not_measured` for PRDM2 --
# its event catalogue stops at ~13,787,067 and never reaches the distal terminal exon. That
# is a gap in THAT catalogue only, and must not be reported as a failure to confirm. The set-level evidence is panels B-D plus figOrthogonalConfirm.
# Do NOT annotate a structural consequence here from `structural_consequence`: that field
# flags each transcript against a gene reference, so 63 of 76 concordant events carry all
# seven flags and it cannot distinguish one event from another. The structure drawn in this
# panel comes from the GENCODE exon coordinates directly.
# Reads 05_genetic_anchoring/_m/coloc/coloc_isoform_events_combined.parquet,
#       08_integration/_m/deep_dive/{prdm2_transcript_exons.tsv,deep_dive_panel.parquet},
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
ANCHOR_COLORS <- c(LBD = "#009E73", PD = "#CC79A7", ALS = "#0072B2")

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
# Panel A - PRDM2 switch vignette (ALS, cortex)
# ===========================================================================
# The colocalizing junction chr1:13816570-13823159 is the last intron of the forms that
# reach PRDM2's distal terminal exon. Transcripts carrying it (MANE …311066 and the
# internal-promoter form …505823) end at 13,823,159+; …413440 lacks it and terminates early
# at 13,788,079. The ALS risk allele A at rs2744682 DECREASES usage of the junction
# (risk_qtl_effect -1.10), i.e. shifts PRDM2 toward the early-terminating form.
ex <- read.delim(file.path(DD_DIR, "prdm2_transcript_exons.tsv"))
TX <- c(ENST00000413440 = "…413440\n(early 3' end)",
        ENST00000505823 = "…505823\n(switch)",
        ENST00000311066 = "…311066\n(MANE Select)")
ex <- ex |> filter(transcript %in% names(TX))
ex$ty <- match(ex$transcript, names(TX))            # 1 early-3'end .. 3 MANE
W_LO <- 13772.0; W_HI <- 13826.5                     # kb window on the 3' half of the gene
ex$s_kb <- ex$start / 1000; ex$e_kb <- ex$end / 1000
exon_win <- ex |> filter(e_kb >= W_LO, s_kb <= W_HI)
intron <- exon_win |> group_by(transcript, ty) |>
  summarise(lo = min(s_kb), hi = max(e_kb), .groups = "drop")
# colocalizing junction, drawn only on the transcripts that actually splice it
JUNC_LO <- 13816.570; JUNC_HI <- 13823.159
JUNC_MID <- (JUNC_LO + JUNC_HI) / 2
junc_df <- data.frame(ty = c(2, 3), x = JUNC_LO, xend = JUNC_HI)

pA <- ggplot() +
  geom_segment(data = intron, aes(x = lo, xend = hi, y = ty, yend = ty),
               colour = "grey70", linewidth = 0.4) +
  geom_segment(data = junc_df, aes(x = x, xend = xend, y = ty, yend = ty),
               colour = "#D55E00", linewidth = 1.2) +
  geom_rect(data = exon_win, aes(xmin = s_kb, xmax = e_kb, ymin = ty - 0.2, ymax = ty + 0.2),
            fill = "grey30", colour = NA) +
  annotate("text", x = JUNC_MID, y = 3.46, size = 2.5, fontface = "italic", colour = "#D55E00",
           label = "colocalizing junction\n(distal terminal exon)") +
  annotate("segment", x = JUNC_MID, xend = JUNC_MID, y = 3.30, yend = 3.06,
           arrow = arrow(length = unit(0.1, "cm")), colour = "#D55E00", linewidth = 0.4) +
  annotate("text", x = W_LO + 1.0, y = 2.58, hjust = 0, size = 2.5, fontface = "bold",
           colour = "grey20", label = "ALS risk allele rs2744682 (A)") +
  annotate("text", x = W_LO + 1.0, y = 2.30, hjust = 0, size = 2.5,
           colour = ANCHOR_COLORS[["ALS"]],
           label = "decreases usage of the junction \u2192 shifts to the early 3' end") +
  scale_x_continuous(limits = c(W_LO - 0.5, W_HI + 0.5),
                     breaks = c(13780, 13800, 13820)) +
  scale_y_continuous(breaks = 1:3, labels = unname(TX[names(TX)]), limits = c(0.6, 3.75)) +
  labs(x = "PRDM2, chr1 position (kb)", y = NULL) +
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
