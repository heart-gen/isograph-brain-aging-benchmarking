# Fig 5: cis-genetic support for isoform choice, and its links to disease loci.
# Panels follow the order the Results cite them (reorganised 2026-09-29):
# (a) within-donor allelic test (BrainSEQ recount): nominal P < 0.05 rate for the
#     lead-heterozygous test against the lead-homozygous null, with the q < 0.05 counts;
#     shared with the supplementary allelic figure through allelic_panels.R;
# (b) the 119 sQTL colocalization nominations by estimator tier (coloc.susie signal level /
#     coloc.abf fallback only), with those whose colocalizing intron lands on a
#     tissue-matched switch-pair isoform filled;
# (c) the 20 signal-level nominations on a switch-pair isoform by max PP4, coloured by trait,
#     filled where the CLPP layer resolves the same gene x trait to a switch pair.
#     (b, c) were CLPP-based until 2026-09-29; colocalization now leads, following the
#     estimator hierarchy susie > abf > CLPP, and CLPP sits in Table S24 as corroboration;
# (d) PRDM2 vignette - the ALS risk allele lowers usage of the junction that joins PRDM2's
#     distal terminal exon, shifting the gene toward an early-terminating form;
# (e) S-LDSC partitioned heritability of the aging switch layer across five traits
#     (single-annotation cis/sQTL/eQTL enrichment); stars are the one-sided coefficient
#     test (tau > 0 given baselineLD), the statistic the text quotes.
# The module-convergence null that was panel E moved to the supplement, where
# figColocConvergence carries it in full. RETRACTION NOTE kept for the record: that slot
# once showed per-module SCZ coloc counts and was read as convergence (retracted
# 2026-08-30) -- a gene can only colocalize if it was coloc-TESTED, and MAGMA-anchored
# modules are enriched for the GWAS that decides which genes enter the test.
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
#     usage is 0.045 -- the early-3'-end form is a real but minority form.
# WHICH early-terminating transcript to draw (resolved 2026-09-20; see reports/pi/08_integration
# section 6b). ENST00000413440 (PRDM2-206) and ENST00000343137 (PRDM2-203) are the SAME FORM,
# not competing candidates. Both start at the internal promoter (~13,749.4 kb), both terminate
# ~28 kb before the colocalizing junction's donor, so NEITHER can carry it and both sit on the
# alternative arm; and their CDS is IDENTICAL (3 coding exons, chr1:13,773,170-13,786,524).
# They differ only by a 56-bp cassette exon at chr1:13,769,118-13,769,173 that PRDM2-203
# includes and PRDM2-206 skips, and by 422 bp of 3' UTR (13,788,079 vs 13,787,657).
# That difference is also why the two assays 'disagree': against MANE, PRDM2-203 contributes 3
# transcript-specific introns to PRDM2-206's 1 (the cassette supplies two), so the short-read
# recount has more unique evidence for it (median minor-form usage 0.108 vs 0.045 over 498
# donors), while ONT assigns it essentially nothing (mean isoform fraction 2.8e-09 vs 0.192)
# and every long-read pair containing it is undetected. Each assay favours the form its own
# resolution can distinguish; they are not disagreeing measurements of one quantity.
# So this panel draws the alternative arm as THE EARLY-TERMINATING FORM with PRDM2-206 as its
# representative -- the only one ONT can attribute, the coordinate `_REFERENCE_ARM` already
# uses for this gene (13,788,079), and the form in the single ONT-confirmed switch-like pair
# against MANE. The legend must say that PRDM2-203 is the same form carrying an additional
# 56-bp 5' UTR cassette. Do NOT redraw this as a choice between two partners, and do not
# describe the anchored event as PRDM2's RIZ1/RIZ2 promoter choice: it is a 3' terminus choice
# nested inside that program and independent of it (the distal junction is carried by both
# RIZ1-type transcripts, CDS from 13,715,606, and by two RIZ2-type ones, CDS from 13,773,170).
# The LIBD PSI arm (`junction_coloc_confirm`) returns `junction_not_measured` for PRDM2 --
# its event catalogue stops at ~13,787,067 and never reaches the distal terminal exon. That
# is a gap in THAT catalogue only, and must not be reported as a failure to confirm. The set-level evidence is panels B-D plus figOrthogonalConfirm.
# Do NOT annotate a structural consequence here from `structural_consequence`: that field
# flags each transcript against a gene reference, so 63 of 76 concordant events carry all
# seven flags and it cannot distinguish one event from another. The structure drawn in this
# panel comes from the GENCODE exon coordinates directly.
# Reads 08_integration/_m/deep_dive/prdm2_transcript_exons.tsv,
#       05_genetic_anchoring/_m/coloc_signal_susie/all_introns/coloc_isoform_events.parquet
#       and manuscript/_m/supp_tables/{tableS22_allelic_imbalance_regions,
#       tableS24_clpp_isoform_events,
#       tableS27_ldsc_partitioned}.csv (LDSC through its table CSV: the ledger parquet is
#       zstd, which an `arrow` built without it cannot read);
#       writes figGeneticAnchoring.{pdf,png}.
# Run: python manuscript/_h/assemble_supp_tables.py && \
#      Rscript manuscript/_h/genetic_anchoring_figure.R
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
SUSIE_DIR <- rel("05_genetic_anchoring", "_m", "coloc_signal_susie", "all_introns")
TAB_DIR <- rel("manuscript", "_m", "supp_tables")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
source(rel("manuscript", "_h", "allelic_panels.R"))

need_tab <- function(name) {
  f <- file.path(TAB_DIR, name)
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun: python manuscript/_h/assemble_supp_tables.py", call. = FALSE)
  }
  read.csv(f)
}

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
# The colocalizing event is an ALTERNATIVE 3' ACCEPTOR at a shared donor, chr1:13,816,570.
# The junction drawn here, chr1:13816570-13823159, is the last intron of the forms that reach
# PRDM2's distal terminal exon; the ALS risk allele A at rs2744682 DECREASES its usage
# (risk_qtl_effect -1.10) and raises an UNANNOTATED proximal acceptor at 13,821,649 (no GENCODE
# transcript has an exon there), i.e. it shifts PRDM2 toward the early-terminating form.
# Transcripts carrying the drawn junction (MANE …311066 and the internal-promoter form …505823)
# end at 13,823,159+; …413440 lacks it and terminates early at 13,788,079.
ex <- read.delim(file.path(DD_DIR, "prdm2_transcript_exons.tsv"))
TX <- c(ENST00000413440 = "PRDM2-206 (…413440)\nearly 3' terminus",
        ENST00000505823 = "PRDM2-214 (…505823)\ndistal 3' end",
        ENST00000311066 = "PRDM2-202 (…311066)\nMANE Select, distal 3' end")
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
  annotate("text", x = W_HI + 0.5, y = 3.46, hjust = 1, size = 2.5, fontface = "italic",
           colour = "#D55E00",
           label = "colocalizing junction\n(distal terminal exon)") +
  annotate("segment", x = JUNC_MID, xend = JUNC_MID, y = 3.30, yend = 3.06,
           arrow = arrow(length = unit(0.1, "cm")), colour = "#D55E00", linewidth = 0.4) +
  annotate("text", x = W_LO + 1.0, y = 2.58, hjust = 0, size = 2.5, fontface = "bold",
           colour = "grey20", label = "ALS risk allele rs2744682 (A)") +
  annotate("text", x = W_LO + 1.0, y = 2.30, hjust = 0, size = 2.5,
           colour = ANCHOR_COLORS[["ALS"]],
           label = "decreases junction usage \u2192 early 3' end") +
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
ld <- need_tab("tableS27_ldsc_partitioned.csv")
lb <- ld |>
  filter(annotation == "aging", model %in% c("cis_only", "sqtl_only", "eqtl_only")) |>
  mutate(annot = recode(model, cis_only = "cis", sqtl_only = "sQTL", eqtl_only = "eQTL"),
         annot = factor(annot, c("cis", "sQTL", "eQTL")),
         trait = factor(toupper(trait), c("AD", "PD", "LBD", "ALS", "SCZ")),
         star = sig_star(coef_p))

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
# Coloc-layer events (reorganised 2026-09-29: colocalization leads, CLPP corroborates).
# One row per (nomination x tissue) from the all-introns signal-level arm; a nomination is
# signal-level when coloc.susie scored it in at least one tissue, and coloc.abf-only
# otherwise (the fallback tier the Methods keep separate).
# ===========================================================================
ev <- as.data.frame(read_parquet(file.path(SUSIE_DIR, "coloc_isoform_events.parquet"))) |>
  mutate(ev_switch = junction_in_switch_pair %in% TRUE, trait = toupper(trait))
# summarise() sees earlier summaries under their new names, so the per-event flag keeps a
# name of its own.
nom <- ev |>
  group_by(gene_name, trait) |>
  summarise(signal = any(estimator == "susie"),
            on_switch_signal = any(ev_switch & estimator == "susie"),
            pp4_switch_signal = suppressWarnings(max(PP4_sQTL[ev_switch & estimator == "susie"])),
            on_switch = any(ev_switch),
            .groups = "drop")
# The numbers the Results quote; fail loudly if the ledger moves under the text.
# 21 signal-level nominations reach a switch pair, 20 through a coloc.susie event itself;
# GPM6A (SCZ) is signal-level in caudate/ACC but switch-mapped only in a coloc.abf tissue.
stopifnot(nrow(nom) == 119, sum(nom$signal) == 33, sum(nom$on_switch) == 56,
          sum(nom$signal & nom$on_switch) == 21, sum(nom$on_switch_signal) == 20)

# CLPP corroboration: the same gene x trait also resolves to a switch pair in the CLPP layer.
clpp <- need_tab("tableS24_clpp_isoform_events.csv") |>
  filter(as.character(junction_in_switch_pair) %in% c("True", "TRUE")) |>
  distinct(gene_name, trait = toupper(trait))
clpp_key <- paste(clpp$gene_name, clpp$trait)

# ===========================================================================
# Panel C - the signal-level nominations that land on a switch-pair isoform (max PP4)
# ===========================================================================
# TMEM175 and SNCA head this panel. Earlier vignette decisions (see the header) rejected both
# on the representative-intron abf and CLPP layers; in the all-introns coloc.susie arm drawn
# here both are signal-level and switch-mapped (TMEM175: cortex and cerebellum; SNCA:
# hypothalamus and cerebellum, not cortex). The legend carries their SMR caveats (TMEM175
# 6/36 supported; SNCA LBD 12 supported vs 11 HEIDI-rejected). PRDM2 is abf-only and so is
# absent here by construction.
pc <- nom |>
  filter(on_switch_signal) |>
  mutate(dup = duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
         label = ifelse(dup, sprintf("%s (%s)", gene_name, trait), gene_name),
         label = reorder(label, pp4_switch_signal),
         trait = factor(trait, intersect(names(TRAIT_COLORS), trait)),
         clpp_too = factor(ifelse(paste(gene_name, trait) %in% clpp_key,
                                  "also CLPP-resolved", "coloc only"),
                           c("also CLPP-resolved", "coloc only")))

pC <- ggplot(pc, aes(pp4_switch_signal, label)) +
  geom_segment(aes(x = 0.8, xend = pp4_switch_signal, yend = label), colour = "grey80",
               linewidth = 0.4) +
  geom_point(aes(colour = trait, shape = clpp_too), size = 2, stroke = 0.7) +
  scale_colour_manual(values = TRAIT_COLORS, name = "trait") +
  scale_shape_manual(values = c("also CLPP-resolved" = 16, "coloc only" = 1), name = NULL) +
  scale_x_continuous(limits = c(0.8, 1), breaks = c(0.8, 0.9, 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = expression("Max sQTL PP"[4] * " (coloc.susie)"), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        # inside, over the empty low-PP4 corner: a side legend starved panel d of width
        # lower right: the bottom rows sit at PP4 0.80-0.90, so that corner is empty
        legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
        legend.box = "vertical", legend.box.just = "right",
        legend.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Panel B (drawn as pD) - the 119 sQTL colocalization nominations by estimator tier, with
# those whose colocalizing intron lands on a tissue-matched switch-pair isoform filled.
# ===========================================================================
TIER <- c("coloc.susie\n(signal level)", "coloc.abf\nfallback only")
RES <- c("On a switch-pair isoform", "Not on a switch pair")
vd <- nom |>
  mutate(support = factor(ifelse(signal, TIER[1], TIER[2]), rev(TIER)),
         res = factor(ifelse(on_switch, RES[1], RES[2]), rev(RES))) |>
  count(support, res, .drop = FALSE)
tot <- vd |> group_by(support) |>
  summarise(n_res = sum(n[res == RES[1]]), n = sum(n), .groups = "drop") |>
  mutate(label = sprintf("%d (%d on switch pair)", n, n_res))

pD <- ggplot(vd, aes(n, support, fill = res)) +
  geom_col(width = 0.68, colour = NA) +
  geom_text(data = tot, aes(n, support, label = label), inherit.aes = FALSE,
            hjust = -0.08, size = 2.6, colour = "grey15") +
  scale_fill_manual(values = setNames(c("#D55E00", "grey80"), RES), breaks = RES,
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.75))) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = sprintf("sQTL coloc nominations (n = %d)", sum(vd$n)), y = NULL) +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Panel A (drawn as pAllelic) - within-donor allelic test, shared with the supplement
# ===========================================================================
pAllelic <- allelic_calibration_panel(need_tab("tableS22_allelic_imbalance_regions.csv"),
                                      theme_pub)

# ===========================================================================
# Assemble in citation order: a allelic | b QTL support; c CLPP genes | (d PRDM2 / e LDSC).
# ===========================================================================
fig <- (pAllelic | pD) / ((pC | (pA / pB)) + plot_layout(widths = c(0.8, 1.2))) +
  plot_layout(heights = c(0.62, 1.38)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figGeneticAnchoring", width = 7.2, height = 8.2)
cat("Done. Output in", FIG_DIR, "\n")
