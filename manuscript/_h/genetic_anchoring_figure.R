# Main real-data figure: genetic anchoring resolves disease variants to isoform switches.
# (A) SNCA cross-disease vignette - LBD and PD risk alleles both raise usage of the SAME
#     alternative-first-exon junction, mapping onto one IsoGraph switch pair (GO-invisible);
# (B) S-LDSC partitioned heritability of the aging switch layer across five traits
#     (single-annot cis/sQTL/eQTL enrichment) - splicing- vs expression-lean by trait;
# (C) the 12 splicing-led resolved coloc genes (max CLPP, coloured by trait span; all GO-invisible);
# (D) verdict breakdown over all 68 colocalized genes.
# Reads 05_genetic_anchoring/_m/coloc/coloc_isoform_events_combined.parquet,
#       05_genetic_anchoring/_m/deep_dive/{snca_transcript_exons.tsv,deep_dive_panel.parquet},
#       05_genetic_anchoring/_m/ldsc/ldsc_partitioned.parquet; writes figGeneticAnchoring.{pdf,png}.
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
DD_DIR  <- rel("real_data", "_m", "deep_dive")
FIG_DIR <- rel("real_data", "_m", "figures")
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
  annotate("segment", x = 89836.55, xend = JUNC_MID + 0.02, y = 3.3, yend = 3.02,
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
  labs(x = "SNCA · chr4 position (kb)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 0.85),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

# ===========================================================================
# Panel B - S-LDSC single-annot enrichment of the aging switch layer, by trait
# ===========================================================================
ld <- as.data.frame(read_parquet(rel("real_data", "ldsc", "_m", "ldsc_partitioned.parquet")))
lb <- ld |>
  filter(annotation == "aging", model %in% c("cis_only", "sqtl_only", "eqtl_only")) |>
  mutate(annot = recode(model, cis_only = "cis", sqtl_only = "sQTL", eqtl_only = "eQTL"),
         annot = factor(annot, c("cis", "sQTL", "eQTL")),
         trait = factor(toupper(trait), c("AD", "PD", "LBD", "ALS", "SCZ")),
         star = sig_star(enrichment_p))

pB <- ggplot(lb, aes(trait, enrichment, fill = annot)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = NA) +
  geom_text(aes(label = star, y = enrichment + 0.15),
            position = position_dodge(width = 0.75), size = 3, vjust = 0.6) +
  scale_fill_manual(values = ANNOT_COLORS, name = NULL) +
  labs(x = NULL, y = "Partitioned h² enrichment\n(aging switch layer)") +
  theme_pub() + theme(legend.position = c(0.5, 0.92), legend.direction = "horizontal")

# ===========================================================================
# Panel C - the 12 resolved splicing-led coloc genes (max CLPP, by trait span)
# ===========================================================================
panel <- as.data.frame(read_parquet(file.path(DD_DIR, "deep_dive_panel.parquet")))
pc <- panel |>
  filter(resolved_to_switch_pair) |>
  mutate(trait_span = ifelse(grepl(",", concordant_traits), "cross-disease", concordant_traits),
         trait_col = ifelse(trait_span == "cross-disease", "ALS", concordant_traits),
         gene = reorder(gene, max_clpp)) |>
  arrange(desc(max_clpp))
# label multi-locus genes
pc$face <- ifelse(pc$multi_locus, "bold", "plain")

pC <- ggplot(pc, aes(max_clpp, gene)) +
  geom_segment(aes(x = 0, xend = max_clpp, yend = gene), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(colour = concordant_traits, shape = multi_locus), size = 2) +
  scale_colour_manual(values = c(TRAIT_COLORS, "ALS,SCZ" = "#0072B2", "LBD,PD" = "#009E73"),
                      name = "concordant trait(s)") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("single", "multi-locus"), name = NULL) +
  labs(x = "max eCAVIAR CLPP", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7, face = pc$face[order(pc$max_clpp)]),
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
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  scale_fill_manual(values = VD_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "colocalized genes (n = 68)", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ===========================================================================
# Assemble: A wide on top; B | C on the middle; D bottom-left with legend space
# ===========================================================================
fig <- (pA) / (pB | pC) / (pD | plot_spacer()) +
  plot_layout(heights = c(1.0, 1.15, 0.7)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figGeneticAnchoring", width = 7.2, height = 8.4)
cat("Done. Output in", FIG_DIR, "\n")
