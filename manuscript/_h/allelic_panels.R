# Shared panel for the within-donor allelic test of isoform choice (BrainSEQ recount).
# Sourced by allelic_imbalance_figure.R (Fig S-allelic) and genetic_anchoring_figure.R
# (Fig 5a), so the main-figure panel and its supplementary detail cannot drift apart.
#
# allelic_calibration_panel(tab): tab is tableS22_allelic_imbalance_regions.csv. Per region,
# the fraction of tests at nominal P < 0.05 for the lead-heterozygous test (gate family, the
# BH-tested set) against the same model fitted on donors HOMOZYGOUS at the lead, which has no
# allelic contrast to find and so should sit at 0.05. The label carries the BH q < 0.05 counts
# the text quotes. The two fractions have their own denominators (fitted gate-family pairs;
# fitted homozygous-null pairs) and the legend must say so.

ALLELIC_REGION_LABELS <- c(caudate = "Caudate", hippocampus = "Hippocampus", dlpfc = "DLPFC")
ALLELIC_HET_COL  <- "#D55E00"
ALLELIC_NULL_COL <- "grey45"

allelic_region_rows <- function(tab) {
  tab <- tab[tab$region %in% names(ALLELIC_REGION_LABELS), ]
  tab$region_lab <- factor(unname(ALLELIC_REGION_LABELS[tab$region]),
                           levels = rev(unname(ALLELIC_REGION_LABELS)))
  tab
}

allelic_calibration_panel <- function(tab, theme_fn) {
  tab <- allelic_region_rows(tab)
  arms <- c("Lead-heterozygous donors", "Homozygous-at-lead null")
  pts <- rbind(
    data.frame(region_lab = tab$region_lab, arm = arms[1],
               frac = tab$pairs_nominal_p05 / tab$pairs_fitted),
    data.frame(region_lab = tab$region_lab, arm = arms[2], frac = tab$hom_null_frac_p05)
  )
  pts$arm <- factor(pts$arm, arms)
  seg <- data.frame(region_lab = tab$region_lab, lo = tab$hom_null_frac_p05,
                    hi = tab$pairs_nominal_p05 / tab$pairs_fitted)
  lab <- data.frame(region_lab = tab$region_lab,
                    label = sprintf("%d pairs, %d genes", tab$pairs_q05, tab$genes_q05))
  x_lab <- max(pts$frac) + 0.02
  ggplot() +
    geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_segment(data = seg, aes(x = lo, xend = hi, y = region_lab, yend = region_lab),
                 colour = "grey75", linewidth = 0.5) +
    geom_point(data = pts, aes(frac, region_lab, colour = arm, shape = arm), size = 2.3) +
    geom_text(data = lab, aes(x_lab, region_lab, label = label), hjust = 0, size = 2.5,
              colour = "grey15") +
    annotate("text", x = x_lab, y = 3.45, label = "q < 0.05", hjust = 0, size = 2.5,
             fontface = "italic", colour = "grey30") +
    scale_colour_manual(values = setNames(c(ALLELIC_HET_COL, ALLELIC_NULL_COL), arms),
                        name = NULL) +
    scale_shape_manual(values = setNames(c(16, 18), arms), name = NULL) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.75))) +
    coord_cartesian(clip = "off") +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    labs(x = "Within-donor allelic tests at P < 0.05", y = NULL) +
    theme_fn() +
    theme(legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
          plot.margin = margin(4, 70, 4, 4, "pt"))
}
