# Real-data figure: genetic anchoring of co-switch modules, and where the effect lives.
#
# REFRAMED 2026-09-19 on the switching-filter re-run. Two claims the earlier version carried
# are no longer true and the figure must not imply them:
#   * the contrast does NOT concentrate in GO-invisible modules -- GO-invisible (1.080,
#     p = 0.093) and GO-visible (1.049, p = 0.126) are both null, and only the
#     phenotype-associated set clears 0.05 for IsoGraph;
#   * the effect is NOT IsoGraph-only -- the matched wgcna_multiplex baseline, fed the same
#     switch+abundance features, shows a STRONGER contrast (1.084, p = 0.0044) than IsoGraph
#     (1.065, p = 0.020), while the switch-only baseline is null and below 1.
# The honest read is a REPRESENTATION effect: both methods that consume switch+abundance
# features show splicing specificity, the switch-only representation does not, and the network
# inference is not what produces it. Panel C is built to show that directly.
#
# (A) raw cis-QTL ORs (both depleted, sQTL > eQTL) -> (B) splicing-specificity contrast by
# module set for IsoGraph (only the phenotype-associated set clears) -> (C) the same contrast
# across all three methods, full 17-analysis meta (filled) with the shared-tissue subset
# (open) beside it. Reads 05_genetic_anchoring/_m/qtl_anchoring_meta/*.parquet, writes
# figQtlSpecificity.{pdf,png} to manuscript/_m/figures/.
#
# ENCODING (must be stated in the caption): point FILL is significance -- solid = p < 0.05,
# hollow = not -- and in panel C point OPACITY is the arm (solid colour = all 17 analyses,
# pale = the shared-tissue subset). Fill and opacity therefore mean different things; a
# hollow pale point is a non-significant shared-tissue estimate. Cells with k < 3 are not
# drawn.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/qtl_specificity_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
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
MT_DIR  <- rel("05_genetic_anchoring", "_m", "qtl_anchoring_meta")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito (shared with the other real-data figures).
METHOD_COLORS <- c(isograph = "#D55E00", wgcna_switch_only = "#0072B2",
                   wgcna_multiplex = "#009E73")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna_switch_only = "WGCNA switch-only",
                   wgcna_multiplex = "WGCNA multiplex")
QTL_COLORS    <- c(sQTL = "#D55E00", eQTL = "#56B4E9")

SET_LABELS <- c(all_modules = "All", pheno_sig_modules = "Pheno-sig",
                go_invisible_modules = "GO-invis.", go_visible_modules = "GO-vis.")
set_factor <- function(x) factor(SET_LABELS[x], levels = unname(SET_LABELS))

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.38, "cm"),
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

meta     <- as.data.frame(read_parquet(file.path(MT_DIR, "qtl_anchoring_meta.parquet")))
contrast <- as.data.frame(read_parquet(file.path(MT_DIR, "qtl_anchoring_meta_contrast.parquet")))
common   <- as.data.frame(read_parquet(file.path(MT_DIR, "qtl_anchoring_meta_contrast_common.parquet")))

sig_star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                        ifelse(p < 0.05, "*", "ns")))

# ---------------------------------------------------------------------------
# Panel A - raw pooled cis-QTL ORs (IsoGraph): both depleted, sQTL > eQTL
# ---------------------------------------------------------------------------
a <- meta |>
  filter(graph_method == "isograph") |>
  transmute(set = set_factor(module_set), xqtl = factor(xqtl_kind, c("eQTL", "sQTL")),
            or = or_fe, lo = or_fe_low, hi = or_fe_high)

pA <- ggplot(a, aes(set, or, colour = xqtl)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5),
                  size = 0.35, linewidth = 0.5) +
  scale_colour_manual(values = QTL_COLORS) +
  labs(x = NULL, y = "cis-QTL odds ratio\n(module vs background)") +
  theme_pub() +  # guides are collected to a single strip at the figure foot
  theme(axis.text.x = element_text(size = 7, angle = 20, hjust = 1))

# ---------------------------------------------------------------------------
# Panel B - splicing-specificity contrast (sQTL OR / eQTL OR), IsoGraph, 17 analyses
# ---------------------------------------------------------------------------
b <- contrast |>
  filter(graph_method == "isograph") |>
  transmute(set = set_factor(module_set), ratio = ratio_fe,
            lo = ratio_fe_low, hi = ratio_fe_high, p = p_fe) |>
  mutate(lab = sprintf("%.2f %s", ratio, sig_star(p)))

# Significant sets are drawn solid, non-significant ones hollow, so a reader cannot take
# the GO-invisible point for a positive result at a glance.
b$fill <- ifelse(b$p < 0.05, METHOD_COLORS[["isograph"]], "white")

pB <- ggplot(b, aes(ratio, set)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi, fill = fill), shape = 21,
                  colour = METHOD_COLORS[["isograph"]], size = 0.4, linewidth = 0.55) +
  geom_text(aes(x = hi, label = lab), hjust = -0.18, size = 2.5) +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(x = "Splicing specificity (sQTL OR / eQTL OR)", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel C - the same contrast across all three methods.
# Filled = the full 17-analysis meta (the primary arm); open = the shared-tissue subset,
# which holds tissue constant but drops to k = 1-5 and clears 0.05 almost nowhere.
# This panel is where the matched-baseline reversal has to be visible: wgcna_multiplex,
# on identical features, sits at or above IsoGraph.
# ---------------------------------------------------------------------------
SETS_C <- c("pheno_sig_modules", "go_invisible_modules", "go_visible_modules")

c_full <- contrast |>
  filter(module_set %in% SETS_C) |>
  transmute(set = set_factor(module_set),
            method = factor(graph_method, names(METHOD_LABELS)),
            ratio = ratio_fe, lo = ratio_fe_low, hi = ratio_fe_high, p = p_fe, k = k,
            arm = "All 17 analyses")

c_common <- common |>
  filter(module_set %in% SETS_C) |>
  transmute(set = set_factor(module_set),
            method = factor(graph_method, names(METHOD_LABELS)),
            ratio = ratio_fe, lo = ratio_fe_low, hi = ratio_fe_high, p = p_fe, k = k,
            arm = "Shared tissues")

# A pooled estimate over a single analysis is not a meta-analysis; the shared-tissue
# wgcna_switch_only / GO-visible cell has k = 1 and a CI spanning 0.67-1.49, which would set
# the x-range for every other point. Drop k < 3 rather than draw it and then argue about it.
K_MIN <- 3
c_df <- bind_rows(c_full, c_common) |>
  filter(k >= K_MIN) |>
  mutate(arm = factor(arm, c("All 17 analyses", "Shared tissues")),
         fill = ifelse(p < 0.05, as.character(METHOD_COLORS[as.character(method)]), "white"))

pC <- ggplot(c_df, aes(ratio, set, colour = method, fill = fill, alpha = arm)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.75),
                  shape = 21, size = 0.3, linewidth = 0.45) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_fill_identity(guide = "none") +
  scale_alpha_manual(values = c(`All 17 analyses` = 1, `Shared tissues` = 0.45),
                     name = NULL, labels = c("All analyses", "Shared tissues")) +
  guides(colour = guide_legend(order = 1, override.aes = list(fill = METHOD_COLORS)),
         alpha  = guide_legend(order = 2, override.aes = list(colour = "grey30",
                                                              fill = "grey30"))) +
  labs(x = "Splicing specificity (sQTL OR / eQTL OR)", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Assemble: A on top, B | C below
# ---------------------------------------------------------------------------
# Panel C carries two legends and three methods, so it gets the full width rather than
# sharing a row with B.
fig <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"),
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left", legend.spacing.y = unit(1, "pt"),
        legend.margin = margin(0, 6, 0, 0, "pt"))

save_fig(fig, "figQtlSpecificity", width = 7.2, height = 5.6)

# Print what the panels assert, so a change of direction is visible in the run log rather
# than only in the rendered figure.
cat("  panel B (IsoGraph, by module set):\n")
for (i in seq_len(nrow(b))) cat(sprintf("    %-14s %.3f (%.3f-%.3f) p=%.3g%s\n",
    as.character(b$set[i]), b$ratio[i], b$lo[i], b$hi[i], b$p[i],
    ifelse(b$p[i] < 0.05, "", "   [n.s.]")))
cat("  panel C (all 17 analyses, pheno-sig):\n")
cc <- c_full[c_full$set == SET_LABELS[["pheno_sig_modules"]], ]
for (i in seq_len(nrow(cc))) cat(sprintf("    %-18s %.3f (%.3f-%.3f) p=%.4g%s\n",
    as.character(cc$method[i]), cc$ratio[i], cc$lo[i], cc$hi[i], cc$p[i],
    ifelse(cc$p[i] < 0.05, "", "   [n.s.]")))
cat("Done. Output in", FIG_DIR, "\n")
