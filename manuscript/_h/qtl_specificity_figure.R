# Main real-data figure: genetic anchoring of IsoGraph co-switch modules.
# (A) raw cis-QTL ORs (both depleted, sQTL > eQTL) -> (B) splicing-specificity contrast
# by module set (concentrates in pheno-sig + GO-invisible) -> (C) matched-baseline method
# effect on shared tissues (only IsoGraph positive). Reads
# 05_genetic_anchoring/_m/qtl_anchoring_meta/*.parquet, writes figQtlSpecificity.{pdf,png} to
# manuscript/_m/figures/.
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
MT_DIR  <- rel("real_data", "_m", "qtl_anchoring_meta")
FIG_DIR <- rel("real_data", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito (shared with the other real-data figures).
METHOD_COLORS <- c(isograph = "#D55E00", wgcna_switch_only = "#0072B2",
                   wgcna_multiplex = "#009E73")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna_switch_only = "WGCNA switch-only",
                   wgcna_multiplex = "WGCNA multiplex")
QTL_COLORS    <- c(sQTL = "#D55E00", eQTL = "#56B4E9")

SET_LABELS <- c(all_modules = "All", pheno_sig_modules = "Phenotype-sig",
                go_invisible_modules = "GO-invisible", go_visible_modules = "GO-visible")
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
  theme_pub() + theme(legend.position = c(0.82, 0.18),
                      axis.text.x = element_text(angle = 25, hjust = 1))

# ---------------------------------------------------------------------------
# Panel B - splicing-specificity contrast (sQTL OR / eQTL OR), IsoGraph, 17 analyses
# ---------------------------------------------------------------------------
b <- contrast |>
  filter(graph_method == "isograph") |>
  transmute(set = set_factor(module_set), ratio = ratio_fe,
            lo = ratio_fe_low, hi = ratio_fe_high, p = p_fe) |>
  mutate(lab = sprintf("%.2f %s", ratio, sig_star(p)))

pB <- ggplot(b, aes(ratio, set)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = METHOD_COLORS[["isograph"]],
                  size = 0.4, linewidth = 0.55) +
  geom_text(aes(x = hi, label = lab), hjust = -0.15, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = "Splicing specificity (sQTL OR / eQTL OR)", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Panel C - matched-baseline method effect on the shared tissues
# ---------------------------------------------------------------------------
c_df <- common |>
  filter(module_set %in% c("pheno_sig_modules", "go_invisible_modules", "go_visible_modules")) |>
  transmute(set = set_factor(module_set),
            method = factor(graph_method, names(METHOD_LABELS)),
            ratio = ratio_fe, lo = ratio_fe_low, hi = ratio_fe_high)

pC <- ggplot(c_df, aes(ratio, set, colour = method)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6),
                  size = 0.32, linewidth = 0.5) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  labs(x = "Splicing specificity, shared tissues", y = NULL) +
  theme_pub() + theme(legend.position = "bottom",
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# ---------------------------------------------------------------------------
# Assemble: A on top, B | C below
# ---------------------------------------------------------------------------
fig <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figQtlSpecificity", width = 7.2, height = 6.0)
cat("Done. Output in", FIG_DIR, "\n")
