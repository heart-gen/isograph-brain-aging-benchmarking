# Supplementary real-data figure: the GO-invisible schizophrenia switch modules carry
# genuine, functionally-consequential isoform switching comparable to the genome-wide
# background -> GO-invisibility reflects GO's gene-level bias, not low module quality.
# (A) per-module driver functional-consequence fractions (CDS / coding-status / biotype /
#     UTR change) against the pooled background line. (B) switch coherence: nearly every
#     member gene carries a real anticorrelated transcript pair, with max switch strength
#     annotated. Reads real_data/brainseq/caudate_sczd/_m/go_invisible_gate.parquet +
#     go_invisible_gate_background.json, writes figGoInvisible.{pdf,png} to
#     real_data/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        real_data/_h/go_invisible_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT     <- find_root()
rel      <- function(...) file.path(ROOT, ...)
GATE_DIR <- rel("real_data", "brainseq", "caudate_sczd", "_m")
FIG_DIR  <- rel("real_data", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

MOD_COL <- "#D55E00"   # module switch evidence (vermillion, shared with the QTL figure)
BG_COL  <- "grey30"    # background reference
# qualitative module identity palette (Okabe-Ito), high contrast on white
MODULE_COLORS <- c(M026 = "#D55E00", M020 = "#E69F00", M010 = "#0072B2",
                   M023 = "#009E73")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.34, "cm"),
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

gate <- as.data.frame(read_parquet(file.path(GATE_DIR, "go_invisible_gate.parquet")))
bg   <- fromJSON(file.path(GATE_DIR, "go_invisible_gate_background.json"))

CONSEQ <- c(frac_cds_changed = "CDS change", frac_coding_status_change = "Coding-status",
            frac_biotype_switch = "Biotype switch", frac_utr_changed = "UTR change")
conseq_factor <- function(x) factor(CONSEQ[x], levels = unname(CONSEQ))
# order modules by phenotype significance (most significant first)
mod_levels <- gate |> arrange(pheno_fdr) |> pull(module)

# ---------------------------------------------------------------------------
# Panel A - functional-consequence fractions vs the pooled background
# ---------------------------------------------------------------------------
long <- gate |>
  select(module, all_of(names(CONSEQ))) |>
  pivot_longer(-module, names_to = "conseq", values_to = "frac") |>
  mutate(conseq = conseq_factor(conseq),
         module = factor(module, mod_levels))
bg_df <- data.frame(conseq = conseq_factor(names(CONSEQ)),
                    frac = unlist(bg[names(CONSEQ)]))

pA <- ggplot(long, aes(conseq, frac)) +
  geom_crossbar(data = bg_df, aes(x = conseq, y = frac, ymin = frac, ymax = frac),
                width = 0.66, linewidth = 0.5, colour = BG_COL, fatten = 0) +
  geom_point(aes(fill = module), shape = 21, colour = "grey20", stroke = 0.3,
             size = 2.6, position = position_dodge(width = 0.62)) +
  scale_fill_manual(values = MODULE_COLORS) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Fraction of driver switches") +
  annotate("text", x = 0.62, y = bg_df$frac[1] + 0.045, label = "background",
           hjust = 0, size = 2.4, colour = BG_COL) +
  theme_pub() + theme(legend.position = c(0.5, 0.12),
                      legend.direction = "horizontal",
                      axis.text.x = element_text(angle = 18, hjust = 1))

# ---------------------------------------------------------------------------
# Panel B - switch coherence: members carrying a real anticorrelated transcript pair
# ---------------------------------------------------------------------------
datB <- gate |>
  mutate(module = factor(module, mod_levels),
         frac_switch = pmin(genes_with_real_switch / n_genes, 1),
         lab = sprintf("max strength %.2f", max_switch_strength)) |>
  arrange(module)
pB <- ggplot(datB, aes(frac_switch, module)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = 0, xend = frac_switch, y = module, yend = module),
               linewidth = 0.45, colour = "grey80") +
  geom_point(colour = MOD_COL, size = 2.8) +
  geom_text(aes(x = frac_switch, label = lab), hjust = 1.08, vjust = -0.9, size = 2.3,
            colour = "grey35") +
  scale_x_continuous(limits = c(0, 1.04), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Members carrying a real isoform switch", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

fig <- (pA | pB) +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figGoInvisible", width = 7.2, height = 3.1)
cat("Done. Output in", FIG_DIR, "\n")
