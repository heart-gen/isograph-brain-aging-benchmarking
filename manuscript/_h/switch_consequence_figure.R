# Supplementary real-data figure: coding consequence of IsoGraph isoform switches.
# (A) across 9 structural classes, switched transcript pairs are enriched for productive
#     UTR and CDS remodeling but not for decay routing (NMD), biotype change, or
#     coding-status loss; (B) the productive-remodeling signal is indistinguishable between
#     GO-invisible and GO-visible switch modules. Reads 06_switch_mechanism/_m/switch_consequence_meta.parquet;
# writes figSwitchConsequence.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/switch_consequence_figure.R
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(patchwork)
})
find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d); d
}
ROOT <- find_root(); rel <- function(...) file.path(ROOT, ...)
FIG_DIR <- rel("real_data", "_m", "figures"); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(size = 7.5), axis.title = element_text(size = 8.5),
          legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
          legend.key.size = unit(0.36, "cm"),
          panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
          panel.grid.major.y = element_blank(), plot.margin = margin(4, 6, 4, 4, "pt"))
}
save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300); cat("  ", name, " saved\n", sep = "")
}
sig_star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "")))

sc <- as.data.frame(read_parquet(rel("real_data", "_m", "switch_consequence_meta.parquet")))
LAB <- c(utr_changed = "UTR remodeled", coding_consequence = "CDS remodeled",
         cds_changed = "CDS length change", first_exon_changed = "First exon change",
         internal_exon_difference = "Internal exon diff.", last_exon_changed = "Last exon change",
         nmd_switch = "NMD routing", biotype_switch = "Biotype switch",
         coding_status_change = "Coding-status loss")

# ---- Panel A: enrichment across consequence classes (all switch genes) ----
a <- sc |> filter(stratum == "all", consequence != "cds_changed") |>
  mutate(lab = LAB[consequence],
         enriched = median_enrichment > 1 & n_enriched_p05 == n_regions,
         star = ifelse(enriched, sprintf("%d/%d %s", n_enriched_p05, n_regions, sig_star(fisher_p)), ""),
         lab = reorder(lab, median_enrichment))
pA <- ggplot(a, aes(median_enrichment, lab, fill = enriched)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_col(width = 0.68, colour = NA) +
  geom_text(aes(label = star), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey72"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Median within-gene enrichment\n(observed / permuted)", y = NULL) +
  theme_pub()

# ---- Panel B: GO-invisible vs GO-visible for the two productive classes ----
b <- sc |> filter(consequence %in% c("utr_changed", "coding_consequence"),
                  stratum %in% c("go_invisible", "go_visible")) |>
  mutate(lab = LAB[consequence],
         stratum = factor(recode(stratum, go_invisible = "GO-invisible", go_visible = "GO-visible"),
                          c("GO-invisible", "GO-visible")))
pB <- ggplot(b, aes(lab, median_enrichment, fill = stratum)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = NA) +
  scale_fill_manual(values = c(`GO-invisible` = "#0072B2", `GO-visible` = "#E69F00"), name = NULL) +
  coord_cartesian(ylim = c(0, 1.45)) +
  labs(x = NULL, y = "Median enrichment") +
  theme_pub() + theme(panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
                      legend.position = c(0.5, 0.92), legend.direction = "horizontal")

fig <- pA + pB + plot_layout(widths = c(1.4, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "figSwitchConsequence", width = 7.2, height = 3.1)
cat("Done.\n")
