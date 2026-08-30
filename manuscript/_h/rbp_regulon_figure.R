# Supplementary real-data figure: RBP regulons of the switch modules.
# (A) recurrent regulators - RBP motifs whose switched-exon enrichment recurs across regions;
# (B) breadth of significant module x RBP enrichments per region, split by GO-invisibility.
# Reads 07_rbp_regulation/_m/rbp/rbp_regulon.parquet; writes figRbpRegulon.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/rbp_regulon_figure.R
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(patchwork)
})
find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d); d
}
ROOT <- find_root(); rel <- function(...) file.path(ROOT, ...)
FIG_DIR <- rel("manuscript", "_m", "figures"); dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
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

rb <- as.data.frame(read_parquet(rel("07_rbp_regulation", "_m", "rbp", "rbp_regulon.parquet")))
sig <- rb |> filter(q < 0.05)

# ---- Panel A: recurrent regulators (regions with a significant enrichment) ----
rec <- sig |> group_by(rbp) |>
  summarise(n_regions = n_distinct(region),
            n_go_invisible = n_distinct(region[go_invisible]), .groups = "drop") |>
  slice_max(n_regions, n = 14, with_ties = FALSE) |>
  arrange(n_regions) |> mutate(rbp = factor(rbp, rbp))
pA <- ggplot(rec, aes(n_regions, rbp)) +
  geom_col(fill = "#009E73", width = 0.72, colour = NA) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)),
                     breaks = scales::breaks_width(2)) +
  labs(x = "Regions with significant\nswitched-exon enrichment (of 10)", y = NULL) +
  theme_pub()

# ---- Panel B: significant module x RBP enrichments per region, by GO-invisibility ----
perreg <- sig |>
  mutate(cls = ifelse(go_invisible, "GO-invisible", "GO-visible")) |>
  count(region, cls) |>
  group_by(region) |> mutate(tot = sum(n)) |> ungroup() |>
  mutate(region = reorder(region, tot),
         cls = factor(cls, c("GO-visible", "GO-invisible")))
pB <- ggplot(perreg, aes(n, region, fill = cls)) +
  geom_col(width = 0.72, colour = NA) +
  scale_fill_manual(values = c(`GO-invisible` = "#0072B2", `GO-visible` = "#E69F00"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Significant module x RBP enrichments (q < 0.05)", y = NULL) +
  theme_pub() + theme(legend.position = c(0.72, 0.18),
                      axis.text.y = element_text(size = 6.8))

fig <- pA + pB + plot_layout(widths = c(1, 1.35)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "figRbpRegulon", width = 7.2, height = 3.4)
cat("Done.\n")
