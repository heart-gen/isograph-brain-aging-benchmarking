# Supplementary real-data figure: clinical consequence of the switch layer.
# (A) switch genes are more loss-of-function constrained (lower gnomAD LOEUF) than genome-wide,
#     in every region; (B) switched (differentially-used) exons carry LOWER ClinVar P/LP density
#     than the gene's constitutive exons - expected alternative-exon biology, robust to a
#     coding-only (CDS) scope; (C) the colocalized splicing-led genes are themselves constrained.
# Reads 06_switch_mechanism/_m/clinical_consequence_meta.parquet, the per-region constraint_summary.parquet,
# and 05_genetic_anchoring/_m/deep_dive/deep_dive_panel.parquet; writes figClinicalConsequence.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/clinical_consequence_figure.R
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
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
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
          panel.grid.major.x = element_blank(), plot.margin = margin(4, 6, 4, 4, "pt"))
}
save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300); cat("  ", name, " saved\n", sep = "")
}

# ---- Panel A: per-region switch-gene vs genome-wide median LOEUF ----
cs <- bind_rows(lapply(
  Sys.glob(rel("02_module_discovery", "*", "*", "_m", "isograph_vae", "clinical_consequence",
               "constraint_summary.parquet")),
  function(f) as.data.frame(read_parquet(f)))) |>
  filter(stratum == "all", !is.na(median_loeuf_switch))
csl <- cs |>
  select(region, switch = median_loeuf_switch, `all genes` = median_loeuf_all) |>
  pivot_longer(-region, names_to = "set", values_to = "loeuf") |>
  mutate(set = factor(set, c("all genes", "switch")))
pA <- ggplot(csl, aes(set, loeuf)) +
  geom_line(aes(group = region), colour = "grey78", linewidth = 0.4) +
  geom_point(aes(colour = set), size = 1.7) +
  scale_colour_manual(values = c(`all genes` = "grey55", switch = "#D55E00"), guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  labs(x = NULL, y = "Median LOEUF\n(lower = more constrained)") +
  theme_pub() + theme(axis.text.x = element_text(size = 8))

# ---- Panel B: switched vs constitutive-exon ClinVar P/LP density ----
cc <- as.data.frame(read_parquet(rel("06_switch_mechanism", "_m", "clinical_consequence_meta.parquet")))
b <- cc |> filter(stratum == "all") |>
  select(scope, switched = median_switched_per_kb, constitutive = median_bg_per_kb, median_ratio) |>
  pivot_longer(c(switched, constitutive), names_to = "exon_set", values_to = "per_kb") |>
  mutate(scope = recode(scope, all_exons = "All exons", cds = "CDS only"),
         exon_set = factor(exon_set, c("switched", "constitutive")))
rat <- cc |> filter(stratum == "all") |>
  mutate(scope = recode(scope, all_exons = "All exons", cds = "CDS only"))
pB <- ggplot(b, aes(exon_set, per_kb, fill = exon_set)) +
  geom_col(width = 0.66, colour = NA) +
  geom_text(data = rat, aes(x = 1.5, y = 17.8, label = sprintf("ratio %.2f", median_ratio)),
            inherit.aes = FALSE, size = 2.6) +
  facet_wrap(~scope) +
  scale_fill_manual(values = c(switched = "#0072B2", constitutive = "grey65"), guide = "none") +
  labs(x = NULL, y = "ClinVar P/LP variants per kb") +
  theme_pub() + theme(strip.background = element_blank(),
                      strip.text = element_text(size = 8, face = "bold"),
                      axis.text.x = element_text(angle = 15, hjust = 1))

# ---- Panel C: LOEUF of colocalized genes by verdict ----
panel <- as.data.frame(read_parquet(rel("05_genetic_anchoring", "_m", "deep_dive", "deep_dive_panel.parquet")))
pc <- panel |> filter(!is.na(loeuf)) |>
  mutate(v = case_when(grepl("splicing-led", verdict) ~ "splicing-led",
                       grepl("not resolved", verdict) ~ "splicing\n(unresolved)",
                       TRUE ~ "expression-led"),
         v = factor(v, c("splicing-led", "splicing\n(unresolved)", "expression-led")))
pC <- ggplot(pc, aes(v, loeuf, colour = v)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA, linewidth = 0.5) +
  geom_jitter(width = 0.14, height = 0, size = 1.1, alpha = 0.7) +
  scale_colour_manual(values = c(`splicing-led` = "#D55E00",
                                 `splicing\n(unresolved)` = "#E69F00",
                                 `expression-led` = "#56B4E9"), guide = "none") +
  labs(x = NULL, y = "gnomAD LOEUF") +
  theme_pub() + theme(axis.text.x = element_text(size = 7))

fig <- pA + pB + pC + plot_layout(widths = c(0.9, 1.25, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "figClinicalConsequence", width = 7.4, height = 3.1)
cat("Done.\n")
