# Supplementary real-data figure: RBP regulons of the switch modules.
# (a) per region, module x RBP cells significant by the hypergeometric test, split by
#     whether they survive the opportunity-adjusted GLM; (b) RBPs whose supported cells
#     recur across regions; (c) ENCODE eCLIP binding capacity, switched vs constitutive.
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

# ---- Panels A/B redrawn 2026-09-29 on the opportunity-adjusted arm ----
# A cell is "supported" when it is significant by the hypergeometric test (q < 0.05) AND
# survives the opportunity-adjusted binomial GLM (qvalue_glm < 0.05, OR > 1). The raw
# test alone is confounded by transcript length, GC and UTR composition, so it is drawn
# only as the denominator. Every count is read from the parquet, including the region total.
n_regions <- n_distinct(rb$region)
rb <- rb |> mutate(supported = q < 0.05 & !is.na(qvalue_glm) & qvalue_glm < 0.05 &
                     odds_ratio > 1)
sig <- rb |> filter(q < 0.05)
sup <- rb |> filter(supported)
stopifnot(nrow(rb) == 11040, nrow(sig) == 384, nrow(sup) == 89)
cat(sprintf("  regulon cells: %d tested in %d regions; %d hypergeometric q<0.05 (%d RBPs); %d survive adjustment (%d RBPs, %d regulons)\n",
            nrow(rb), n_regions, nrow(sig), n_distinct(sig$rbp), nrow(sup),
            n_distinct(sup$rbp), nrow(distinct(sup, region, module_id))))

# ---- Panel A: per region, raw hits and those that survive adjustment ----
PRETTY_REGION <- c(
  frontal_cortex_ba9 = "Frontal ctx BA9", anterior_cingulate_cortex_ba24 = "ACC BA24",
  cerebellar_hemisphere = "Cerebellar hem.", caudate_basal_ganglia = "Caudate (GTEx)",
  cortex = "Cortex", hypothalamus = "Hypothalamus", amygdala = "Amygdala",
  dlpfc = "DLPFC", caudate = "Caudate", hippocampus = "Hippocampus",
  cerebellum = "Cerebellum", caudate_sczd = "Caudate, SCZD",
  putamen_basal_ganglia = "Putamen", substantia_nigra = "Substantia nigra",
  nucleus_accumbens_basal_ganglia = "N. accumbens",
  spinal_cord_cervical_c_1 = "Spinal cord C1")
relabel <- function(x) ifelse(x %in% names(PRETTY_REGION), PRETTY_REGION[x], x)
CLS <- c("Survives opportunity adjustment", "Hypergeometric only")
perreg <- sig |>
  mutate(region = unname(relabel(region)),
         cls = factor(ifelse(supported, CLS[1], CLS[2]), rev(CLS))) |>
  count(region, cls) |>
  group_by(region) |> mutate(tot = sum(n)) |> ungroup() |>
  mutate(region = reorder(region, tot))
pA <- ggplot(perreg, aes(n, region, fill = cls)) +
  geom_col(width = 0.72, colour = NA) +
  scale_fill_manual(values = setNames(c("#D55E00", "grey80"), CLS), breaks = CLS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Module x RBP cells, hypergeometric q < 0.05", y = NULL) +
  theme_pub() + theme(legend.position = "bottom",
                      legend.margin = margin(0, 0, 0, 0),
                      axis.text.y = element_text(size = 6.8))

# ---- Panel B: RBPs whose supported cells recur across regions ----
rec <- sup |> group_by(rbp) |>
  summarise(n_regions = n_distinct(region), .groups = "drop") |>
  slice_max(n_regions, n = 14, with_ties = FALSE) |>
  arrange(n_regions) |> mutate(rbp = factor(rbp, rbp))
pB <- ggplot(rec, aes(n_regions, rbp)) +
  geom_col(fill = "#D55E00", width = 0.72, colour = NA) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)), breaks = scales::breaks_width(1)) +
  labs(x = sprintf("Regions with a supported cell (of %d)", n_regions), y = NULL) +
  theme_pub()

# ---- Panel C: ENCODE eCLIP binding CAPACITY at switched vs constitutive exons ----
# The motif panels above are sequence predictions. This panel adds the only measured
# binding layer available, and its framing is deliberately narrow.
#
# What it shows: for each nominated RBP, the fraction of its regulon genes with an eCLIP
# peak over the SWITCHED exon interval vs over that gene's CONSTITUTIVE exons. The
# within-gene contrast is what makes it interpretable -- peak-dense factors blanket the
# transcriptome, so an unpaired "supported" rate saturates at ~100% and means nothing.
#
# What it does NOT show: factor-specific occupancy of these regulons in brain. ENCODE
# eCLIP is HepG2/K562, the median switched-constitutive gap is only 0.024 (GC-matched
# arm; the deprecated flat-background arm gave 0.013), and the
# neuronal-CLIP program that would have shown brain occupancy is an honest null. This is
# binding CAPACITY at alternative-exon sequence, and the legend must say so.
#
# One two-sided exact McNemar per RBP over its unique nominated regulon genes (deduped
# across regions), BH across the testable RBPs (35 in this run) -- NOT one test per region x module,
# which would count the same region-invariant binding fact up to five times.
bs <- as.data.frame(read_parquet(
  rel("07_rbp_regulation", "_m", "rbp", "rbp_binding_support.parquet")))

n_pref <- sum(bs$preferential, na.rm = TRUE)
n_supp <- sum(bs$binding_supported, na.rm = TRUE)
cat(sprintf("  eCLIP: %d/%d preferential, %d/%d binding-supported, median gap %.3f\n",
            n_pref, nrow(bs), n_supp, nrow(bs), median(bs$rate_diff, na.rm = TRUE)))

bc <- bs |>
  slice_max(rate_diff, n = 14, with_ties = FALSE) |>
  mutate(supported = ifelse(binding_supported, "q < 0.05", "not significant"),
         supported = factor(supported, c("q < 0.05", "not significant")),
         rbp = factor(rbp, rev(rbp)))

pC <- ggplot(bc, aes(rate_diff, rbp, colour = supported)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = median(bs$rate_diff, na.rm = TRUE),
             linewidth = 0.3, linetype = "dotted", colour = "grey45") +
  geom_segment(aes(x = 0, xend = rate_diff, yend = rbp), linewidth = 0.5) +
  geom_point(size = 1.8) +
  annotate("text", x = median(bs$rate_diff, na.rm = TRUE), y = 14.7,
           label = "median", hjust = -0.12, vjust = 0.5, size = 2.1, colour = "grey45") +
  scale_colour_manual(values = c(`q < 0.05` = "#D55E00",
                                 `not significant` = "#999999"), name = NULL) +
  scale_x_continuous(limits = c(0, max(bc$rate_diff) * 1.08), breaks = scales::breaks_width(0.05)) +
  coord_cartesian(clip = "off") +
  labs(x = "eCLIP binding rate,\nswitched - constitutive exons", y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 6.8),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

fig <- (pA | (pB / pC)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "figRbpRegulon", width = 7.2, height = 6.4)
cat("Done.\n")
