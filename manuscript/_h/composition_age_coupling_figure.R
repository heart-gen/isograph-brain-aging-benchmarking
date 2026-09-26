# Supplementary real-data figure: how age-coupled is the estimated composition that the
# adjustment conditions on?
#
# The composition contrast (figCompositionRobustness) shows switch-unique genes collapsing
# in the two GTEx cortical regions and surviving in limbic/striatal ones. The manuscript
# reports that both ways -- composition as confounder, composition as part of aging -- but
# a reader cannot weigh the two without seeing how hard the covariates and the age term are
# competing for the same variance. This figure shows exactly that, and is deliberately
# built so it does NOT rescue the result: panel C makes plain that age-coupling alone does
# not predict which regions collapse.
#
# (A) Spearman rho of each MuSiC cell-type proportion against donor age, per analysis;
#     asterisks mark FDR < 0.05 within analysis. Cell types absent from a region's
#     reference panel are left blank rather than imputed.
# (B) the raw scatter behind the strongest coupling in each anatomical class, so the rho
#     is not read as an effect size it is not.
# (C) gene-level persistence of the switch-unique set against the region's median |rho|.
#     Cortical panels are the most age-coupled, yet ACC BA24 and DLPFC retain while BA9
#     and cortex do not -- coupling bounds the over-adjustment risk, it does not settle it.
#
# Reads 03_module_characterization/_m/composition_age_{coupling,persistence}.csv and
#       composition_age_samples.csv.gz
# (CSV, not parquet, so the figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figCompositionAgeCoupling.{pdf,png}.
# Run: bash 03_module_characterization/_h/04b.composition_age_coupling.sh
suppressPackageStartupMessages({
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
COMP    <- rel("03_module_characterization", "_m")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, matched to composition_robustness_figure.R so the same regions read the same
# colour across the two composition supplements.
CLASS_COLORS <- c(`Limbic / striatal` = "#0072B2",
                  `Cortical`          = "#D55E00",
                  `Disease (SCZD)`    = "#CC79A7")
COHORT_SHAPE <- c(BrainSEQ = 16, GTEx = 17)
# Diverging fill for rho: blue-white-red, symmetric about zero.
RHO_LOW <- "#2166AC"; RHO_MID <- "#F7F7F7"; RHO_HIGH <- "#B2182B"

PRETTY <- c(
  "SCZD (caudate)"                       = "Caudate, SCZD",
  "aging caudate"                        = "Caudate",
  "aging hippocampus"                    = "Hippocampus",
  "aging DLPFC"                          = "DLPFC",
  "GTEx amygdala"                        = "Amygdala",
  "GTEx hippocampus"                     = "Hippocampus",
  "GTEx caudate_basal_ganglia"           = "Caudate",
  "GTEx putamen_basal_ganglia"           = "Putamen",
  "GTEx nucleus_accumbens_basal_ganglia" = "N. accumbens",
  "GTEx anterior_cingulate_cortex_ba24"  = "ACC BA24",
  "GTEx frontal_cortex_ba9"              = "Frontal ctx BA9",
  "GTEx cortex"                          = "Cortex")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.36, "cm"),
      strip.text         = element_text(size = 7.8),
      strip.background   = element_blank(),
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

coup <- read.csv(file.path(COMP, "composition_age_coupling.csv"))
samp <- read.csv(gzfile(file.path(COMP, "composition_age_samples.csv.gz")))
pers <- read.csv(file.path(COMP, "composition_age_persistence.csv"))

label_of <- function(x) unname(PRETTY[x])
cohort_tag <- function(x) ifelse(x == "BrainSEQ", "BrainSEQ", "GTEx")
row_label <- function(label, cohort) sprintf("%s (%s)", label_of(label), cohort_tag(cohort))

# Order rows by class, then by how age-coupled the panel is: the reader should be able to
# see the cortical block sitting at the coupled end without hunting for it.
ord <- pers |>
  mutate(row_lab = row_label(label, ifelse(grepl("^GTEx", label), "GTEx", "BrainSEQ"))) |>
  arrange(factor(class, names(CLASS_COLORS)), median_abs_rho_age) |>
  select(label, row_lab, class, median_abs_rho_age, persistence,
         comp_unique_base, n_overlap)

# ---------------------------------------------------------------------------
# Panel A - rho(proportion, age) per analysis x cell type
# ---------------------------------------------------------------------------
hm <- coup |>
  inner_join(select(ord, label, row_lab), by = "label") |>
  mutate(row_lab = factor(row_lab, levels = ord$row_lab),
         star = ifelse(!is.na(fdr_age) & fdr_age < 0.05, "*", ""))
# Cell-type order: most consistently age-coupled on the right.
ct_ord <- hm |> group_by(cell_type) |>
  summarise(m = median(abs(rho_age), na.rm = TRUE), .groups = "drop") |>
  arrange(m)
hm$cell_type <- factor(hm$cell_type, levels = ct_ord$cell_type)

lim <- max(abs(hm$rho_age), na.rm = TRUE)
pA <- ggplot(hm, aes(cell_type, row_lab, fill = rho_age)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 2.8, colour = "grey15", vjust = 0.75) +
  scale_fill_gradient2(low = RHO_LOW, mid = RHO_MID, high = RHO_HIGH, midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey93",
                       name = expression(rho~"(proportion, age)")) +
  labs(x = NULL, y = NULL,
       caption = "* FDR < 0.05 within analysis; grey = cell type absent from the reference panel") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.caption = element_text(size = 6.4, colour = "grey40", hjust = 0),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Panel B - the raw scatter behind the strongest coupling in each class
# ---------------------------------------------------------------------------
# A rho is not an effect size. Showing the points for the single most age-coupled cell
# type per anatomical class keeps panel A from being read as a bigger claim than it is.
top_per_class <- pers |>
  group_by(class) |> slice_max(abs(top_rho_age), n = 1, with_ties = FALSE) |> ungroup()
# Short facet labels: three narrow panels cannot carry the cohort-qualified row label.
datB <- samp |>
  inner_join(select(top_per_class, label, class, top_cell_type, top_rho_age),
             by = c("label", "cell_type" = "top_cell_type")) |>
  mutate(facet = sprintf("%s (%s)\n%s, rho = %.2f",
                         label_of(label), cohort, cell_type, top_rho_age))

pB <- ggplot(datB, aes(age, proportion, colour = class)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  # Stacked, not side by side: three panels across cannot carry a readable strip label.
  facet_wrap(~ facet, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = CLASS_COLORS, guide = "none") +
  scale_x_continuous(breaks = c(25, 50, 75)) +
  labs(x = "Age (years)", y = "Estimated proportion") +
  theme_pub() +
  theme(strip.text = element_text(size = 6.9, lineheight = 1.1),
        axis.text.x = element_text(size = 6.8),
        panel.spacing.y = unit(5, "pt"))

# ---------------------------------------------------------------------------
# Panel C - persistence vs age-coupling: the honest, non-monotone panel
# ---------------------------------------------------------------------------
# Median |rho| rather than the count of significant cell types: reference panels differ in
# size between regions (6 vs 8 cell types), so a count would compare granularity.
datC <- ord |>
  filter(comp_unique_base > 0, is.finite(persistence)) |>
  mutate(cohort = ifelse(grepl("^GTEx", label), "GTEx", "BrainSEQ"),
         cohort = factor(cohort, names(COHORT_SHAPE)))

pC <- ggplot(datC, aes(median_abs_rho_age, persistence)) +
  geom_point(aes(colour = class, shape = cohort), size = 2) +
  ggrepel::geom_text_repel(aes(label = row_lab, colour = class), size = 2.1,
                           min.segment.length = 0, segment.size = 0.22,
                           segment.colour = "grey60", box.padding = 0.55,
                           point.padding = 0.3, max.overlaps = Inf, seed = 13,
                           max.time = 2, max.iter = 20000,
                           show.legend = FALSE) +
  scale_colour_manual(values = CLASS_COLORS, name = NULL) +
  scale_shape_manual(values = COHORT_SHAPE, name = NULL) +
  scale_y_continuous(limits = c(-0.16, 1.16), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(x = "Median |rho| of cell-type proportion with age",
       y = "Unadjusted switch-unique genes still\nswitch-unique after adjustment") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape = guide_legend(order = 2, nrow = 2)) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey92"))

fig <- (pA / ((pB | pC) + plot_layout(widths = c(0.78, 1)))) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figCompositionAgeCoupling", width = 7.4, height = 7.2)
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
