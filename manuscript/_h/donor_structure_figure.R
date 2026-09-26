# Supplementary real-data figure: 16 regions, 17 analyses, 844 donors -- not 17 cohorts.
#
# The manuscript states that samples from different regions can come from the same donor
# and that the SCZD caudate analysis shares controls with the BrainSEQ aging caudate
# analysis. Until now that was a sentence a reader had to take on trust, while the results
# text repeatedly reports agreement "across regions". This figure shows the sharing, so a
# reviewer can see exactly how much between-region agreement is re-measurement of one
# donor pool and how much is independent.
#
# The three identifier levels are kept apart deliberately. A donor (BrNum / SUBJID) spans
# regions; a tissue sample (sample_id, the RNum for BrainSEQ) does not; and an analysis is
# not a region -- BrainSEQ caudate is ONE region carrying TWO analyses, the aging analysis
# over its 238 controls and the schizophrenia analysis over those same 238 control samples
# plus 152 patients. The figure is therefore drawn per region (16 columns), with the
# caudate double-use called out rather than drawn as a second region.
#
# (A) donor x region incidence: one column per region, one row per donor, filled where that
#     donor contributed tissue. Donors are ordered by the set of regions they appear in, so
#     the GTEx block structure reads as a block rather than as noise.
# (B) how many regions each donor contributed, by cohort.
# (C) pairwise shared-donor counts between regions (lower triangle) -- the number the
#     "replicates across regions" claims should be read against.
#
# Reads 02_module_discovery/_m/donor_structure/{region_overlap,analysis_donor_counts,
# region_donor_counts,donor_analysis_counts,donor_incidence}.csv (CSV, not parquet, so the
# figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figDonorStructure.{pdf,png}.
# Run: bash 02_module_discovery/_h/03a.donor_structure.sh
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
DS      <- rel("02_module_discovery", "_m", "donor_structure")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, matched to the other real-data supplements.
COHORT_COL <- c(BrainSEQ = "#E69F00", GTEx = "#0072B2")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.36, "cm"),
      strip.text         = element_text(size = 8),
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

ov       <- read.csv(file.path(DS, "region_overlap.csv"))
an       <- read.csv(file.path(DS, "analysis_donor_counts.csv"))
per_d    <- read.csv(file.path(DS, "donor_analysis_counts.csv"))
inc      <- read.csv(file.path(DS, "donor_incidence.csv"))

# Region order: the order the manuscript introduces them, i.e. the analysis-table order.
REGIONS <- unique(an$region)
inc$region <- factor(inc$region, levels = REGIONS)

# Regions carrying more than one analysis get a marker, so the caudate double-use is
# visible on the axis instead of being implied by the caption alone.
multi <- an |> count(region, name = "n_analyses") |> filter(n_analyses > 1)
region_lab <- function(x) ifelse(x %in% multi$region, paste0(x, " *"), as.character(x))

# ---------------------------------------------------------------------------
# Panel A - donor x region incidence
# ---------------------------------------------------------------------------
# Order donors by their membership pattern so the block structure reads as blocks. The key
# is the sorted set of regions a donor contributed to, which groups identical patterns.
pattern <- inc |>
  group_by(cohort, donor_id) |>
  summarise(key = paste(sort(as.integer(region)), collapse = ","),
            n = n(), .groups = "drop") |>
  arrange(cohort, desc(n), key) |>
  mutate(row = row_number())
incA <- inc |> inner_join(select(pattern, donor_id, cohort, row),
                          by = c("donor_id", "cohort"))

pA <- ggplot(incA, aes(region, row, fill = cohort)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = COHORT_COL, name = NULL) +
  scale_x_discrete(labels = region_lab) +
  scale_y_reverse(expand = expansion(mult = 0.01)) +
  labs(x = NULL, y = "Donors (ordered by membership pattern)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.6),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "top")

# ---------------------------------------------------------------------------
# Panel B - regions per donor
# ---------------------------------------------------------------------------
pB <- ggplot(per_d, aes(factor(n_regions), fill = cohort)) +
  geom_bar(width = 0.8) +
  facet_wrap(~ cohort, scales = "free", ncol = 1) +
  scale_fill_manual(values = COHORT_COL, guide = "none") +
  labs(x = "Regions the donor contributed", y = "Donors") +
  theme_pub()

# ---------------------------------------------------------------------------
# Panel C - pairwise shared donors between regions
# ---------------------------------------------------------------------------
# Lower triangle only: the matrix is symmetric, and the diagonal is each region's own donor
# count, which panel A already carries.
idx <- setNames(seq_along(REGIONS), REGIONS)
datC <- ov |>
  filter(region_a != region_b, idx[region_a] < idx[region_b]) |>
  mutate(region_a = factor(region_a, levels = REGIONS),
         region_b = factor(region_b, levels = rev(REGIONS)))

pC <- ggplot(datC, aes(region_a, region_b, fill = n_shared_donors)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#F7F7F7", high = "#0072B2", name = "Shared\ndonors") +
  scale_x_discrete(labels = region_lab) +
  scale_y_discrete(labels = region_lab) +
  labs(x = NULL, y = NULL,
       caption = paste("Blank = no shared donors; BrainSEQ and GTEx donor ids are never",
                       "matched across cohorts.\n* region carrying two analyses:",
                       "BrainSEQ caudate is analysed for aging (238 controls) and for",
                       "schizophrenia\n(those same 238 control samples plus 152 patients).")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.4),
        axis.text.y = element_text(size = 6.4),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.caption = element_text(size = 6.2, colour = "grey40", hjust = 0))

fig <- ((pA | pB) + plot_layout(widths = c(1.35, 1))) / pC +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figDonorStructure", width = 7.4, height = 8)
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
