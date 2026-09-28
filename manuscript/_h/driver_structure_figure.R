# Structural classes of the top driver isoform switches in age-associated IsoGraph modules.
#
# This was the closing panel of figTrustFunnel until the module-reproducibility subsection
# took that figure over. The numbers belong with the transcript-evidence results -- they
# say the drivers are bona fide isoform switches (they change coding sequence, UTR
# structure, biotype or coding status) rather than abundance artefacts -- so they keep
# their own display item instead of being dropped.
#
# The bars are mean fractions over age-significant trusted modules with the standard error
# across modules, not a test: there is no null here, and the enrichment of these classes
# against within-gene random pairs is the separate switch-consequence analysis
# (switch_consequence_figure.R). Read this as a description of the drivers, not evidence
# that the classes are over-represented.
#
# Reads 04_module_trust/_m/stability/module_trust_tables/driver_structure.csv (CSV, not
# parquet, so the figure builds with an `arrow` compiled without zstd).
# Writes manuscript/_m/figures/figDriverStructure.{pdf,png}.
# Run: bash 04_module_trust/_h/04e.module_trust_tables.sh
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
TBL     <- rel("04_module_trust", "_m", "stability", "module_trust_tables")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ISOGRAPH_COL <- "#D55E00"   # Okabe-Ito vermillion, as everywhere else in the project

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}

f <- file.path(TBL, "driver_structure.csv")
if (!file.exists(f)) {
  stop("missing ", f, "\nRun: bash 04_module_trust/_h/04e.module_trust_tables.sh",
       call. = FALSE)
}
dd <- read.csv(f) |>
  filter(is.finite(mean_frac)) |>
  mutate(switch_class = factor(switch_class, switch_class))

n_mod <- max(dd$n_modules, na.rm = TRUE)

p <- ggplot(dd, aes(switch_class, mean_frac)) +
  geom_col(fill = ISOGRAPH_COL, width = 0.68, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, mean_frac - se_frac),
                    ymax = pmin(1, mean_frac + se_frac)),
                width = 0.2, linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Mean driver fraction",
       caption = sprintf(paste("Age-significant trusted IsoGraph modules (n = %d).",
                               "\nError bars: standard error across modules."),
                         n_mod)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.caption = element_text(size = 6.3, colour = "grey40", hjust = 0))

ggsave(file.path(FIG_DIR, "figDriverStructure.pdf"), p, width = 3.8, height = 3.1,
       units = "in", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "figDriverStructure.png"), p, width = 3.8, height = 3.1,
       units = "in", dpi = 300)
cat("  figDriverStructure saved\n")
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
