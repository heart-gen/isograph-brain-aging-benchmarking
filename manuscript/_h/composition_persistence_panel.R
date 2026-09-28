# Standalone gene-level persistence panel, sized to drop into main Fig. 2 in place of the
# count-ratio panel.
#
# This is panel B of figCompositionRobustness rendered on its own. It is a separate script
# rather than a crop of that figure because the panel there is one cell of a patchwork: it
# borrows panel A's legend and inherits the sheet's margins, so lifting it leaves a panel
# with no scale of its own.
#
# What it plots, and why it is not the panel it replaces: the fraction of the UNADJUSTED
# switch-unique genes that are still switch-unique after composition adjustment
# (n_overlap / switch_unique_base). The main figure currently plots switch_unique_adj /
# switch_unique_base under the label "Retained fraction after adjustment". That is a ratio
# of two counts, not a retained fraction -- adjustment adds genes as well as dropping them,
# so the adjusted set is not nested in the unadjusted one and the ratio can exceed 1 (GTEx
# putamen 2.00, BrainSEQ DLPFC 1.27, against true persistence of 1.00 and 0.36). Everything
# here has a subset as its numerator, so it is bounded by 1 by construction.
#
# No legend, by design: the class colours are the same Okabe-Ito values as the paired-count
# panel beside it, which carries the legend for both.
#
# Reads manuscript/_m/supp_tables/tableS13_composition_adjustment.csv (CSV, not parquet, so
# the panel builds with an `arrow` compiled without zstd; regenerate it with
# manuscript/_h/assemble_supp_tables.py).
# Writes manuscript/_m/figures/figCompositionPersistence.{pdf,png}.
#
# Default size is 3.6 x 5.2 in -- the area this panel occupied in the 7.2 x 5.2 in
# supplementary sheet, so it drops into the composite at unchanged type size.
# Run: Rscript manuscript/_h/composition_persistence_panel.R [width_in] [height_in]
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args   <- commandArgs(trailingOnly = TRUE)
WIDTH  <- if (length(args) >= 1) as.numeric(args[1]) else 3.6
HEIGHT <- if (length(args) >= 2) as.numeric(args[2]) else 5.2

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
TBL     <- rel("manuscript", "_m", "supp_tables", "tableS13_composition_adjustment.csv")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, identical to figCompositionRobustness and to the paired-count panel of Fig. 2.
CLASS_COLORS <- c(`Limbic / striatal` = "#0072B2",
                  `Cortical`          = "#D55E00",
                  `Disease (SCZD)`    = "#CC79A7")

CLASS_OF <- c(
  "SCZD (caudate)"                       = "Disease (SCZD)",
  "aging caudate"                        = "Limbic / striatal",
  "aging hippocampus"                    = "Limbic / striatal",
  "aging DLPFC"                          = "Cortical",
  "GTEx amygdala"                        = "Limbic / striatal",
  "GTEx hippocampus"                     = "Limbic / striatal",
  "GTEx caudate_basal_ganglia"           = "Limbic / striatal",
  "GTEx putamen_basal_ganglia"           = "Limbic / striatal",
  "GTEx nucleus_accumbens_basal_ganglia" = "Limbic / striatal",
  "GTEx anterior_cingulate_cortex_ba24"  = "Cortical",
  "GTEx frontal_cortex_ba9"              = "Cortical",
  "GTEx cortex"                          = "Cortical")

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

# The five GTEx regions with no defensibly matched Tran/LIBD snRNA reference were not
# deconvolved rather than forced against a mismatched panel. They are drawn as named rows:
# a reader has to be able to see this covers 12 of 17 analyses, not all of them.
NOT_DECONVOLVED <- c("Cerebellum", "Cerebellar hem.", "Hypothalamus",
                     "Spinal cord C1", "Substantia nigra")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text  = element_text(size = 7.5),
      axis.title = element_text(size = 8.5),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
      plot.margin = margin(4, 6, 4, 4, "pt")
    )
}

if (!file.exists(TBL)) {
  stop("missing ", TBL, "\nRun: python manuscript/_h/assemble_supp_tables.py", call. = FALSE)
}
dat <- read.csv(TBL, check.names = FALSE) |>
  mutate(class = factor(unname(CLASS_OF[region]), names(CLASS_COLORS)),
         label = paste0(unname(PRETTY[region]), " (", ifelse(cohort == "GTEx",
                                                             "GTEx", "BrainSEQ"), ")"))
stopifnot(all(!is.na(dat$class)))

# An analysis with zero switch-unique genes before adjustment has nothing for adjustment to
# remove: 0/0 is undefined here, not a loss, so it is dropped from this panel (the paired
# count panel keeps it as an honest zero).
pb <- dat |>
  filter(switch_unique_base > 0, is.finite(n_overlap)) |>
  mutate(persist = n_overlap / switch_unique_base,
         ann = paste0(n_overlap, "/", switch_unique_base, " kept, +", n_new)) |>
  select(label, class, persist, ann) |>
  arrange(persist)
pb_all <- bind_rows(pb, data.frame(label = paste0(NOT_DECONVOLVED, " (GTEx)"),
                                   class = NA_character_, persist = NA_real_,
                                   ann = NA_character_))
pb_all$label <- factor(pb_all$label, levels = pb_all$label)

p <- ggplot(pb_all, aes(persist, label)) +
  geom_segment(aes(x = 0, xend = persist, yend = label, colour = class),
               linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(colour = class), size = 1.8, na.rm = TRUE) +
  # Rows near 1 would push their annotation off the panel, so flip it inside the segment.
  geom_text(data = subset(pb_all, !is.na(persist)),
            aes(x = ifelse(persist > 0.6, persist - 0.03, persist + 0.03), label = ann,
                hjust = ifelse(persist > 0.6, 1, 0)),
            size = 2.1, colour = "grey35") +
  geom_text(data = subset(pb_all, is.na(persist)),
            aes(x = 0.02, label = "not deconvolved"),
            hjust = 0, size = 2.2, colour = "grey45") +
  scale_colour_manual(values = CLASS_COLORS, na.value = "grey70", guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Unadjusted switch-unique genes still\nswitch-unique after adjustment",
       y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8))

ggsave(file.path(FIG_DIR, "figCompositionPersistence.pdf"), p,
       width = WIDTH, height = HEIGHT, units = "in", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "figCompositionPersistence.png"), p,
       width = WIDTH, height = HEIGHT, units = "in", dpi = 300)
cat(sprintf("  figCompositionPersistence saved (%.2f x %.2f in)\n", WIDTH, HEIGHT))
cat("Done. Output in ", FIG_DIR, "\n", sep = "")
