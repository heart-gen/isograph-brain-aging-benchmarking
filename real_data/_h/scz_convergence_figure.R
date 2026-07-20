# Real-data panel (Fig 4E / S-real-8): schizophrenia-risk loci converge on age-sensitive
# co-switching modules that carry named candidate RBP trans-regulators.
#   Each row is an age-sensitive (SCZ-GWAS-enriched) co-switch module, defined out-of-cohort
#   in the aging caudate fits; x = number of independent SCZ-colocalized loci that land in it
#   (convergence), point sized by GO-invisible membership and coloured by aging cohort, with
#   the module's significant candidate RBP regulons (rbp_regulon combined scope, q<0.05)
#   printed alongside. Set-level convergence: 15/32 SCZ coloc loci fall in age-sensitive
#   modules vs 25% background (hypergeometric P=0.0058).
# Reads  real_data/_m/scz_age_projection/convergence.parquet
# Writes real_data/_m/figures/figSczConvergence.{pdf,png}
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        real_data/_h/scz_convergence_figure.R
suppressPackageStartupMessages({
  library(arrow)
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
FIG_DIR <- rel("real_data", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared house style (matches genetic_anchoring_figure.R).
theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.34, "cm"),
      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.y = element_blank(),
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
sig_star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                        ifelse(p < 0.05, "*", "")))

# Okabe-Ito: aging cohort in which the module was defined out-of-cohort.
COHORT_COL <- c(`GTEx caudate` = "#0072B2", `BrainSEQ caudate` = "#E69F00")
COHORT_LAB <- c(gtex_caudate_bg = "GTEx caudate", brainseq_caudate = "BrainSEQ caudate")

# Compact regulator label: RBP names only (drop q), cap at 4 + ellipsis.
compact_reg <- function(s) {
  if (is.na(s) || !nzchar(trimws(s))) return("")
  nm <- trimws(sub("\\(.*?\\)", "", strsplit(s, ";")[[1]]))
  nm <- nm[nzchar(nm)]
  if (length(nm) > 4) paste0(paste(nm[1:4], collapse = ", "), ", …")
  else paste(nm, collapse = ", ")
}

conv <- read_parquet(rel("real_data", "_m", "scz_age_projection", "convergence.parquet")) |>
  mutate(
    cohort = unname(COHORT_LAB[source]),
    reg    = vapply(candidate_rbp_regulators, compact_reg, character(1)),
    # module label carries its SCZ-GWAS (MAGMA) significance star
    ylab   = paste0(cohort, " ", module_id, " ", sig_star(scz_p))
  ) |>
  arrange(n_coloc_loci, scz_p)                       # bottom -> top ordering
conv$ylab <- factor(conv$ylab, levels = conv$ylab)

XMAX <- 13
pE <- ggplot(conv, aes(n_coloc_loci, ylab, colour = cohort)) +
  geom_segment(aes(x = 0, xend = n_coloc_loci, yend = ylab), linewidth = 0.6) +
  geom_point(aes(size = n_go_invisible)) +
  geom_text(aes(label = reg), hjust = 0, nudge_x = 0.35, size = 2.35,
            colour = "grey25") +
  scale_colour_manual(values = COHORT_COL, name = "Aging cohort") +
  scale_size_continuous(range = c(1.6, 4.2), breaks = c(0, 3, 6, 9),
                        name = "GO-invisible\nmembers") +
  scale_x_continuous(limits = c(0, XMAX), breaks = 0:6,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Independent schizophrenia-colocalized loci in module", y = NULL,
       caption = paste(
         "Module labels *, **, *** : SCZ MAGMA P < 0.05, 0.01, 0.001.",
         "15/32 SCZ colocalized loci fall in age-sensitive modules (~2x the",
         "background rate); hypergeometric P = 0.0058.", sep = "\n")) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.4)),
         size   = guide_legend(order = 2)) +
  theme(legend.position = "right", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30",
                                    margin = margin(6, 0, 0, 0)))

save_fig(pE, "figSczConvergence", width = 6.4, height = 2.9)
cat("done\n")
