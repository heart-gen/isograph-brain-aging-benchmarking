# Figure 1: what an isoform switch is, what IsoGraph does with it, and that it works
# where the truth is known.
#
# The previous Fig 1 was the 24-panel synthetic benchmark grid. That is the right content
# for a methods supplement and the wrong opener for a general-audience genomics journal:
# nothing anywhere in the repository drew what an isoform switch IS, so a reader met the
# paper's central object for the first time as a boxplot axis label.
#
# This rebuild is schematic-first. The exhaustive parameter-resolved panels are unchanged
# and stay in the supplement (S1-S12), built in place by
# isograph_benchmark/figures/synthetic_benchmark.R; two summary panels are retained here.
#
# (A) DEFINITION. One gene, two isoforms, usage crossing over with age while total
#     gene-level abundance stays flat -- the event a DGE pipeline cannot see.
# (B) WHAT IS MISSED. The DGE x DTU quadrant map: abundance pipelines cover the top row,
#     IsoGraph's contribution is the DTU-without-DGE quadrant.
# (C) METHOD. Per-gene switch features -> VAE latent -> gene-gene graph -> co-switch
#     modules.
# (D) Module recovery on synthetic ground truth, six core scenarios, six main methods.
# (E) Switch-gene detection on the same runs.
#
# Panels A-C are drawn from explicit coordinates -- they are illustrations of a
# definition, and carry no data. Panels D-E are read from the committed benchmark ledger.
#
# Reads 01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet.
# Writes manuscript/_m/figures/figConceptOverview.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/concept_overview_figure.R
suppressPackageStartupMessages({
  library(jsonlite); library(tibble)
  library(arrow)
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
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, shared with every other figure in the paper.
ISO_A   <- "#0072B2"   # isoform A (declining)
ISO_B   <- "#D55E00"   # isoform B (rising)
TOTAL   <- "#555555"   # gene-level total
HILITE  <- "#D55E00"   # the quadrant / module IsoGraph adds
MUTED   <- "grey78"

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.36, "cm"),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}
theme_blank <- function() {
  theme_void(base_size = 8.5) + theme(plot.margin = margin(2, 4, 2, 4, "pt"))
}

save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300)
  cat("  ", name, " saved\n", sep = "")
}

# ===========================================================================
# Panel A - the definition: usage crosses over, total abundance does not
# ===========================================================================
age <- seq(20, 80, length.out = 200)
# A logistic crossover in usage; the total is deliberately flat so the panel makes the
# single point it exists to make.
u_b <- 1 / (1 + exp(-(age - 50) / 7))
defn <- bind_rows(
  data.frame(age = age, value = 1 - u_b, series = "Isoform A usage"),
  data.frame(age = age, value = u_b,     series = "Isoform B usage"))
tot <- data.frame(age = age, value = rep(0.5, length(age)),
                  series = "Total gene abundance (scaled)")

pA <- ggplot() +
  geom_line(data = defn, aes(age, value, colour = series), linewidth = 0.85) +
  geom_line(data = tot, aes(age, value, colour = series), linewidth = 0.7,
            linetype = "22") +
  annotate("segment", x = 50, xend = 50, y = 0, yend = 1.0,
           linewidth = 0.3, colour = "grey70") +
  annotate("text", x = 50, y = 1.06, label = "switch", size = 2.5, colour = "grey35") +
  scale_colour_manual(
    values = c(`Isoform A usage` = ISO_A, `Isoform B usage` = ISO_B,
               `Total gene abundance (scaled)` = TOTAL),
    breaks = c("Isoform A usage", "Isoform B usage", "Total gene abundance (scaled)"),
    name = NULL) +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linetype = c("solid", "solid", "22")))) +
  labs(x = "Age (years)", y = "Fraction of gene output") +
  theme_pub() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6.8), legend.key.width = unit(0.5, "cm"))

# ===========================================================================
# Panel B - what abundance pipelines cover, and what is left
# ===========================================================================
# Structural device is a 2x2 because the content genuinely is a 2x2: the two axes are
# independent, and the paper's claim is about exactly one of the four cells.
quad <- data.frame(
  x = c(1, 2, 1, 2), y = c(2, 2, 1, 1),
  dge = c("no", "yes", "no", "yes"),
  dtu = c("no", "no", "yes", "yes"),
  lab = c("no change",
          "abundance shift\n(DGE pipelines)",
          "isoform switch,\nflat abundance",
          "both"),
  hit = c(FALSE, FALSE, TRUE, FALSE))

pB <- ggplot(quad, aes(x, y)) +
  geom_tile(aes(fill = hit), colour = "white", linewidth = 1.4, width = 0.98, height = 0.98) +
  geom_text(aes(label = lab, fontface = ifelse(hit, "bold", "plain"),
                colour = hit), size = 2.5, lineheight = 0.95) +
  annotate("text", x = 0.36, y = 1, label = "DTU\nyes", size = 2.4,
           colour = "grey30", lineheight = 0.9) +
  annotate("text", x = 0.36, y = 2, label = "DTU\nno", size = 2.4,
           colour = "grey30", lineheight = 0.9) +
  annotate("text", x = 1, y = 0.36, label = "DGE no", size = 2.4, colour = "grey30") +
  annotate("text", x = 2, y = 0.36, label = "DGE yes", size = 2.4, colour = "grey30") +
  annotate("text", x = 1, y = 1.34, label = "IsoGraph's layer",
           size = 2.4, colour = HILITE, fontface = "italic") +
  scale_fill_manual(values = c(`FALSE` = "grey93", `TRUE` = "#F6E0D2"), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = HILITE), guide = "none") +
  coord_cartesian(xlim = c(0.22, 2.55), ylim = c(0.22, 2.6), clip = "off") +
  theme_blank()

# ===========================================================================
# Panel C - the method, as one left-to-right row
# ===========================================================================
steps <- data.frame(
  x = c(1, 2, 3, 4),
  lab = c("per-gene\nswitch features",
          "VAE\nlatent",
          "gene-gene\ngraph",
          "co-switch\nmodules"),
  hit = c(FALSE, FALSE, FALSE, TRUE))

pC <- ggplot(steps, aes(x, 1)) +
  geom_tile(aes(fill = hit), colour = "white", linewidth = 1.2,
            width = 0.86, height = 0.72) +
  geom_text(aes(label = lab, colour = hit,
                fontface = ifelse(hit, "bold", "plain")),
            size = 2.4, lineheight = 0.95) +
  annotate("segment", x = c(1.45, 2.45, 3.45), xend = c(1.55, 2.55, 3.55),
           y = 1, yend = 1, linewidth = 0.5, colour = "grey45",
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  scale_fill_manual(values = c(`FALSE` = "grey93", `TRUE` = "#F6E0D2"), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = HILITE), guide = "none") +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.55, 1.45), clip = "off") +
  theme_blank()

# ===========================================================================
# Panels D/E - it works where the truth is known
# ===========================================================================
raw <- as.data.frame(read_parquet(
  rel("01_synthetic_benchmark", "01_synthetic", "_m", "synthetic_results.parquet")))

MAIN_METHODS <- c("isograph_baseline", "isograph_latent", "isograph_graph",
                  "isograph_vae", "isograph_spearman_leiden", "wgcna_gene")
METHOD_SHORT <- c(isograph_baseline = "Baseline", isograph_latent = "Latent",
                  isograph_graph = "Graph", isograph_vae = "IsoGraph VAE",
                  isograph_spearman_leiden = "Spearman-Leiden", wgcna_gene = "WGCNA")
CORE_SCEN <- c(idealized_switching = "Idealized\nswitching",
               noise_stress = "Noise\nstress",
               feature_space_interactions = "Feature\ninteractions",
               non_switching_background = "Non-switching\nbackground",
               unequal_isoform_abundance = "Unequal\nabundance",
               negative_control_noise = "Negative\ncontrol")

b <- raw |>
  filter(run_method %in% MAIN_METHODS, run_scenario %in% names(CORE_SCEN)) |>
  mutate(method = factor(unname(METHOD_SHORT[run_method]), unname(METHOD_SHORT)),
         scenario = factor(unname(CORE_SCEN[run_scenario]), unname(CORE_SCEN)),
         is_vae = run_method == "isograph_vae")

mk_box <- function(df, yvar, ylab) {
  ggplot(df, aes(method, .data[[yvar]], fill = is_vae)) +
    geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.35, linewidth = 0.28,
                 width = 0.7, colour = "grey25") +
    facet_wrap(~ scenario, nrow = 1) +
    scale_fill_manual(values = c(`FALSE` = MUTED, `TRUE` = HILITE), guide = "none") +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
    labs(x = NULL, y = ylab) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.4),
          strip.text = element_text(size = 7, lineheight = 0.9),
          strip.background = element_blank(),
          panel.spacing.x = unit(3, "pt"))
}

pD <- mk_box(b, "metrics_module_recovery", "Module recovery\n(AUC; 1 = perfect)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pE <- mk_box(b, "metrics_switch_gene_detection_rate",
             "Switch-gene detection\n(recall)") +
  theme(strip.text = element_blank())

cat(sprintf("  benchmark panels: %d runs, %d methods, %d scenarios\n",
            nrow(b), dplyr::n_distinct(b$method), dplyr::n_distinct(b$scenario)))

# ===========================================================================
# Panel F - what the real data showed, as observed vs its own comparator
# ===========================================================================
# The PI asked Fig 1 to preview the paper's RESULTS, not only its method. The honest way
# to do that in one panel is to put every real-data claim on a common fraction scale
# beside the thing it must beat -- a matched WGCNA baseline where the question is "does
# the method add anything", a permutation/abundance-matched null where the question is
# "is this above chance". Claims whose natural unit is not a fraction (the switch-unique
# partial-R2 ratio, the within-donor cis-control counts) are named in the caption rather
# than forced onto this axis.
#
# Every value is read from a result file. Nothing here is typed in by hand, because this
# panel is the one a reader checks first.

# -- trusted-module rate, both methods ---------------------------------------
stab_dir <- rel("04_module_trust", "_m", "stability", "module_trust")
trusted_rate <- function(meth) {
  fs <- list.files(stab_dir, pattern = paste0("^module_stability__.*__", meth, "\\.parquet$"),
                   full.names = TRUE)
  d <- bind_rows(lapply(fs, function(f) as.data.frame(read_parquet(f))))
  sum(d$trusted) / nrow(d)
}

# -- within-cohort split-half age sign concordance ---------------------------
within_conc <- function(meth) {
  fs <- list.files(stab_dir, pattern = paste0("^within_cohort__.*__", meth, "\\.parquet$"),
                   full.names = TRUE)
  d <- bind_rows(lapply(fs, function(f) as.data.frame(read_parquet(f))))
  d <- d[d$both_age_sig, ]
  sum(d$sign_concordant) / nrow(d)
}

# -- cross-cohort eigengene projection ---------------------------------------
proj <- as.data.frame(read_parquet(rel("04_module_trust", "_m", "stability",
                                       "eigengene_projection",
                                       "eigengene_projection_summary.parquet")))
proj_rate <- function(meth, dir) {
  r <- proj[proj$method == meth & proj$direction == dir, ]
  r$sign_match / r$n_age_testable
}

# -- long-read confirmation, observed vs its abundance-matched null ----------
lr <- function(f, field) {
  j <- jsonlite::fromJSON(rel("06_switch_mechanism", "_m", "switch_orthogonal_confirm", f))
  j$matched_null$switch_like_rate[[field]]
}

ev <- tibble::tribble(
  ~claim,                                            ~observed,                             ~comparator,                           ~ctype,
  "Modules chance-trusted",                          trusted_rate("isograph"),              trusted_rate("wgcna"),                 "WGCNA baseline",
  "Split-half age sign concordance",                 within_conc("isograph"),               within_conc("wgcna"),                  "WGCNA baseline",
  "Aging transfers, BrainSEQ\u2192GTEx",               proj_rate("isograph", "brainseq_to_gtex"), proj_rate("wgcna", "brainseq_to_gtex"), "WGCNA baseline",
  "Aging transfers, GTEx\u2192BrainSEQ",               proj_rate("isograph", "gtex_to_brainseq"), proj_rate("wgcna", "gtex_to_brainseq"), "WGCNA baseline",
  "Long-read switch-like, anchored pairs",           lr("anchored_summary.json", "observed"),   lr("anchored_summary.json", "null_mean"),   "Matched null",
  "Long-read switch-like, all pairs",                lr("global_null_summary.json", "observed"), lr("global_null_summary.json", "null_mean"), "Matched null"
) |>
  mutate(claim = factor(claim, rev(claim)))

ev_long <- bind_rows(
  transmute(ev, claim, value = observed,   what = "IsoGraph / observed"),
  transmute(ev, claim, value = comparator, what = ctype)
) |>
  mutate(what = factor(what, c("IsoGraph / observed", "WGCNA baseline", "Matched null")))

pF <- ggplot(ev, aes(y = claim)) +
  geom_segment(aes(x = comparator, xend = observed, yend = claim),
               colour = "grey70", linewidth = 0.45) +
  geom_point(data = ev_long, aes(x = value, colour = what), size = 1.9) +
  geom_text(aes(x = observed, label = sprintf("%.2f", observed)),
            vjust = -1.05, size = 2.2, colour = HILITE) +
  scale_colour_manual(values = c(`IsoGraph / observed` = HILITE,
                                 `WGCNA baseline` = "#0072B2",
                                 `Matched null` = "grey45"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "Fraction (rate, or concordant modules)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.8),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
        legend.position = "bottom", legend.text = element_text(size = 6.8),
        legend.key.size = unit(0.3, "cm"))

cat("  results panel:\n")
for (i in seq_len(nrow(ev))) cat(sprintf("    %-38s %.3f vs %.3f (%s)\n",
    ev$claim[i], ev$observed[i], ev$comparator[i], ev$ctype[i]))

# ===========================================================================
# Assemble: schematic row on top (the definition), benchmark rows beneath
# ===========================================================================
# F sits last: definition (A-C), synthetic validation (D-E), then what the real data
# showed (F). A reader who stops after Fig 1 should still know what was found.
design <- "AAABBB
           AAABBB
           CCCCCC
           DDDDDD
           DDDDDD
           EEEEEE
           EEEEEE
           FFFFFF
           FFFFFF"
fig <- wrap_plots(pA, pB, pC, pD, pE, pF, design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

# Nature Communications caps a figure at 180 x 247 mm, i.e. 7.09 x 9.72 in. 9.5 leaves a
# margin for the caption block without shrinking the benchmark panels further.
save_fig(fig, "figConceptOverview", width = 7.09, height = 9.5)
cat("Done. Output in", FIG_DIR, "\n")
