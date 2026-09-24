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
# (A) DEFINITION. One gene, two isoforms, usage shifting continuously with age while
#     total gene-level abundance stays flat -- the event a DGE pipeline cannot see. The
#     shift is deliberately NOT a dominance reversal (0.80 -> 0.60 against 0.20 -> 0.40):
#     a switch is a change in usage, and requiring the minor isoform to overtake the major
#     one would define the object more narrowly than the method measures it.
# (B) WHAT IS MISSED. The DGE x DTU quadrant map: abundance pipelines cover the top row,
#     IsoGraph's contribution is the DTU-without-DGE quadrant.
# (C) METHOD. Per-gene switch features -> VAE latent -> gene-gene graph -> co-switch
#     modules.
# (D) Module recovery on synthetic ground truth, six core scenarios, six main methods.
# (E) Switch-gene detection on the same runs.
# (F) Synthetic ground-truth genetics: module recovery on identical features, the matched
#     IsoGraph-versus-WGCNA separator.
# (G) The same runs: recovery of the planted cis-variants and the null-variant false
#     positive rate, comparable across methods -- the genetic signal is detectable by all
#     of them once a module is recovered at all; what differs is the module.
#
# Panels A-C are drawn from explicit coordinates -- they are illustrations of a
# definition, and carry no data. Panels D-G are read from the committed benchmark ledger.
# (The real-data summary panel that closed this figure until 2026-09-22 now opens
# figTrustFunnel; its long-read rows moved to figOrthogonalConfirm.)
#
# Reads 01_synthetic_benchmark/01_synthetic/_m/synthetic_results.parquet.
# Writes manuscript/_m/figures/figConceptOverview.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/concept_overview_figure.R
suppressPackageStartupMessages({
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
# Panel A - the definition: usage shifts continuously, total abundance does not
# ===========================================================================
age <- seq(20, 80, length.out = 200)
# A smooth, monotone shift in usage: isoform A 0.80 -> 0.60, isoform B 0.20 -> 0.40. The
# two isoforms never trade dominance -- the earlier draft crossed them at 0.5, which reads
# as if a switch REQUIRED a reversal (and even a 0.75 -> 0.45 / 0.25 -> 0.55 shift would
# still technically cross). The total is deliberately flat so the panel makes the single
# point it exists to make.
ramp <- 1 / (1 + exp(-(age - 50) / 8))          # 0 -> 1, centred on age 50
USAGE_A <- c(from = 0.80, to = 0.60)
USAGE_B <- c(from = 0.20, to = 0.40)
u_a <- USAGE_A[["from"]] + (USAGE_A[["to"]] - USAGE_A[["from"]]) * ramp
u_b <- USAGE_B[["from"]] + (USAGE_B[["to"]] - USAGE_B[["from"]]) * ramp
defn <- bind_rows(
  data.frame(age = age, value = u_a, series = "Isoform A usage"),
  data.frame(age = age, value = u_b, series = "Isoform B usage"))
tot <- data.frame(age = age, value = rep(1.0, length(age)),
                  series = "Total gene abundance (scaled)")

pA <- ggplot() +
  geom_line(data = defn, aes(age, value, colour = series), linewidth = 0.85) +
  geom_line(data = tot, aes(age, value, colour = series), linewidth = 0.7,
            linetype = "22") +
  # the size of the shift, drawn once at the old end so the reader sees it is a change
  # in usage and not a change of the dominant isoform
  annotate("segment", x = 79, xend = 79, y = USAGE_A[["to"]], yend = USAGE_A[["from"]],
           linewidth = 0.3, colour = "grey45",
           arrow = arrow(length = unit(0.04, "in"), ends = "both", type = "closed")) +
  annotate("text", x = 21, y = 0.915, hjust = 0, size = 2.2, colour = "grey35",
           label = sprintf("usage shifts by %.2f with age; the dominant isoform is unchanged",
                           USAGE_A[["from"]] - USAGE_A[["to"]])) +
  scale_colour_manual(
    values = c(`Isoform A usage` = ISO_A, `Isoform B usage` = ISO_B,
               `Total gene abundance (scaled)` = TOTAL),
    breaks = c("Isoform A usage", "Isoform B usage", "Total gene abundance (scaled)"),
    name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1)) +
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

pD <- mk_box(b, "metrics_module_recovery", "Module recovery\n(best-match Jaccard; 1 = perfect)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pE <- mk_box(b, "metrics_switch_gene_detection_rate",
             "Switch-gene detection\n(recall)") +
  theme(strip.text = element_blank())

cat(sprintf("  benchmark panels: %d runs, %d methods, %d scenarios\n",
            nrow(b), dplyr::n_distinct(b$method), dplyr::n_distinct(b$scenario)))

# ===========================================================================
# Panels F/G - synthetic ground-truth genetics: the matched IsoGraph-vs-WGCNA separator
# ===========================================================================
# The genetic_anchoring scenario plants one bi-allelic cis-variant per co-switching module
# (dosage shifts that module's switch latent) plus an equal number of unlinked null
# variants, and runs four methods on the IDENTICAL switch+abundance feature matrix. The
# per-allele effect is modest, so the genotype -> usage signal is only recoverable by
# pooling a recovered module into an eigengene. Two things are therefore read from the
# same runs, and they must be shown together or the panel overclaims:
#
#   F  module recovery -- the decisive separator. WGCNA's correlation network dilutes the
#      planted modules on the same input; the VAE, its multiplex variant and the
#      Spearman-Leiden control all recover them.
#   G  genetic-signal detection -- NOT a separator. A diluted-but-nonzero module still
#      clears the Bonferroni threshold at n = 200, so every method recovers ~all planted
#      variants at a calibrated null-variant rate. The claim the pair supports is "only the
#      switch-aware network recovers the MODULE the genetic signal rides on", not "only
#      IsoGraph finds the genetics".
#
# Every number is read from the ledger; the paired test is recomputed here on the
# dataset_id pairing so the annotation cannot drift from the data it sits above.
GEN_METHODS <- c(isograph_vae = "IsoGraph VAE", isograph_vae_multiplex = "IsoGraph multiplex",
                 isograph_spearman_leiden = "Spearman-Leiden", wgcna_gene = "WGCNA")
gen <- raw |>
  filter(run_scenario == "genetic_anchoring", run_method %in% names(GEN_METHODS),
         is.finite(metrics_module_recovery)) |>
  mutate(method = factor(unname(GEN_METHODS[run_method]), unname(GEN_METHODS)),
         is_vae = run_method == "isograph_vae")

# paired one-sided test, IsoGraph VAE > WGCNA, on the datasets both methods completed
pair <- gen |>
  filter(run_method %in% c("isograph_vae", "wgcna_gene")) |>
  select(run_dataset_id, run_method, metrics_module_recovery) |>
  pivot_wider(names_from = run_method, values_from = metrics_module_recovery) |>
  filter(is.finite(isograph_vae), is.finite(wgcna_gene))
wt <- wilcox.test(pair$isograph_vae, pair$wgcna_gene, paired = TRUE, alternative = "greater")
delta <- mean(pair$isograph_vae - pair$wgcna_gene)
f_lab <- sprintf("IsoGraph VAE vs WGCNA, same features\n\u0394 = %+.2f, paired P = %s, n = %d",
                 delta, formatC(wt$p.value, format = "g", digits = 2), nrow(pair))

gen_box <- function(df, yvar, ylab) {
  ggplot(df, aes(method, .data[[yvar]], fill = is_vae)) +
    geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.35, linewidth = 0.28,
                 width = 0.7, colour = "grey25") +
    scale_fill_manual(values = c(`FALSE` = MUTED, `TRUE` = HILITE), guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.4),
          strip.text = element_text(size = 7, lineheight = 0.9),
          strip.background = element_blank(),
          panel.spacing.x = unit(6, "pt"))
}

pF <- gen_box(gen, "metrics_module_recovery", "Module recovery\n(best-match Jaccard; 1 = perfect)") +
  annotate("text", x = 0.55, y = 1.2, hjust = 0, vjust = 1, size = 2.1,
           lineheight = 0.95, colour = "grey25", label = f_lab) +
  scale_y_continuous(limits = c(0, 1.22), breaks = c(0, 0.5, 1))

gen_long <- gen |>
  select(method, is_vae, metrics_genetic_anchor_recall, metrics_genetic_anchor_fpr) |>
  pivot_longer(starts_with("metrics_"), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(recode(metric,
                                metrics_genetic_anchor_recall = "Planted cis-variants recovered\n(recall)",
                                metrics_genetic_anchor_fpr    = "Null variants flagged\n(false-positive rate)"),
                         c("Planted cis-variants recovered\n(recall)",
                           "Null variants flagged\n(false-positive rate)")))
pG <- gen_box(gen_long, "value", "Fraction of variants") +
  geom_hline(data = data.frame(metric = factor("Null variants flagged\n(false-positive rate)",
                                               levels(gen_long$metric)), y = 0.05),
             aes(yintercept = y), linetype = "22", linewidth = 0.3, colour = "grey45") +
  facet_wrap(~ metric, nrow = 1) +
  scale_y_continuous(limits = c(0, 1.22), breaks = c(0, 0.5, 1))

gen_sum <- gen |> group_by(run_method) |>
  summarise(n = n(), module_recovery = mean(metrics_module_recovery),
            recall = mean(metrics_genetic_anchor_recall),
            fpr = mean(metrics_genetic_anchor_fpr),
            best_r2 = mean(metrics_genetic_anchor_best_r2_mean), .groups = "drop")
cat("  synthetic-genetics panels (genetic_anchoring scenario):\n")
for (i in seq_len(nrow(gen_sum))) cat(sprintf(
  "    %-26s n=%3d  module recovery %.3f  recall %.3f  FPR %.3f  best R2 %.3f\n",
  gen_sum$run_method[i], gen_sum$n[i], gen_sum$module_recovery[i], gen_sum$recall[i],
  gen_sum$fpr[i], gen_sum$best_r2[i]))
cat(sprintf("    paired IsoGraph VAE - WGCNA module recovery: delta %+.3f, one-sided P %s, n %d\n",
            delta, formatC(wt$p.value, format = "g", digits = 3), nrow(pair)))

# ===========================================================================
# Assemble: schematic row on top (the definition), benchmark rows beneath
# ===========================================================================
# Definition (A-C), synthetic validation on the core grid (D-E), then the synthetic
# genetics arm (F-G): the one benchmark that isolates network inference from feature
# representation, and the ground-truth warrant for the real-data genetic anchoring.
design <- "AAABBB
           AAABBB
           CCCCCC
           DDDDDD
           DDDDDD
           EEEEEE
           EEEEEE
           FFGGGG
           FFGGGG"
fig <- wrap_plots(pA, pB, pC, pD, pE, pF, pG, design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

# Nature Communications caps a figure at 180 x 247 mm, i.e. 7.09 x 9.72 in. 9.5 leaves a
# margin for the caption block without shrinking the benchmark panels further.
save_fig(fig, "figConceptOverview", width = 7.09, height = 9.5)
cat("Done. Output in", FIG_DIR, "\n")
