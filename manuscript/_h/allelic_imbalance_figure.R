# Within-donor allelic test of isoform choice (supplementary figure).
#
# The BrainSEQ allele-aware junction recount compares the two haplotypes of a donor who is
# heterozygous at the switch-QTL lead: both share a nucleus, a cell mixture and an
# environment, so a difference in isoform-specific junction use between them is a cis effect
# that cell composition cannot produce. The figure carries the checks the claim rests on.
#
# (a) Same panel as Fig 5a: nominal P < 0.05 rate, lead-heterozygous test vs the same model
#     on lead-homozygous donors, with the BH q < 0.05 counts.
# (b) Quantile-quantile plot per region, all fitted pairs, both arms; lambda_GC from the
#     stage summary. Observed -log10 P is capped for display (the top caudate pair is 1e-116).
# (c) Agreement of the beta-binomial GLMM with the model-free stratified score test.
# (d) Sign agreement of the within-donor effect with the between-donor dosage association
#     from the same reads, over all fitted pairs and over q < 0.05 pairs.
# (e) Module-level cis control: per module, the fraction of fitted genes with a q < 0.05
#     pair against the module's age-association rank within its trait; rho and the
#     permutation P (module labels permuted, sizes held fixed). DLPFC has too few modules
#     per trait to rank and is not drawn. Descriptive: 7 and 14 modules.
#
# Reads 06_switch_mechanism/_m/ase_junction_switch/{<region>/allelic_test.parquet,
#       module_cis_control.parquet} and manuscript/_m/supp_tables/
#       {tableS22_allelic_imbalance_regions,tableS23_module_cis_control}.csv
#       (run manuscript/_h/assemble_supp_tables.py first).
# Writes manuscript/_m/figures/figAllelicImbalance.{pdf,png}.
# Run: Rscript manuscript/_h/allelic_imbalance_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
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
ASE     <- rel("06_switch_mechanism", "_m", "ase_junction_switch")
TAB_DIR <- rel("manuscript", "_m", "supp_tables")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
source(rel("manuscript", "_h", "allelic_panels.R"))

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(9, "pt"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 7.5, face = "bold"),
      panel.grid.major   = element_blank(),
      plot.margin        = margin(4, 6, 4, 4, "pt")
    )
}

need <- function(f) {
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun: python manuscript/_h/assemble_supp_tables.py", call. = FALSE)
  }
  read.csv(f)
}
tab <- need(file.path(TAB_DIR, "tableS22_allelic_imbalance_regions.csv"))
mod <- need(file.path(TAB_DIR, "tableS23_module_cis_control.csv"))
reg <- allelic_region_rows(tab)
REG_LEVELS <- unname(ALLELIC_REGION_LABELS)

tests <- bind_rows(lapply(names(ALLELIC_REGION_LABELS), function(r) {
  read_parquet(file.path(ASE, r, "allelic_test.parquet"),
               col_select = c("status", "gate_family", "pval", "pval_hom_null", "beta",
                              "score_z", "qval", "sign_agrees_between")) |>
    as.data.frame() |>
    filter(status == "fitted") |>
    mutate(region_lab = factor(ALLELIC_REGION_LABELS[[r]], REG_LEVELS))
}))

# (a) -----------------------------------------------------------------------
pA <- allelic_calibration_panel(tab, theme_pub)

# (b) QQ, both arms, all fitted pairs (the set the stage's lambda_GC is computed on) --------
ARMS <- c("Lead-heterozygous donors", "Homozygous-at-lead null")
Y_CAP <- 40
qq_arm <- function(p, arm) {
  p <- sort(p[!is.na(p) & p > 0])
  n <- length(p)
  data.frame(expected = -log10(ppoints(n)), observed = pmin(-log10(p), Y_CAP), arm = arm)
}
qq <- tests |>
  group_by(region_lab) |>
  group_modify(~ rbind(qq_arm(.x$pval, ARMS[1]), qq_arm(.x$pval_hom_null, ARMS[2]))) |>
  ungroup() |>
  mutate(arm = factor(arm, ARMS))
# the stage summary's lambda_GC is the reference; the recomputed one must agree with it
lam_chk <- tests |> group_by(region_lab) |>
  summarise(het = median(qchisq(pval, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1),
            hom = median(qchisq(pval_hom_null, 1, lower.tail = FALSE), na.rm = TRUE) /
              qchisq(0.5, 1))
chk <- merge(lam_chk, reg, by = "region_lab")
stopifnot(all(abs(chk$het - chk$het_lambda_gc) < 0.01),
          all(abs(chk$hom - chk$hom_null_lambda_gc) < 0.01))
lam_lab <- reg |>
  transmute(region_lab = factor(as.character(region_lab), REG_LEVELS),
            label = sprintf("λ het %.2f\nλ null %.2f", het_lambda_gc, hom_null_lambda_gc))
pB <- ggplot(qq, aes(expected, observed, colour = arm)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 0.5, alpha = 0.6, shape = 16) +
  geom_text(data = lam_lab, aes(x = 0.05, y = Y_CAP, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 2.4, lineheight = 0.9, colour = "grey20") +
  facet_wrap(~ region_lab, nrow = 1) +
  scale_colour_manual(values = setNames(c(ALLELIC_HET_COL, ALLELIC_NULL_COL), ARMS),
                      guide = "none") +
  scale_y_continuous(limits = c(0, Y_CAP)) +
  labs(x = expression(Expected ~ -log[10] * italic(P)),
       y = expression(Observed ~ -log[10] * italic(P))) +
  theme_pub()

# (c) GLMM beta vs score-test Z ---------------------------------------------------------
rho_lab <- reg |> transmute(region_lab = factor(as.character(region_lab), REG_LEVELS),
                            label = sprintf("ρ = %.2f", score_vs_glmm_spearman))
# 160 fitted pairs have no defined score statistic (score_z NA) and are not drawn.
pC <- ggplot(filter(tests, !is.na(score_z)), aes(score_z, beta)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_point(size = 0.4, alpha = 0.35, shape = 16, colour = "grey25") +
  geom_text(data = rho_lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.15,
            vjust = 1.3, size = 2.5, colour = "grey15") +
  facet_wrap(~ region_lab, nrow = 1) +
  labs(x = "Score-test Z", y = expression("GLMM " * beta * " (ALT vs REF haplotype)")) +
  theme_pub()

# (d) within- vs between-donor sign agreement -------------------------------------------
SETS <- c("All fitted pairs", "q < 0.05 pairs")
sg <- rbind(
  data.frame(region_lab = reg$region_lab, set = SETS[1], frac = reg$sign_agree_between_all,
             n = reg$pairs_fitted),
  data.frame(region_lab = reg$region_lab, set = SETS[2], frac = reg$sign_agree_between_q05,
             n = reg$pairs_q05)
) |> mutate(set = factor(set, SETS))
pD <- ggplot(sg, aes(frac, region_lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
  geom_line(aes(group = region_lab), colour = "grey75", linewidth = 0.5) +
  geom_point(aes(fill = set), shape = 21, size = 2.3, stroke = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), vjust = -0.9, size = 2.4,
            colour = "grey20") +
  scale_fill_manual(values = setNames(c("white", ALLELIC_HET_COL), SETS), name = NULL) +
  scale_x_continuous(limits = c(0.4, 1), labels = scales::percent) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.8))) +
  labs(x = "Same sign as the between-donor\ndosage association", y = NULL) +
  theme_pub() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

# (e) module cis-control rate vs module age-association rank ----------------------------
mm <- mod |>
  filter(!is.na(age_rank_within_trait)) |>
  mutate(region_lab = factor(unname(ALLELIC_REGION_LABELS[region]), REG_LEVELS))
mlab <- reg |>
  filter(!is.na(module_cis_rho)) |>
  transmute(region_lab = factor(as.character(region_lab), REG_LEVELS),
            label = sprintf("ρ = %.2f\nperm P = %.2f\n%d modules", module_cis_rho,
                            module_cis_perm_p, module_cis_modules_ranked))
mm <- mm |> filter(region_lab %in% mlab$region_lab) |> droplevels()
mlab <- mlab |> droplevels()
pE <- ggplot(mm, aes(age_rank_within_trait, cis_rate)) +
  geom_point(aes(size = genes_fitted), shape = 21, fill = ALLELIC_HET_COL, colour = "black",
             stroke = 0.25, alpha = 0.85) +
  geom_text(data = mlab, aes(x = 1.02, y = 1.02, label = label), hjust = 1, vjust = 1,
            size = 2.4, lineheight = 0.9, colour = "grey15") +
  facet_wrap(~ region_lab, nrow = 1) +
  scale_size_area(max_size = 4, breaks = c(10, 50, 100), name = "Genes fitted") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("weakest", "", "strongest")) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent) +
  labs(x = "Module age-association rank (within trait)",
       y = "Genes with a q < 0.05 pair") +
  theme_pub() +
  theme(legend.position = "right", panel.spacing.x = unit(1.6, "lines"))

# Reading order follows the text: rate and calibration (a, b), estimator agreement (c),
# between-donor agreement (d), module-level aggregation (e).
fig <- ((pA | pB) + plot_layout(widths = c(1, 2.1))) /
  ((pC | pD) + plot_layout(widths = c(1.9, 1))) / pE +
  plot_layout(heights = c(1, 1, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

W <- 7.1; H <- 9.2
ggsave(file.path(FIG_DIR, "figAllelicImbalance.pdf"), fig, width = W, height = H,
       units = "in", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "figAllelicImbalance.png"), fig, width = W, height = H,
       units = "in", dpi = 300)
cat("  figAllelicImbalance saved\n")
