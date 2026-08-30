# Publication figure for the IsoGraph module trust funnel (real data).
# Q1 stability -> Q2 driver reproducibility -> Q3 cross-cohort aging replication ->
# Q4 structural-switch drivers. Reads 03_module_trust/_m/stability/module_trust/*.parquet and
# writes figTrustFunnel.{pdf,png} to 03_module_trust/_m/stability/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/trust_funnel_figure.R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT    <- find_root()
rel     <- function(...) file.path(ROOT, ...)
MT_DIR  <- rel("03_module_trust", "_m", "stability", "module_trust")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito: IsoGraph vermillion, WGCNA reddish purple (shared with synthetic_benchmark.R).
METHOD_COLORS <- c(isograph = "#D55E00", wgcna = "#CC79A7")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna = "WGCNA")
TRUST_COLORS  <- c(`TRUE` = "#D55E00", `FALSE` = "grey75")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.38, "cm"),
      strip.text         = element_text(size = 8.3),
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

# Concatenate every per-region file for a funnel stage (method encoded in the filename).
load_stage <- function(prefix) {
  files <- list.files(MT_DIR, pattern = paste0("^", prefix, ".*\\.parquet$"), full.names = TRUE)
  files <- files[!grepl("_pooled__", files)]   # pooled tables are empty by the n=3 floor
  bind_rows(lapply(files, function(f) {
    df <- as.data.frame(read_parquet(f))
    method <- if (grepl("__wgcna\\.parquet$", f)) "wgcna" else "isograph"
    df$method <- method
    df
  }))
}

REGION_LABELS <- c(
  caudate = "Caudate", hippocampus = "Hippocampus", dlpfc = "DLPFC",
  caudate_basal_ganglia = "Caudate (GTEx)", frontal_cortex_ba9 = "Frontal ctx (GTEx)"
)
relabel_region <- function(x) ifelse(x %in% names(REGION_LABELS), REGION_LABELS[x], x)

# ---------------------------------------------------------------------------
# Panel A - Q1 stability: per-module co-assignment density vs permutation null
# ---------------------------------------------------------------------------
stab <- load_stage("module_stability__")
trust_counts <- stab |>
  group_by(method) |>
  summarise(trusted = sum(trusted), M = n(), .groups = "drop") |>
  mutate(lab = sprintf("%s trusted\n%d / %d", METHOD_LABELS[method], trusted, M))
null_band <- stab |> group_by(method) |>
  summarise(null = median(null_mean, na.rm = TRUE), .groups = "drop")

pA <- ggplot(stab, aes(method, coassign_density)) +
  geom_violin(aes(fill = method), colour = NA, alpha = 0.18, scale = "width") +
  geom_jitter(aes(colour = trusted), width = 0.18, height = 0, size = 0.5, alpha = 0.7) +
  geom_crossbar(data = null_band, aes(x = method, y = null, ymin = null, ymax = null),
                width = 0.55, linewidth = 0.4, colour = "grey25", linetype = "dashed") +
  geom_text(data = trust_counts, aes(x = method, y = 1.02, label = lab),
            size = 2.5, vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  scale_colour_manual(values = TRUST_COLORS, labels = c(`TRUE` = "Trusted (FDR<0.05)",
                                                        `FALSE` = "Not trusted")) +
  scale_x_discrete(labels = METHOD_LABELS) +
  coord_cartesian(ylim = c(0, 1.18), clip = "off") +
  labs(x = NULL, y = "Co-assignment density") +
  theme_pub() + theme(legend.position = "bottom",
                      plot.margin = margin(12, 6, 4, 4, "pt"))

# ---------------------------------------------------------------------------
# Panel B - Q2 driver-loading reproducibility (IsoGraph only; WGCNA has no tx drivers)
# ---------------------------------------------------------------------------
wc <- load_stage("within_cohort__") |> filter(method == "isograph", is.finite(driver_load_rho))
wc$region_lab <- factor(relabel_region(wc$region),
                        levels = relabel_region(names(REGION_LABELS)))
wc_med <- wc |> group_by(region_lab) |>
  summarise(med = median(driver_load_rho), .groups = "drop")

pB <- ggplot(wc, aes(region_lab, driver_load_rho)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_violin(fill = METHOD_COLORS[["isograph"]], colour = NA, alpha = 0.22, scale = "width") +
  geom_boxplot(width = 0.16, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_text(data = wc_med, aes(region_lab, 1.06, label = sprintf("%.2f", med)),
            size = 2.4) +
  coord_cartesian(ylim = c(-0.6, 1.12)) +
  labs(x = NULL, y = expression("Driver-loading "*rho*" (split-half)")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------------
# Panel C - Q3 cross-cohort aging concordance (paired BrainSEQ vs GTEx age effect).
# Plots the linear/covariate-free arm; the linear-vs-spline asymmetry is reported in
# REPLICATION_PERMUTATION.md rather than by switching the panel's model.
# ---------------------------------------------------------------------------
rep <- load_stage("module_aging_replication__") |>
  filter(is.finite(age_effect_bs), is.finite(age_effect_gtex))
rep$method_lab <- METHOD_LABELS[rep$method]
# Permutation P for the arm this panel actually plots: the linear/Pearson statistic
# against the `matching` null (the stricter of the two nulls -- it holds every age
# statistic fixed and permutes only which GTEx module each BrainSEQ module is matched to).
# Read from the stats JSON rather than transcribed, so the annotation cannot drift from
# replication_permutation.py's output.
perm_p <- function(method) {
  f <- file.path(MT_DIR, sprintf(
    "replication_permutation__%s__pearson__matching__stats.json", method))
  if (!file.exists(f)) return(NA_real_)
  jsonlite::fromJSON(f)$p_emp
}
# "concordant", not "replicate": the pre-registered decision rule in
# REPLICATION_PERMUTATION.md forbids the word "replication" for this count, which does not
# separate from the matching null under the covariate-adjusted model.
rep_counts <- rep |> group_by(method) |>
  summarise(rep = sum(replicates), M = n(), .groups = "drop") |>
  mutate(p_emp = vapply(method, perm_p, numeric(1)),
         lab = sprintf("%s: %d / %d concordant%s", METHOD_LABELS[method], rep, M,
                       ifelse(is.finite(p_emp),
                              sprintf("  (perm P = %.3f)", p_emp), "")))
lim <- max(abs(c(rep$age_effect_bs, rep$age_effect_gtex)), na.rm = TRUE)

pC <- ggplot(rep, aes(age_effect_bs, age_effect_gtex)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted", colour = "grey60") +
  geom_point(aes(colour = method, alpha = replicates, size = replicates)) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.3), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.7), guide = "none") +
  annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1, size = 2.4, lineheight = 0.95,
           label = paste(rep_counts$lab, collapse = "\n")) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "BrainSEQ age effect", y = "GTEx age effect") +
  theme_pub() + theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Panel D - Q4 structural class of driver isoform switches (age-significant modules)
# ---------------------------------------------------------------------------
comp <- load_stage("module_complementarity__") |> filter(method == "isograph", age_sig)
SWITCH_COLS <- c(drv_cds_changed = "CDS change", drv_utr_changed = "UTR change",
                 drv_biotype_switch = "Biotype switch",
                 drv_coding_status_change = "Coding-status change")
dd <- comp |>
  select(region, all_of(names(SWITCH_COLS))) |>
  pivot_longer(-region, names_to = "switch", values_to = "frac") |>
  filter(is.finite(frac)) |>
  group_by(switch) |>
  summarise(mean_frac = mean(frac), se = sd(frac) / sqrt(n()), .groups = "drop") |>
  mutate(switch = factor(SWITCH_COLS[switch], levels = unname(SWITCH_COLS)))

pD <- ggplot(dd, aes(switch, mean_frac)) +
  geom_col(fill = METHOD_COLORS[["isograph"]], width = 0.68, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, mean_frac - se), ymax = pmin(1, mean_frac + se)),
                width = 0.2, linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Mean driver fraction") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------------
# Assemble: (A | B) / (C | D), full width
# ---------------------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figTrustFunnel", width = 7.2, height = 6.6)
cat("Done. Output in", FIG_DIR, "\n")
