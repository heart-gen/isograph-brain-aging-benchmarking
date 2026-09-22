# Supplementary real-data figure: BrainSEQ -> GTEx module aging replication, and why it is
# reported as a limitation rather than as a method result.
#
# The two cohorts are quantified by different pipelines -- BrainSEQ transcripts by Salmon,
# GTEx by RSEM -- so a module whose age effect fails to carry across them has failed a test
# that confounds quantification with biology. The PI ruled on 2026-09-19 that this is not a
# fair replication test; within-cohort split-half concordance carries the reproducibility
# claim (figTrustFunnel panel D) and this figure documents the cross-cohort attempt.
#
# (A) paired age effects, BrainSEQ vs GTEx, for matched modules of both methods.
# (B) the observed concordant-module count against its permutation null: neither method
#     separates from chance, so the panel reports a null result for both.
# Reads 04_module_trust/_m/stability/module_trust/{module_aging_replication__*,
# replication_permutation__*}.parquet|json; writes figCrossCohortReplication.{pdf,png}.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/crosscohort_replication_figure.R
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(ggplot2); library(patchwork); library(jsonlite)
})

find_root <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, ".here")) && d != dirname(d)) d <- dirname(d)
  d
}
ROOT   <- find_root(); rel <- function(...) file.path(ROOT, ...)
MT_DIR <- rel("04_module_trust", "_m", "stability", "module_trust")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

METHOD_COLORS <- c(isograph = "#D55E00", wgcna = "#CC79A7")
METHOD_LABELS <- c(isograph = "IsoGraph", wgcna = "WGCNA")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(size = 7.5), axis.title = element_text(size = 8.5),
          legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
          legend.key.size = unit(0.36, "cm"),
          panel.grid.major = element_line(linewidth = 0.3, colour = "grey92"),
          plot.margin = margin(4, 6, 4, 4, "pt"))
}
save_fig <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300)
  cat("  ", name, " saved\n", sep = "")
}

load_stage <- function(prefix) {
  files <- list.files(MT_DIR, pattern = paste0("^", prefix, ".*\\.parquet$"),
                      full.names = TRUE)
  files <- files[!grepl("_pooled__", files)]
  bind_rows(lapply(files, function(f) {
    df <- as.data.frame(read_parquet(f))
    df$method <- if (grepl("__wgcna\\.parquet$", f)) "wgcna" else "isograph"
    df
  }))
}

# ---------------------------------------------------------------------------
# Panel A - paired age effects across the two cohorts
# ---------------------------------------------------------------------------
rep <- load_stage("module_aging_replication__") |>
  filter(is.finite(age_effect_bs), is.finite(age_effect_gtex))

perm <- function(method, field) {
  f <- file.path(MT_DIR, sprintf("replication_permutation__%s__pearson__matching__stats.json",
                                 method))
  if (!file.exists(f)) return(NA) else fromJSON(f)[[field]]
}
counts <- rep |> group_by(method) |>
  summarise(rep = sum(replicates), M = n(), .groups = "drop") |>
  mutate(p_emp = vapply(method, function(m) as.numeric(perm(m, "p_emp")), numeric(1)),
         lab = sprintf("%s  %d/%d, P = %s", METHOD_LABELS[method], rep, M,
                       ifelse(is.finite(p_emp), sprintf("%.2f", p_emp), "NA")))
lab_a <- paste(c("Concordant modules", counts$lab), collapse = "\n")
lim <- max(abs(c(rep$age_effect_bs, rep$age_effect_gtex)), na.rm = TRUE) * 1.04

pA <- ggplot(rep, aes(age_effect_bs, age_effect_gtex)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted",
              colour = "grey60") +
  geom_point(aes(colour = method, alpha = replicates, size = replicates)) +
  scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.3), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.7), guide = "none") +
  annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1, size = 2.3, lineheight = 0.95,
           label = lab_a) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "BrainSEQ age effect (Salmon)", y = "GTEx age effect (RSEM)") +
  theme_pub() + theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# Panel B - observed count against the matching permutation null
# The matching null is the stricter of the two: it holds every age statistic fixed and
# permutes only which GTEx module each BrainSEQ module is paired with.
# ---------------------------------------------------------------------------
nulls <- bind_rows(lapply(names(METHOD_LABELS), function(m) {
  f <- file.path(MT_DIR, sprintf("replication_permutation__%s__pearson__matching.parquet", m))
  if (!file.exists(f)) return(NULL)
  data.frame(method = m, T = as.data.frame(read_parquet(f))$T)
}))
obs <- data.frame(method = names(METHOD_LABELS)) |>
  mutate(T_obs = vapply(method, function(m) as.numeric(perm(m, "T_obs")), numeric(1)),
         p_emp = vapply(method, function(m) as.numeric(perm(m, "p_emp")), numeric(1)),
         lab = sprintf("observed %d\nP = %.2f", T_obs, p_emp),
         facet = METHOD_LABELS[method])
nulls$facet <- METHOD_LABELS[nulls$method]

pB <- ggplot(nulls, aes(T)) +
  geom_bar(aes(y = after_stat(prop)), fill = "grey80", width = 0.85) +
  geom_vline(data = obs, aes(xintercept = T_obs, colour = method), linewidth = 0.7) +
  geom_text(data = obs, aes(x = T_obs, y = 0.42, label = lab, colour = method),
            hjust = -0.12, size = 2.4, lineheight = 0.95) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = METHOD_COLORS, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_width(2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Concordant modules under the matching null", y = "Fraction of permutations") +
  theme_pub() + theme(strip.background = element_blank(),
                      strip.text = element_text(size = 8, face = "bold"),
                      panel.grid.major.x = element_blank())

fig <- (pA | pB) + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figCrossCohortReplication", width = 7.2, height = 3.4)
cat(sprintf("  cross-cohort: %s\n", paste(counts$lab, collapse = "; ")))
cat("Done. Output in", FIG_DIR, "\n")
