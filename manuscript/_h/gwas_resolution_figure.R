# Supplementary real-data figure: module size and MAGMA significance, by clustering
# resolution and by method.
#
# REWRITTEN 2026-09-19. Under the legacy expression filter this figure carried "resolution
# 5.0 removes the giant-module GWAS artifact". On the switching filter that claim REVERSED
# (PI, 2026-09-16) and production moved to resolution 2.0, so the figure now reports what is
# actually true: neither resolution escapes giant modules, and gene-level WGCNA is worse than
# either on the same gene universe.
#
# It also fixes a data bug. After the resolution switch the canonical
# `magma_results_combined.parquet` holds the PRODUCTION (res 2.0) results under the backend
# label `isograph_vae`, while the res-5.0 arm moved to `..._res5.parquet`
# (`isograph_vae_res5`). The old script read the canonical file as "res 5.0" and the _res2
# file as "res 2.0", so it plotted production twice under two different labels.
#
# (A) MAGMA significance vs module size, per group. (B) significant-module counts split by
#     giant (>=900 genes) vs not. (C) the production (res 2.0) significant sets.
#     Reads 05_genetic_anchoring/_m/gwas/magma_results_combined{,_res5}.parquet,
#     writes figGwasResolution.{pdf,png} to manuscript/_m/figures/.
# Run: /ocean/projects/bio260021p/shared/opt/envs/rnaseq/bin/Rscript \
#        manuscript/_h/gwas_resolution_figure.R
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
GWAS_DIR <- rel("05_genetic_anchoring", "_m", "gwas")
FIG_DIR  <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

GIANT     <- 900          # giant-module gene-count cutoff
FDR_THR   <- 0.05
SIG_COL   <- "#D55E00"    # vermillion = significant
NS_COL    <- "grey75"
GIANT_COL <- "#999999"    # grey  = giant module
SMALL_COL <- "#0072B2"    # blue  = non-giant module
TRAIT_COLORS <- c(SCZ = "#D55E00", BP = "#0072B2", PD = "#009E73", AD = "#E69F00",
                  MDD = "#CC79A7", Stroke = "#56B4E9")

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_blank(),
      legend.key.size    = unit(0.34, "cm"),
      strip.text         = element_text(size = 8, face = "bold"),
      strip.background   = element_rect(fill = "grey92", colour = NA),
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

# canonical = production (res 2.0); the res-5.0 arm is a disclosed sensitivity
prod <- as.data.frame(read_parquet(file.path(GWAS_DIR, "magma_results_combined.parquet")))
r5   <- as.data.frame(read_parquet(file.path(GWAS_DIR, "magma_results_combined_res5.parquet")))
stopifnot("isograph_vae" %in% prod$backend, "isograph_vae_res5" %in% r5$backend)

GROUP_LEVELS <- c("IsoGraph res 2.0\n(production)", "IsoGraph res 5.0", "WGCNA gene")
prep <- function(df, backend, group) {
  df |>
    filter(backend == !!backend) |>
    transmute(group = factor(group, GROUP_LEVELS),
              variable = FULL_NAME, trait, ngenes = NGENES, P, FDR,
              sig = FDR < FDR_THR, giant = NGENES >= GIANT)
}
# WGCNA gene modules are resolution-independent; take them once from the res-5.0 file.
dat <- bind_rows(
  prep(prod, "isograph_vae",      "IsoGraph res 2.0\n(production)"),
  prep(r5,   "isograph_vae_res5", "IsoGraph res 5.0"),
  prep(prod, "wgcna_gene",        "WGCNA gene")   # resolution-independent
)

# ---------------------------------------------------------------------------
# Panel A - module significance vs module size; sig hits in giant modules = artifact
# ---------------------------------------------------------------------------
datA <- dat |> mutate(nlp = -log10(pmax(P, 1e-300)),
                      status = ifelse(sig, "FDR < 0.05", "n.s."))
pA <- ggplot(datA, aes(ngenes, nlp)) +
  geom_vline(xintercept = GIANT, linewidth = 0.3, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = status), size = 0.7, alpha = 0.7) +
  scale_colour_manual(values = c("FDR < 0.05" = SIG_COL, "n.s." = NS_COL)) +
  scale_x_log10(breaks = c(30, 100, 300, 900, 3000, 12000),
                labels = c("30", "100", "300", "900", "3k", "12k")) +
  facet_wrap(~group, nrow = 1) +
  labs(x = "Module size (number of genes, log scale)",
       y = expression(-log[10] ~ "MAGMA" ~ italic(P))) +
  theme_pub() + theme(legend.position = c(0.07, 0.86),
                      axis.text.x = element_text(angle = 0))

# ---------------------------------------------------------------------------
# Panel B - significant-module counts split by giant vs non-giant
# ---------------------------------------------------------------------------
datB <- dat |>
  filter(sig) |>
  mutate(size_class = ifelse(giant, "Giant (>= 900 genes)", "< 900 genes")) |>
  count(group, size_class) |>
  mutate(size_class = factor(size_class, c("Giant (>= 900 genes)", "< 900 genes")))
tot <- datB |> group_by(group) |> summarise(n = sum(n), .groups = "drop")
pB <- ggplot(datB, aes(group, n, fill = size_class)) +
  geom_col(width = 0.68) +
  geom_text(data = datB |> filter(size_class == "Giant (>= 900 genes)" & n > 0),
            aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", size = 2.6) +
  geom_text(data = tot, aes(group, n, label = n), inherit.aes = FALSE,
            vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = c("Giant (>= 900 genes)" = GIANT_COL,
                               "< 900 genes" = SMALL_COL),
                    guide = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Significant modules\n(FDR < 0.05)") +
  theme_pub() + theme(legend.position = c(0.5, 0.93), legend.key.size = unit(0.3, "cm"),
                      axis.text.x = element_text(size = 6.4, lineheight = 0.85))

# ---------------------------------------------------------------------------
# Panel C - the production (res 2.0) significant sets
# ---------------------------------------------------------------------------
shorten <- function(v) {
  p <- strsplit(v, "__")
  reg <- vapply(p, function(x) gsub("_", " ", x[2]), character(1))
  mod <- vapply(p, function(x) x[length(x)], character(1))
  paste0(reg, " ", mod)
}
datC <- dat |>
  filter(group == "IsoGraph res 2.0\n(production)", sig) |>
  mutate(label = shorten(variable), nlf = -log10(FDR)) |>
  arrange(nlf)
# one module can be significant for several traits, so the label is not unique: order on a
# unique key and print the shared label on the axis (the old code made it a factor level
# directly, which errored out with "factor level is duplicated")
datC$key <- factor(make.unique(paste(datC$label, datC$trait)),
                   levels = make.unique(paste(datC$label, datC$trait)))
pC <- ggplot(datC, aes(nlf, key, colour = trait)) +
  geom_segment(aes(x = 0, xend = nlf, y = key, yend = key),
               linewidth = 0.4, colour = "grey80") +
  geom_point(aes(size = ngenes)) +
  scale_y_discrete(labels = setNames(datC$label, datC$key)) +
  scale_colour_manual(values = TRAIT_COLORS, name = "GWAS trait") +
  scale_size_area(max_size = 4.2, breaks = c(300, 900, 2000), name = "Module genes") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = expression(-log[10] ~ "FDR"), y = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2.4)),
         size = guide_legend(order = 2)) +
  theme_pub() + theme(legend.position = "right",
                      legend.title = element_text(size = 7.5),
                      legend.spacing.y = unit(0.02, "cm"),
                      axis.text.y = element_text(size = 5.6),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"))

fig <- pA / (pB | pC) +
  plot_layout(heights = c(1, 2.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "figGwasResolution", width = 7.2, height = 8.6)
cat("Done. Output in", FIG_DIR, "\n")
