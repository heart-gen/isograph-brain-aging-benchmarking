# Supplementary figure: two corroborations of the colocalization layer (added 2026-09-29).
# (a) SMR/HEIDI on the primary confirmatory probe of each colocalizing sQTL gene, by trait:
#     how many colocalizing probes SMR also supports with HEIDI not rejecting a single
#     shared variant, against those HEIDI rejects, those SMR leaves non-significant and
#     those with no instrument. `no instrument` is untested, never negative.
# (b) BrainSEQ switch-axis (S_g) against abundance-axis (A_g) colocalization per analysis:
#     the discordant genes and the exact McNemar P. BrainSEQ is the discovery cohort, so
#     this is same-cohort anchoring on IsoGraph's own axes, not replication; it reads
#     beside the GTEx sQTL/eQTL contrast (Table S26), which leans the same way.
# Reads manuscript/_m/supp_tables/{tableS30_smr_heidi,tableS31_brainseq_axis_coloc_contrast}.csv
# (CSV, because the ledger parquets are zstd); writes figColocCorroboration.{pdf,png}.
# Run: python manuscript/_h/assemble_supp_tables.py && \
#      Rscript manuscript/_h/coloc_corroboration_figure.R
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
TAB_DIR <- rel("manuscript", "_m", "supp_tables")
FIG_DIR <- rel("manuscript", "_m", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

need_tab <- function(name) {
  f <- file.path(TAB_DIR, name)
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun: python manuscript/_h/assemble_supp_tables.py", call. = FALSE)
  }
  read.csv(f)
}

theme_pub <- function(base_size = 8.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text          = element_text(size = 7.5),
      axis.title         = element_text(size = 8.5),
      legend.text        = element_text(size = 7.5),
      legend.title       = element_text(size = 8),
      legend.key.size    = unit(0.36, "cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey88"),
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

TRAITS <- c("AD", "PD", "LBD", "ALS", "SCZ")

# ===========================================================================
# Panel A - SMR/HEIDI on the colocalizing primary sQTL probes
# ===========================================================================
AGREE <- c(coloc_and_smr_heidi_not_rejected = "SMR significant, HEIDI not rejected",
           coloc_and_smr_heidi_rejected     = "SMR significant, HEIDI rejected",
           coloc_and_smr_heidi_untestable   = "SMR significant, HEIDI untestable",
           coloc_smr_not_significant        = "SMR not significant",
           coloc_no_instrument              = "No instrument (untested)")
AGREE_COL <- setNames(c("#D55E00", "#56B4E9", "#9ecae1", "grey55", "grey85"), AGREE)
smr <- need_tab("tableS30_smr_heidi.csv") |>
  filter(modality == "sQTL", agreement %in% names(AGREE)) |>
  mutate(trait = factor(toupper(trait), rev(TRAITS)),
         agreement = factor(unname(AGREE[agreement]), rev(AGREE))) |>
  count(trait, agreement, .drop = FALSE)
tot_a <- smr |> group_by(trait) |>
  summarise(n_sup = sum(n[agreement == AGREE[[1]]]), n = sum(n), .groups = "drop")
cat(sprintf("  SMR: %d of %d colocalizing primary sQTL probes supported, %d HEIDI-rejected\n",
            sum(tot_a$n_sup), sum(tot_a$n),
            sum(smr$n[smr$agreement == AGREE[[2]]])))

pA <- ggplot(smr, aes(n, trait, fill = agreement)) +
  geom_col(width = 0.68, colour = NA) +
  geom_text(data = tot_a, aes(n, trait, label = sprintf("%d / %d", n_sup, n)),
            inherit.aes = FALSE, hjust = -0.1, size = 2.5, colour = "grey15") +
  scale_fill_manual(values = AGREE_COL, breaks = unname(AGREE), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.18))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Colocalizing primary sQTL probes", y = NULL) +
  theme_pub() + theme(legend.position = "right")

# ===========================================================================
# Panel B - BrainSEQ switch- vs abundance-axis colocalization
# ===========================================================================
AX <- c(switch_only = "Switch axis only (S_g)", abundance_only = "Abundance axis only (A_g)")
bc <- need_tab("tableS31_brainseq_axis_coloc_contrast.csv") |>
  mutate(label = ifelse(analysis == "POOLED", "Pooled",
                        sprintf("%s (%s)", toupper(trait),
                                ifelse(grepl("^brainseq-sczd", analysis), "SCZD", "aging"))))
pooled <- bc |> filter(analysis == "POOLED")
cat(sprintf("  BrainSEQ pooled: switch-only %d vs abundance-only %d, McNemar P = %.2g\n",
            pooled$switch_only, pooled$abundance_only, pooled$mcnemar_p))
bl <- bc |>
  filter(analysis != "POOLED") |>
  mutate(label = factor(label, rev(label))) |>
  pivot_longer(c(switch_only, abundance_only), names_to = "axis", values_to = "n") |>
  mutate(axis = factor(unname(AX[axis]), rev(unname(AX))))
plab <- bc |> filter(analysis != "POOLED") |>
  mutate(label = factor(label, levels(bl$label)),
         x = pmax(switch_only, abundance_only),
         p = sprintf("P = %s", formatC(mcnemar_p, format = "g", digits = 2)))

pB <- ggplot(bl, aes(n, label, fill = axis)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = NA) +
  geom_text(data = plab, aes(x, label, label = p), inherit.aes = FALSE,
            hjust = -0.15, size = 2.4, colour = "grey25") +
  scale_fill_manual(values = setNames(c("#D55E00", "#56B4E9"), unname(AX)),
                    breaks = unname(AX), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.35))) +
  labs(x = sprintf("Genes colocalizing on one axis only (pooled %d vs %d)",
                   pooled$switch_only, pooled$abundance_only),
       y = NULL) +
  theme_pub() + theme(legend.position = "top", legend.justification = "left",
                      legend.margin = margin(0, 0, 0, 0))

fig <- (pA / pB) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "figColocCorroboration", width = 7.0, height = 5.2)
cat("Done. Output in", FIG_DIR, "\n")
