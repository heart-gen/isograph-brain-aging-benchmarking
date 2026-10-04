# Figure 6c: exact nominated transcript pairs supported by junction fragments.
# Standalone, untagged panel for assembly beside figOrthogonalConfirm_{a,b}.
# No internal title. Same 82-mm column, font sizes and Okabe–Ito colors as panel a.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1]) else normalizePath(getwd())
input <- if (length(args) >= 2) args[2] else file.path(root, "06_switch_mechanism/_m/junction_pair_corroboration")
output <- if (length(args) >= 3) args[3] else file.path(root, "manuscript/_m/figures")
dir.create(output, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(file.path(input, "primary_pair_results.csv"))
g <- read.csv(file.path(input, "primary_gene_region_summary.csv"))
s <- fromJSON(file.path(input, "summary.json"))
d <- d |> filter(measurement_status == "measured")
g <- g |> filter(n_measured_pairs > 0) |> arrange(gene_name, region)
if (nrow(g) == 0) stop("No measurable primary pairs; do not make an empty positive panel")
stopifnot(nrow(d) == s$n_measured_pair_regions,
          all(d$n_both <= d$n_covered), all(d$fraction_both >= 0 & d$fraction_both <= 1))
g$key <- paste(g$gene_name, g$region, sep = ":")
g$y <- rev(seq_len(nrow(g)))
d$key <- paste(d$gene_name, d$region, sep = ":")
d <- d |> left_join(g[, c("key", "y")], by = "key") |>
  group_by(key) |> arrange(pair_key, .by_group = TRUE) |>
  mutate(y_point = y + if (n() == 1) 0 else seq(-0.17, 0.17, length.out = n())) |>
  ungroup()
g$pairs_label <- sprintf("%d/%d", g$n_measured_pairs, g$n_target_pairs)
g$donor_label <- ifelse(g$min_covered_donors == g$max_covered_donors,
                       as.character(g$min_covered_donors),
                       sprintf("%d–%d", g$min_covered_donors, g$max_covered_donors))
COL <- c(dlpfc = "#D55E00", hippocampus = "#0072B2", caudate = "#009E73")
LABEL <- c(dlpfc = "DLPFC", hippocampus = "Hippocampus", caudate = "Caudate")
base_theme <- theme_classic(base_size = 8.5) + theme(
  axis.text = element_text(size = 7.5, colour = "grey15"),
  axis.title = element_text(size = 8.5),
  axis.line.y = element_blank(), axis.ticks.y = element_blank(),
  legend.position = "bottom", legend.title = element_blank(),
  legend.text = element_text(size = 7), legend.key.width = unit(0.35, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey30", lineheight = 1.1),
  plot.caption.position = "plot", plot.margin = margin(7, 3, 3, 3, "pt"))
cap <- sprintf("%d nominated genes; %d region-matched; %d measured.\nPoints: pairs; diamonds: medians; lines: pair ranges.",
               s$n_nominated_genes, s$n_primary_region_genes, s$n_measured_genes)
p <- ggplot() +
  geom_segment(data = g, aes(x = 0, xend = 1, y = y, yend = y),
               colour = "grey91", linewidth = 0.3) +
  geom_segment(data = g, aes(x = min_fraction_both, xend = max_fraction_both,
                            y = y, yend = y, colour = region), linewidth = 0.65) +
  geom_point(data = d, aes(fraction_both, y_point, colour = region), size = 1, alpha = 0.75) +
  geom_point(data = g, aes(median_fraction_both, y, colour = region), shape = 18, size = 2.5) +
  geom_text(data = g, aes(x = 1.12, y = y, label = pairs_label), size = 2.35, colour = "grey25") +
  geom_text(data = g, aes(x = 1.48, y = y, label = donor_label), size = 2.35, colour = "grey25") +
  annotate("text", x = c(1.12, 1.48), y = nrow(g) + 0.7,
           label = c("Pairs", "Donors"), size = 2.35, colour = "grey25") +
  scale_colour_manual(values = COL, labels = LABEL, breaks = names(COL)) +
  scale_y_continuous(breaks = g$y, labels = g$gene_name,
                     limits = c(0.4, nrow(g) + 1), expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                     limits = c(-0.015, 1.7), expand = expansion(mult = 0)) +
  labs(x = "Donors supporting both alternatives\n(fraction of covered donors)", y = NULL, caption = cap) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  base_theme
# Keep the numeric axis under the data, not the annotation columns.
p <- p + theme(axis.line.x = element_blank()) +
  annotate("segment", x = 0, xend = 1, y = 0.4, yend = 0.4, linewidth = 0.3)
height_mm <- max(69, 42 + 4.5 * nrow(g))
name <- file.path(output, "figJunctionPairCorroboration_c")
ggsave(paste0(name, ".pdf"), p, width = 82, height = height_mm, units = "mm", device = cairo_pdf)
ggsave(paste0(name, ".png"), p, width = 82, height = height_mm, units = "mm", dpi = 400)
writeLines(capture.output(sessionInfo()), file.path(input, "plot_sessionInfo.txt"))
cat(sprintf("Saved %s.pdf (82 × %.1f mm); %d pair–region points\n", name, height_mm, nrow(d)))
