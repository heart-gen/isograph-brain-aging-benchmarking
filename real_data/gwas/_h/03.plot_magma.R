#!/usr/bin/env Rscript
# Visualize MAGMA gene-set enrichment results for IsoGraph vs WGCNA modules.
#
# Input:  real_data/gwas/_m/results/{trait}_noMHC_{backend}.gsa.out
# Output: real_data/gwas/_m/figures/
#   magma_enrichment_dotplot.pdf   -- dot plot: module × GWAS trait, size=BETA_STD, color=-log10(P)
#   magma_enrichment_top20.pdf     -- top 20 modules per trait ranked by P-value
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
})

script_dir <- dirname(normalizePath(if (interactive()) getwd() else commandArgs()[4], mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
res_dir <- file.path(project_root, "real_data", "gwas", "_m", "results")
out_dir <- file.path(project_root, "real_data", "gwas", "_m", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TRAITS   <- c("scz", "mdd", "bp", "ad", "pd", "stroke")
TRAIT_LABELS <- c(
    scz = "SCZ",
    mdd = "MDD",
    bp = "BP",
    ad = "AD",
    pd = "PD",
    stroke = "Stroke"
)
# IsoGraph backend dir is selectable (default canonical isograph_vae). With a
# non-canonical backend (e.g. isograph_vae_res5) the combined parquet + figures
# get a matching suffix so the canonical artifacts are never overwritten.
ISOGRAPH_BACKEND <- Sys.getenv("MAGMA_ISOGRAPH_BACKEND", "isograph_vae")
SUFFIX <- if (ISOGRAPH_BACKEND == "isograph_vae") "" else sub("^isograph_vae", "", ISOGRAPH_BACKEND)
BACKENDS <- c(ISOGRAPH_BACKEND, "wgcna_gene")
FDR_THRESHOLD <- 0.05
TOP_N <- 20

# ── Load all .gsa.out files ────────────────────────────────────────────────────
load_gsa <- function() {
    rows <- list()
    for (backend in BACKENDS) {
        for (trait in TRAITS) {
            path <- file.path(res_dir, sprintf("%s_noMHC_%s.gsa.out", trait, backend))
            if (!file.exists(path)) next
            df <- read.table(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
            df$trait   <- TRAIT_LABELS[[trait]]
            df$backend <- backend
            rows[[length(rows) + 1]] <- df
        }
    }
    if (length(rows) == 0) stop("No .gsa.out files found in ", res_dir)
    do.call(rbind, rows)
}

gsa <- load_gsa()
gsa$trait <- factor(gsa$trait, levels = unname(TRAIT_LABELS))
cat(sprintf("Loaded %d rows from %d files\n", nrow(gsa), length(unique(paste(gsa$backend, gsa$trait)))))

# ── FDR correction across all tests per backend ───────────────────────────────
gsa <- gsa %>%
    group_by(backend) %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%
    ungroup()

# Save combined results
write_parquet(gsa, file.path(file.path(project_root, "real_data", "gwas", "_m"),
                             sprintf("magma_results_combined%s.parquet", SUFFIX)))

# ── Dot plot: top modules per backend × trait ─────────────────────────────────
top_modules <- gsa %>%
    filter(FDR <= FDR_THRESHOLD) %>%
    group_by(backend, trait) %>%
    slice_min(P, n = TOP_N) %>%
    ungroup()

plot_backend <- function(df, backend_label) {
    if (nrow(df) == 0) {
        return(ggplot() +
            labs(title = backend_label, subtitle = sprintf("No modules at FDR≤%.2f", FDR_THRESHOLD)) +
            theme_void())
    }

    # Keep only modules significant in ≥1 trait
    sig_sets <- unique(df$VARIABLE)
    df_plot <- df %>% filter(VARIABLE %in% sig_sets) %>%
        mutate(label = gsub("__", " / ", VARIABLE),
               nlp = -log10(P))

    ggplot(df_plot, aes(x = trait, y = reorder(label, nlp), color = nlp, size = BETA_STD)) +
        geom_point(alpha = 0.85) +
        scale_color_gradient(
            low = "#FEE0D2", high = "#CB181D",
            name = expression(-log[10](P))
        ) +
        scale_size_continuous(name = "Std. beta", range = c(2, 8)) +
        labs(title = backend_label,
             x = "GWAS trait", y = NULL) +
        theme_bw(base_size = 9) +
        theme(axis.text.y = element_text(size = 7),
              panel.grid.minor = element_blank())
}

p1 <- plot_backend(top_modules %>% filter(backend == ISOGRAPH_BACKEND), "IsoGraph VAE modules")
p2 <- plot_backend(top_modules %>% filter(backend == "wgcna_gene"),    "Gene-level WGCNA modules")

fig <- p1 | p2
dotplot_path <- file.path(out_dir, sprintf("magma_enrichment_dotplot%s.pdf", SUFFIX))
ggsave(dotplot_path, fig, width = 14, height = 10)
cat("Saved:", dotplot_path, "\n")

# ── Top-20 ranked plot per trait ──────────────────────────────────────────────
plot_top20 <- function(trait_name) {
    df <- gsa %>% filter(trait == trait_name) %>%
        group_by(backend) %>% slice_min(P, n = TOP_N) %>% ungroup()
    if (nrow(df) == 0) return(NULL)

    df$label <- paste0(gsub(ISOGRAPH_BACKEND, "ISO", gsub("wgcna_gene", "WGCNA", df$backend)),
                       ": ", gsub("__", " / ", df$VARIABLE))
    df$label <- factor(df$label, levels = df$label[order(df$P, decreasing = TRUE)])

    fill_vals <- setNames(c("#2166AC", "#B2182B"), c(ISOGRAPH_BACKEND, "wgcna_gene"))
    fill_labs <- setNames(c("IsoGraph VAE", "WGCNA gene"), c(ISOGRAPH_BACKEND, "wgcna_gene"))
    ggplot(df, aes(x = -log10(P), y = label, fill = backend)) +
        geom_col(alpha = 0.8) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        scale_fill_manual(values = fill_vals, labels = fill_labs, name = "Method") +
        labs(title = paste("Top modules:", trait_name),
             x = expression(-log[10](P)), y = NULL) +
        theme_bw(base_size = 9) +
        theme(axis.text.y = element_text(size = 7))
}

trait_plots <- lapply(unique(gsa$trait), plot_top20)
trait_plots <- Filter(Negate(is.null), trait_plots)
if (length(trait_plots) > 0) {
    fig2 <- Reduce(`/`, trait_plots)
    top20_path <- file.path(out_dir, sprintf("magma_enrichment_top20%s.pdf", SUFFIX))
    ggsave(top20_path, fig2, width = 10, height = 4 * length(trait_plots))
    cat("Saved:", top20_path, "\n")
}
