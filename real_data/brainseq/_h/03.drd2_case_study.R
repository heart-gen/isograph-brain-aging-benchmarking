#!/usr/bin/env Rscript
# DRD2 isoform-switch case study: IsoGraph vs gene-level WGCNA in BrainSEQ caudate.
#
# Biological claim: IsoGraph's switch coordinate for DRD2 captures the D2L/D2S
# balance and places DRD2 in a module enriched for dopaminergic/synaptic genes.
# Gene-level WGCNA reflects only total DRD2 abundance and assigns DRD2 to a
# broader, less specific module — consistent with the published finding that
# D2L and D2S junctions don't cluster together (Nat Neurosci 2022, doi:10.1038/
# s41593-022-01182-7).
#
# Outputs: real_data/brainseq/caudate_sczd/_m/drd2_case_study/
#   drd2_module_comparison.parquet   -- module assignments and sizes
#   isograph_module_enrich.parquet   -- GO enrichment for DRD2 IsoGraph module
#   wgcna_module_enrich.parquet      -- GO enrichment for DRD2 WGCNA module
#   drd2_case_study.pdf              -- comparison figure
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(gprofiler2)
    library(ggplot2)
    library(patchwork)
})

script_dir <- dirname(normalizePath(if (interactive()) getwd() else commandArgs()[4], mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)

# ── Paths ──────────────────────────────────────────────────────────────────────
iso_dir   <- file.path(project_root, "real_data", "brainseq", "caudate_sczd", "_m", "isograph_vae")
wgcna_dir <- file.path(project_root, "real_data", "brainseq", "caudate_sczd", "_m", "wgcna_gene")
out_dir   <- file.path(project_root, "real_data", "brainseq", "caudate_sczd", "_m", "drd2_case_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tx_path <- file.path(project_root, "inputs", "bundles", "brainseq_sczd", "caudate", "transcripts.parquet")

# DRD2 ENSG ID (Gencode v26)
DRD2_GENE <- "ENSG00000149295"
TOP_N_ENRICH <- 15

# ── Load data ──────────────────────────────────────────────────────────────────
iso_mods   <- read_parquet(file.path(iso_dir, "modules.parquet"))
wgcna_mods <- read_parquet(file.path(wgcna_dir, "modules.parquet"))
tx_tbl     <- read_parquet(tx_path)

# ── DRD2 module lookup ─────────────────────────────────────────────────────────
# Verify DRD2 is in each module set
iso_drd2   <- iso_mods[iso_mods$gene_id == DRD2_GENE, ]
wgcna_drd2 <- wgcna_mods[wgcna_mods$gene_id == DRD2_GENE, ]

if (nrow(iso_drd2) == 0)   stop("DRD2 not found in IsoGraph modules. Check gene ID.")
if (nrow(wgcna_drd2) == 0) stop("DRD2 not found in WGCNA modules. Check gene ID.")

iso_mod_id   <- iso_drd2$module_id[1]
wgcna_mod_id <- wgcna_drd2$module_id[1]

cat(sprintf("DRD2 module: IsoGraph = %s | WGCNA = %s\n", iso_mod_id, wgcna_mod_id))

# ── Module gene lists ──────────────────────────────────────────────────────────
iso_genes   <- iso_mods$gene_id[iso_mods$module_id == iso_mod_id]
wgcna_genes <- wgcna_mods$gene_id[wgcna_mods$module_id == wgcna_mod_id]

# Strip version suffixes for gprofiler (ENSG00000149295.14 → ENSG00000149295)
strip_ver <- function(ids) sub("\\.[0-9]+$", "", ids)
iso_genes_clean   <- strip_ver(iso_genes)
wgcna_genes_clean <- strip_ver(wgcna_genes)

# ── Summary table ──────────────────────────────────────────────────────────────
comparison <- data.frame(
    method     = c("IsoGraph", "WGCNA_gene"),
    module_id  = c(iso_mod_id, wgcna_mod_id),
    n_genes    = c(length(iso_genes), length(wgcna_genes)),
    n_overlap  = rep(length(intersect(iso_genes_clean, wgcna_genes_clean)), 2),
    stringsAsFactors = FALSE
)
write_parquet(comparison, file.path(out_dir, "drd2_module_comparison.parquet"))
print(comparison)

# ── GO enrichment (gprofiler2) ─────────────────────────────────────────────────
run_enrichment <- function(genes, label) {
    cat(sprintf("Running enrichment for %s (%d genes)...\n", label, length(genes)))
    res <- gost(
        query = genes,
        organism = "hsapiens",
        sources = c("GO:BP", "GO:MF", "KEGG", "REAC"),
        correction_method = "fdr",
        significant = TRUE,
        evcodes = FALSE,
        user_threshold = 0.05
    )
    if (is.null(res)) {
        cat(sprintf("  No significant terms for %s\n", label))
        return(data.frame())
    }
    res$result[order(res$result$p_value), c("source", "term_id", "term_name", "p_value", "intersection_size")]
}

iso_enrich   <- run_enrichment(iso_genes_clean, "IsoGraph")
wgcna_enrich <- run_enrichment(wgcna_genes_clean, "WGCNA_gene")

if (nrow(iso_enrich) > 0)   write_parquet(iso_enrich,   file.path(out_dir, "isograph_module_enrich.parquet"))
if (nrow(wgcna_enrich) > 0) write_parquet(wgcna_enrich, file.path(out_dir, "wgcna_module_enrich.parquet"))

# ── Figure ─────────────────────────────────────────────────────────────────────
plot_enrich <- function(enrich_df, title, n = TOP_N_ENRICH, color_hex) {
    if (nrow(enrich_df) == 0) return(ggplot() + labs(title = title, subtitle = "No significant terms"))
    df <- head(enrich_df, n)
    df$term_name <- factor(df$term_name, levels = rev(df$term_name))
    ggplot(df, aes(x = -log10(p_value), y = term_name, size = intersection_size)) +
        geom_point(color = color_hex, alpha = 0.8) +
        scale_size_continuous(name = "Overlap size", range = c(2, 8)) +
        labs(title = title,
             x = expression(-log[10](FDR)),
             y = NULL) +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank())
}

p1 <- plot_enrich(iso_enrich,   sprintf("IsoGraph: %s module (%d genes)", iso_mod_id, length(iso_genes)),   color_hex = "#2166AC")
p2 <- plot_enrich(wgcna_enrich, sprintf("WGCNA: %s module (%d genes)", wgcna_mod_id, length(wgcna_genes)), color_hex = "#B2182B")

fig <- p1 / p2 +
    plot_annotation(
        title = "DRD2 isoform-switch module: IsoGraph vs gene-level WGCNA",
        subtitle = sprintf("BrainSEQ caudate (Control + SCZD) | DRD2 = %s", DRD2_GENE),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )

ggsave(file.path(out_dir, "drd2_case_study.pdf"), fig, width = 10, height = 12)
cat("Saved:", file.path(out_dir, "drd2_case_study.pdf"), "\n")

# ── DRD2 transcript summary ────────────────────────────────────────────────────
drd2_tx <- tx_tbl[grepl(DRD2_GENE, tx_tbl$gene_id, fixed = TRUE), ]
cat(sprintf("\nDRD2 transcripts in bundle: %d\n", nrow(drd2_tx)))
if (nrow(drd2_tx) > 0) print(drd2_tx[, intersect(c("transcript_id", "transcript_name", "transcript_type"), names(drd2_tx))])
