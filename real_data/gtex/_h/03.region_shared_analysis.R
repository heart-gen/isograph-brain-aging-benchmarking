#!/usr/bin/env Rscript
# Cross-region aging module overlap analysis for GTEx v11 brain.
# For each backend (isograph_vae, wgcna_gene), computes:
#   1. Pairwise Jaccard similarity of module gene sets across 13 regions
#   2. "Shared" modules (Jaccard >= 0.20 in >= half the region pairs)
#   3. Significantly age-associated modules per region (FDR <= 0.10)
# Outputs: real_data/gtex/_m/region_overlap.parquet, age_summary.parquet
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- dirname(normalizePath(if (interactive()) getwd() else .script_path, mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
out_dir <- file.path(project_root, "real_data", "gtex", "_m")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GTEX_REGIONS <- c(
    "amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
    "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
    "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
    "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra"
)
BACKENDS <- c("isograph_vae", "wgcna_gene")
FDR_THRESHOLD <- 0.10

# ── Load modules per region ────────────────────────────────────────────────────
load_modules <- function(backend) {
    result <- list()
    for (region in GTEX_REGIONS) {
        path <- file.path(project_root, "real_data", "gtex", region, "_m", backend, "modules.parquet")
        if (!file.exists(path)) next
        m <- read_parquet(path)
        result[[region]] <- split(m$gene_id, m$module_id)
    }
    result
}

# ── Jaccard similarity ─────────────────────────────────────────────────────────
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

pairwise_jaccard <- function(modules_by_region) {
    regions <- names(modules_by_region)
    rows <- list()
    for (i in seq_along(regions)) {
        for (j in seq_along(regions)) {
            if (i >= j) next
            r1 <- regions[i]; r2 <- regions[j]
            mods1 <- modules_by_region[[r1]]
            mods2 <- modules_by_region[[r2]]
            # Best Jaccard match for each module in r1
            for (m1 in names(mods1)) {
                jvals <- sapply(mods2, function(g2) jaccard(mods1[[m1]], g2))
                best_j <- max(jvals)
                best_m2 <- names(which.max(jvals))
                rows[[length(rows) + 1]] <- data.frame(
                    region1 = r1, module1 = m1,
                    region2 = r2, module2 = best_m2,
                    jaccard = best_j,
                    n_genes1 = length(mods1[[m1]]),
                    n_overlap = length(intersect(mods1[[m1]], mods2[[best_m2]])),
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    do.call(rbind, rows)
}

# ── Age association summary ────────────────────────────────────────────────────
age_summary <- function(backend) {
    rows <- list()
    for (region in GTEX_REGIONS) {
        path <- file.path(project_root, "real_data", "gtex", region, "_m", backend, "age_linear.parquet")
        if (!file.exists(path)) next
        lin <- read_parquet(path)
        lin$region <- region
        lin$backend <- backend
        rows[[length(rows) + 1]] <- lin
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

# ── Run ────────────────────────────────────────────────────────────────────────
for (backend in BACKENDS) {
    cat("Backend:", backend, "\n")
    mods <- load_modules(backend)
    if (length(mods) < 2) {
        cat("  Not enough regions with results; skipping.\n")
        next
    }

    jac <- pairwise_jaccard(mods)
    out_path <- file.path(out_dir, paste0("region_jaccard_", backend, ".parquet"))
    write_parquet(jac, out_path)
    cat(sprintf("  Pairwise Jaccard: %d pairs → %s\n", nrow(jac), out_path))

    age_df <- age_summary(backend)
    if (!is.null(age_df)) {
        sig <- age_df[!is.na(age_df$fdr) & age_df$fdr <= FDR_THRESHOLD, ]
        cat(sprintf("  Age-significant modules (FDR<=%.2f): %d\n", FDR_THRESHOLD, nrow(sig)))
        out_age <- file.path(out_dir, paste0("age_summary_", backend, ".parquet"))
        write_parquet(age_df, out_age)
    }
}
