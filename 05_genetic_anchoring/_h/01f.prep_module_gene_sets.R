#!/usr/bin/env Rscript
# Prepare MAGMA gene-set files from IsoGraph and WGCNA modules.
#
# Input:  02_module_discovery/{gtex,brainseq}/<region>/_m/{backend}/modules.parquet
# Output: 05_genetic_anchoring/_m/gwas/gene_sets/
#   isograph_vae_gene_sets.txt   -- one line per module: "region__M000  ENTREZID1 ENTREZID2 ..."
#   wgcna_gene_gene_sets.txt     -- same for gene-level WGCNA modules
suppressPackageStartupMessages({
    library(arrow)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(dplyr)
})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- dirname(normalizePath(if (interactive()) getwd() else .script_path, mustWork = FALSE))
# The wrappers export ISOGRAPH_BENCHMARK_ROOT; fall back to the repo root two
# levels above <stage>/_h/ so a direct Rscript call still resolves correctly.
project_root <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = "")
if (!nzchar(project_root)) project_root <- file.path(script_dir, "..", "..")
project_root <- normalizePath(project_root, mustWork = TRUE)
if (!dir.exists(file.path(project_root, "isograph_benchmark"))) {
    stop("project_root is not the repo root: ", project_root)
}
out_dir <- file.path(project_root, "05_genetic_anchoring", "_m", "gwas", "gene_sets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_GENES <- 10
# The IsoGraph backend dir is selectable so a non-canonical resolution (e.g.
# isograph_vae_res5) can be analysed without clobbering the canonical
# isograph_vae gene sets -- the backend name is embedded in the output filename.
# wgcna_gene is resolution-independent, so only the canonical run emits it: a non-canonical
# backend (isograph_vae_res2) can then run beside the canonical one without rewriting its files.
ISOGRAPH_BACKEND <- Sys.getenv("MAGMA_ISOGRAPH_BACKEND", "isograph_vae")
BACKENDS <- if (ISOGRAPH_BACKEND == "isograph_vae") c(ISOGRAPH_BACKEND, "wgcna_gene") else ISOGRAPH_BACKEND

# Region collections: (dataset_label, region_dir, results_root)
COLLECTIONS <- list(
    list(label = "gtex",    root = file.path(project_root, "02_module_discovery", "gtex"),
         regions = c("amygdala", "anterior_cingulate_cortex_ba24", "caudate_basal_ganglia",
                     "cerebellar_hemisphere", "cerebellum", "cortex", "frontal_cortex_ba9",
                     "hippocampus", "hypothalamus", "nucleus_accumbens_basal_ganglia",
                     "putamen_basal_ganglia", "spinal_cord_cervical_c_1", "substantia_nigra")),
    list(label = "brainseq", root = file.path(project_root, "02_module_discovery", "brainseq"),
         regions = c("caudate", "hippocampus", "dlpfc", "caudate_sczd"))
)

# Map Ensembl gene IDs (with or without version) to Entrez IDs
ensembl_to_entrez <- function(ensembl_ids) {
    clean_ids <- sub("\\.[0-9]+$", "", ensembl_ids)  # strip version suffix
    entrez <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = clean_ids,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
    )
    as.integer(na.omit(unique(entrez)))
}

safe_label <- function(x) {
    x <- gsub("[^A-Za-z0-9_]", "_", x)
    x <- gsub("_+", "_", x)
    gsub("^_|_$", "", x)
}

build_gene_set_file <- function(backend) {
    all_lines <- character(0)
    n_sets <- 0L

    for (coll in COLLECTIONS) {
        for (region in coll$regions) {
            mod_path <- file.path(coll$root, region, "_m", backend, "modules.parquet")
            if (!file.exists(mod_path)) {
                cat(sprintf("  %s / %s / %s: no modules.parquet, skipped\n",
                            coll$label, region, backend))
                next
            }

            mods <- read_parquet(mod_path)
            module_ids <- unique(mods$module_id)
            cat(sprintf("  %s / %s / %s: %d modules\n", coll$label, region, backend, length(module_ids)))

            for (mod_id in module_ids) {
                gene_ids <- mods$gene_id[mods$module_id == mod_id]
                entrez <- ensembl_to_entrez(gene_ids)
                if (length(entrez) < MIN_GENES) next

                set_name <- paste0(safe_label(coll$label), "__",
                                   safe_label(region), "__",
                                   safe_label(mod_id))
                all_lines <- c(all_lines, paste(c(set_name, entrez), collapse = "\t"))
                n_sets <- n_sets + 1L
            }
        }
    }

    # An empty gene set file is never a valid result: MAGMA would silently fall
    # back to whatever stale file is already on disk. Fail instead of writing it.
    if (n_sets == 0L) {
        stop("no gene sets built for backend ", backend,
             " -- no modules.parquet found under ", project_root)
    }

    out_path <- file.path(out_dir, paste0(backend, "_gene_sets.txt"))
    writeLines(all_lines, out_path)
    cat(sprintf("Wrote %d gene sets → %s\n", n_sets, out_path))
    invisible(out_path)
}

for (backend in BACKENDS) {
    cat(sprintf("\n=== Backend: %s ===\n", backend))
    build_gene_set_file(backend)
}
