## =============================================================================
## GTEx cell-type composition (MuSiC) — aging replication arm for the composition
## confound. Deconvolves the GTEx brain bulk bundles against the same Tran/LIBD snRNA
## references used for the BrainSEQ run (sex_context_brain/cell_proportion_estimate),
## producing sample_id-keyed proportions the IsoGraph `celltype_composition fractions
## gtex-aging` step joins into the incremental_association --composition test.
##
## Only regions with a defensibly matched reference are deconvolved (mapping mirrors the
## committed BrainSEQ precedent: striatum -> NAc, cortex -> DLPFC, plus AMY/sACC/HPC direct).
## Cerebellum, hypothalamus, spinal cord and substantia nigra have no matched Tran panel and
## are intentionally skipped rather than deconvolved against a mismatched reference.
##
## Prereq (login node): python -m isograph_benchmark.real_data.celltype_composition export-gtex
##   -> real_data/gtex/_m/composition/inputs/<region>_bulk.parquet
## Usage: Rscript 09.gtex_music_deconv.R <region>
## =============================================================================
suppressPackageStartupMessages({
    library("arrow")
    library("MuSiC")
    library("dplyr")
    library("SingleCellExperiment")
    library("DeconvoBuddies")
})

SEED <- 13
set.seed(SEED)

SN_DIR <- "/ocean/projects/bio260021p/shared/resources/libd-data/single-cell"
## GTEx region -> (reference rda, object name, striatal?)  striatal adds MSN D1/D2 mapping.
REF_MAP <- list(
    "amygdala"                          = list("SCE_AMY-n5_tran-etal.rda",   "sce.amy.tran",   FALSE),
    "anterior_cingulate_cortex_ba24"    = list("SCE_sACC-n5_tran-etal.rda",  "sce.sacc.tran",  FALSE),
    "frontal_cortex_ba9"                = list("SCE_DLPFC-n3_tran-etal.rda", "sce.dlpfc.tran", FALSE),
    "cortex"                            = list("SCE_DLPFC-n3_tran-etal.rda", "sce.dlpfc.tran", FALSE),
    "hippocampus"                       = list("SCE_HPC-n3_tran-etal.rda",  "sce.hpc.tran",   FALSE),
    "caudate_basal_ganglia"             = list("SCE_NAc-n8_tran-etal.rda",  "sce.nac.tran",   TRUE),
    "putamen_basal_ganglia"             = list("SCE_NAc-n8_tran-etal.rda",  "sce.nac.tran",   TRUE),
    "nucleus_accumbens_basal_ganglia"   = list("SCE_NAc-n8_tran-etal.rda",  "sce.nac.tran",   TRUE)
)

## --- cell type -> board-level label (mirrors sex_context_brain celltype_mapping.R) ------
map_celltype_to_board <- function(cellType, striatal = FALSE) {
    mapped <- case_when(
        grepl("^drop", cellType) ~ NA_character_,
        grepl("^Astro", cellType) ~ "Astro",
        grepl("^Excit", cellType) ~ "Excit",
        grepl("^Inhib", cellType) ~ "Inhib",
        grepl("^Oligo", cellType) ~ "Oligo",
        cellType %in% c("Micro", "Micro_resting") ~ "Micro",
        cellType %in% c("Macrophage", "Tcell") ~ "Immune",
        cellType == "Mural" ~ "Mural",
        cellType %in% c("OPC", "OPC_COP") ~ "OPC",
        TRUE ~ NA_character_
    )
    if (striatal) {
        mapped <- case_when(
            !is.na(mapped) ~ mapped,
            grepl("^MSN.D1", cellType) ~ "D1-SPN",
            grepl("^MSN.D2", cellType) ~ "D2-SPN",
            TRUE ~ NA_character_
        )
    }
    mapped
}

load_reference <- function(region, verbose = TRUE) {
    spec <- REF_MAP[[region]]
    e <- new.env()
    load(file.path(SN_DIR, spec[[1]]), envir = e)
    sce <- get(spec[[2]], envir = e)
    sce$cell_type <- map_celltype_to_board(sce$cellType, striatal = spec[[3]])
    sce <- sce[, !is.na(sce$cell_type)]
    sce$cell_type <- droplevels(factor(sce$cell_type))
    rownames(sce) <- rowData(sce)$gene_name
    if (verbose) { cat("Reference", spec[[2]], "cell types:\n"); print(table(sce$cell_type)) }
    sce
}

load_bulk_matrix <- function(region, comp_dir, verbose = TRUE) {
    fn <- file.path(comp_dir, "inputs", paste0(region, "_bulk.parquet"))
    if (!file.exists(fn)) stop("missing bulk parquet (run export-gtex first): ", fn)
    df <- as.data.frame(arrow::read_parquet(fn))
    gene_names <- make.unique(as.character(df$gene_name))
    samples <- setdiff(colnames(df), c("gene_id", "gene_name"))
    m <- t(as.matrix(df[, samples, drop = FALSE]))      # samples x genes
    colnames(m) <- gene_names
    rownames(m) <- samples
    if (verbose) cat("Bulk", region, ":", nrow(m), "samples x", ncol(m), "genes\n")
    m
}

generate_marker_genes <- function(sce, bulk_matrix, file_suffix) {
    ratios <- get_mean_ratio(sce, cellType_col = "cell_type", assay_name = "logcounts",
                             gene_ensembl = "gene_id", gene_name = "gene_name")
    write.csv(ratios, file = paste0("marker_stats_genes.", file_suffix, ".csv"),
              row.names = FALSE)
    ratios |>
        filter(MeanRatio.rank <= 25, gene_name %in% colnames(bulk_matrix)) |>
        pull(gene_name)
}

run_music <- function(bulk_matrix, sce, marker_genes) {
    common <- intersect(intersect(rownames(sce), colnames(bulk_matrix)), marker_genes)
    if (length(common) < 50) stop("Too few common genes: ", length(common))
    cat("MuSiC on", length(common), "common marker genes\n")
    res <- music_prop(bulk.mtx = t(bulk_matrix[, common]), sc.sce = sce[common, ],
                      clusters = "cell_type", samples = "donor", verbose = FALSE)
    res$Est.prop.weighted
}

## --- main --------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
region <- if (length(args) > 0) args[1] else stop("region argument required")
if (!region %in% names(REF_MAP))
    stop("region ", region, " has no matched Tran reference; covered: ",
         paste(names(REF_MAP), collapse = ", "))

repo <- Sys.getenv("ISOGRAPH_BENCHMARK_ROOT", unset = getwd())
comp_dir <- file.path(repo, "real_data", "gtex", "_m", "composition")
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
setwd(comp_dir)   # marker csv + proportions tsv land here

cat("==== GTEx MuSiC deconvolution:", region, "====\n")
sce         <- load_reference(region)
bulk_matrix <- load_bulk_matrix(region, comp_dir)
marker_genes <- generate_marker_genes(sce, bulk_matrix, file_suffix = paste0("gtex-", region))
cat("Marker genes:", length(marker_genes), "\n")

est_prop <- run_music(bulk_matrix, sce, marker_genes)

music_long <- est_prop |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tidyr::pivot_longer(!sample_id, names_to = "cell_type", values_to = "proportion")

data.table::fwrite(music_long,
                   file = paste0("music-proportions-gtex-", region, ".tsv"), sep = "\t")
cat("Wrote music-proportions-gtex-", region, ".tsv (",
    length(unique(music_long$sample_id)), " samples, ",
    length(unique(music_long$cell_type)), " cell types)\n", sep = "")

message("=== Reproducibility ===")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
