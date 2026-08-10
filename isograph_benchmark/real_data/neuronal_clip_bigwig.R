#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: neuronal_clip_bigwig.R WINDOWS.tsv FILES.tsv OUTPUT.tsv")
}

window_path <- args[[1L]]
file_path <- args[[2L]]
output_path <- args[[3L]]

windows <- read.delim(window_path, stringsAsFactors = FALSE, check.names = FALSE)
files <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
required_windows <- c("window_id", "chrom", "start", "end", "strand")
required_files <- c("sample_id", "assay_role", "replicate", "strand", "path")
if (!all(required_windows %in% colnames(windows))) {
  stop("window manifest is missing required columns")
}
if (!all(required_files %in% colnames(files))) {
  stop("bigWig manifest is missing required columns")
}

window_ranges <- GRanges(
  seqnames = windows$chrom,
  ranges = IRanges(start = windows$start + 1L, end = windows$end),
  strand = windows$strand,
  window_id = windows$window_id
)

extract_file <- function(file_row) {
  track_path <- file_row[["path"]]
  if (!file.exists(track_path)) {
    stop(sprintf("missing bigWig: %s", track_path))
  }
  expected_strand <- if (file_row[["strand"]] == "pos") "+" else "-"
  selected <- which(as.character(strand(window_ranges)) == expected_strand)
  if (length(selected) == 0L) {
    return(data.frame())
  }
  track <- import(BigWigFile(track_path))
  if (!("score" %in% colnames(mcols(track)))) {
    stop(sprintf("bigWig lacks score values: %s", track_path))
  }
  score <- as.numeric(track$score)
  score[!is.finite(score)] <- 0
  library_total <- sum(width(track) * score)
  query <- window_ranges[selected]
  hits <- findOverlaps(query, track, ignore.strand = TRUE)
  weighted <- numeric(length(query))
  if (length(hits) > 0L) {
    overlap_width <- width(pintersect(query[queryHits(hits)], track[subjectHits(hits)]))
    contribution <- overlap_width * score[subjectHits(hits)]
    weighted <- rowsum(contribution, queryHits(hits), reorder = FALSE)
    expanded <- numeric(length(query))
    expanded[as.integer(rownames(weighted))] <- weighted[, 1L]
    weighted <- expanded
  }
  data.frame(
    window_id = as.character(mcols(query)$window_id),
    sample_id = as.character(file_row[["sample_id"]]),
    assay_role = as.character(file_row[["assay_role"]]),
    replicate = as.integer(file_row[["replicate"]]),
    file_strand = as.character(file_row[["strand"]]),
    mean_signal = weighted / width(query),
    library_total = library_total,
    bigwig_path = track_path,
    stringsAsFactors = FALSE
  )
}

results <- vector("list", nrow(files))
for (index in seq_len(nrow(files))) {
  message(sprintf("extracting %d/%d: %s", index, nrow(files), files$sample_id[[index]]))
  results[[index]] <- extract_file(as.list(files[index, , drop = FALSE]))
}
result <- do.call(rbind, results)
write.table(result, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
