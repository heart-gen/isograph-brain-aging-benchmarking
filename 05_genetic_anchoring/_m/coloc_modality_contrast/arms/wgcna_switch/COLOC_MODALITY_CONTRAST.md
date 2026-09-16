# Per-gene sQTL-vs-eQTL colocalization contrast

Does a disease association colocalize with a switch gene's **splicing** QTL more
than with the **same gene's expression** QTL, at the same locus, in the same GTEx
brain tissue? This puts the splicing-specificity claim on per-gene footing; the
`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.

Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics
(not credible sets), against the per-locus GWAS summary statistics already built
for the CLPP layer. Because the comparison is made within a gene, the module
membership that selected the gene cancels between the two modalities.

**Gene pool for this run (`--arm wgcna_switch`):** genes in the **matched WGCNA (switch-only) phenotype-significant modules**, built on the same switch features as IsoGraph.

Loci are held FIXED across arms -- same analyses, same `loci_testable.tsv`, same
per-locus GWAS, same exclusions, same testability gate (a GTEx brain QTL credible
set). Arms differ only in which genes are tested, which is the only channel
through which a module-detection method can affect a within-gene statistic.
Compare arms with `--stage compare` (`arm_comparison.parquet`).

## Primary result

Primary arm: `p12 = 1e-5` (coloc default), colocalized at `PP4 >= 0.80`, best of
the 13 GTEx brain tissues, GTEx's own grouped-permutation representative intron,
cells with `>= 100` shared SNPs in **both** modalities.

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P | q | median cond PP4 sQTL | eQTL | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | ad | 1049 | 24 | 30 | 15 | 21 | 0.405 | 0.608 | 0.238 | 0.245 | 0.113 |
| aging__als | als | 518 | 15 | 12 | 10 | 7 | 0.629 | 0.629 | 0.288 | 0.285 | 0.602 |
| aging__lbd | lbd | 192 | 1 | 3 | 1 | 3 | 0.625 | 0.629 | 0.27 | 0.27 | 0.154 |
| aging__pd | pd | 462 | 6 | 13 | 2 | 9 | 0.0654 | 0.196 | 0.246 | 0.25 | 0.803 |
| aging__scz | scz | 2376 | 51 | 79 | 29 | 57 | 0.00335 | 0.0201 | 0.262 | 0.248 | 0.0356 |
| brainseq-sczd__scz | scz | 969 | 30 | 42 | 18 | 30 | 0.111 | 0.223 | 0.289 | 0.272 | 0.178 |
| **POOLED** | ALL | 5566 | 127 | 179 | 75 | 127 | 0.000311 | - | 0.262 | 0.257 | 0.413 |

`splicing-only` / `expression-only` are the **discordant** genes -- the ones the
McNemar test is computed on. Concordant genes (both or neither) carry no
information about which modality colocalizes and are excluded by construction.

## Two statistics, and why both are reported

* **Binary (McNemar).** Counts genes where exactly one modality colocalizes. It is
  the statistic a reader can check against the gene list.
* **Continuous (paired Wilcoxon on `PP4/(PP3+PP4)`).** The conditional posterior
  asks whether the two traits share a causal variant *given that each has one in
  the window*. eQTLs are better powered than sQTLs in GTEx, so unconditional PP4
  would favour expression for reasons unrelated to disease; conditioning divides
  that difference out. Unconditional PP4 is in `contrast.parquet` as
  `wilcoxon_pp4_p` / `median_pp4_*` so the size of the asymmetry stays visible.

## Sensitivity arms

| arm | p12 | PP4 call | min shared SNPs | pooled genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|---|---|
| primary | 1e-05 | 0.8 | 100 | 5566 | 75 | 127 | 0.000311 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 5566 | 52 | 100 | 0.000122 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 5566 | 18 | 32 | 0.0649 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 5566 | 175 | 264 | 2.51e-05 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 5556 | 75 | 127 | 0.000311 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 126 | 2 | 4 | 0.688 |
| aging__ad | NA | 777 | 13 | 12 | 1 |
| aging__ad | True | 146 | 0 | 5 | 0.0625 |
| aging__als | False | 46 | 1 | 0 | 1 |
| aging__als | NA | 398 | 7 | 6 | 1 |
| aging__als | True | 74 | 2 | 1 | 1 |
| aging__lbd | False | 26 | 0 | 0 | n/a |
| aging__lbd | NA | 129 | 1 | 2 | 1 |
| aging__lbd | True | 37 | 0 | 1 | 1 |
| aging__pd | False | 59 | 0 | 2 | 0.5 |
| aging__pd | NA | 347 | 2 | 3 | 1 |
| aging__pd | True | 56 | 0 | 4 | 0.125 |
| aging__scz | False | 305 | 3 | 8 | 0.227 |
| aging__scz | NA | 1697 | 24 | 39 | 0.0769 |
| aging__scz | True | 374 | 2 | 10 | 0.0386 |
| brainseq-sczd__scz | False | 3 | 0 | 0 | n/a |
| brainseq-sczd__scz | NA | 891 | 17 | 27 | 0.174 |
| brainseq-sczd__scz | True | 75 | 1 | 3 | 0.625 |

## Scope and limits

* 428,667 `coloc.abf` fits; 69,343 paired (gene, locus, tissue) cells;
  5,566 (gene, locus) pairs after collapsing tissues by maximum.
* `coloc.abf` assumes a **single causal variant** per trait per window. Where two
  independent causal variants sit in one window it under-calls sharing, for both
  arms alike. The CLPP layer (`03b.coloc_clpp.R`) does not make this assumption and
  remains the estimator of record for *whether* a gene colocalizes at all.
* Only loci passing the GWAS window threshold enter, so this says nothing about
  sub-threshold signal.
* The rsID bridge (GTEx v8 WGS lookup) covers 98.6% of v11 variants. Dropped
  variants are dropped from both arms of a gene identically.

Regenerate: `05_genetic_anchoring/_h/04d.coloc_modality_prep.sh`,
`05a.coloc_modality_abf.sh` (array over 13 tissues), `06a.coloc_modality_meta.sh`.
