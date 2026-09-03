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
| aging__ad | ad | 1114 | 18 | 29 | 9 | 20 | 0.0614 | 0.123 | 0.24 | 0.238 | 0.907 |
| aging__als | als | 516 | 15 | 16 | 9 | 10 | 1 | 1 | 0.284 | 0.288 | 0.832 |
| aging__lbd | lbd | 214 | 1 | 3 | 1 | 3 | 0.625 | 0.75 | 0.263 | 0.274 | 0.291 |
| aging__pd | pd | 508 | 5 | 13 | 1 | 9 | 0.0215 | 0.0645 | 0.242 | 0.243 | 0.631 |
| aging__scz | scz | 2472 | 54 | 94 | 28 | 68 | 5.46e-05 | 0.000328 | 0.265 | 0.248 | 0.0208 |
| brainseq-sczd__scz | scz | 551 | 15 | 23 | 11 | 19 | 0.2 | 0.301 | 0.288 | 0.273 | 0.432 |
| **POOLED** | ALL | 5375 | 108 | 178 | 59 | 129 | 3.59e-07 | - | 0.261 | 0.253 | 0.15 |

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
| primary | 1e-05 | 0.8 | 100 | 5375 | 59 | 129 | 3.59e-07 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 5375 | 44 | 98 | 6.82e-06 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 5375 | 12 | 32 | 0.00366 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 5375 | 160 | 257 | 2.34e-06 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 5366 | 59 | 129 | 3.59e-07 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 52 | 0 | 1 | 1 |
| aging__ad | NA | 812 | 4 | 16 | 0.0118 |
| aging__ad | True | 250 | 5 | 3 | 0.727 |
| aging__als | False | 25 | 1 | 0 | 1 |
| aging__als | NA | 379 | 4 | 6 | 0.754 |
| aging__als | True | 112 | 4 | 4 | 1 |
| aging__lbd | False | 6 | 0 | 0 | n/a |
| aging__lbd | NA | 162 | 0 | 2 | 0.5 |
| aging__lbd | True | 46 | 1 | 1 | 1 |
| aging__pd | False | 21 | 0 | 0 | n/a |
| aging__pd | NA | 371 | 0 | 5 | 0.0625 |
| aging__pd | True | 116 | 1 | 4 | 0.375 |
| aging__scz | False | 103 | 2 | 2 | 1 |
| aging__scz | NA | 1741 | 21 | 45 | 0.00427 |
| aging__scz | True | 628 | 5 | 21 | 0.00249 |
| brainseq-sczd__scz | False | 12 | 0 | 0 | n/a |
| brainseq-sczd__scz | NA | 512 | 11 | 17 | 0.345 |
| brainseq-sczd__scz | True | 27 | 0 | 2 | 0.5 |

## Scope and limits

* 409,875 `coloc.abf` fits; 64,178 paired (gene, locus, tissue) cells;
  5,375 (gene, locus) pairs after collapsing tissues by maximum.
* `coloc.abf` assumes a **single causal variant** per trait per window. Where two
  independent causal variants sit in one window it under-calls sharing, for both
  arms alike. The CLPP layer (`10.coloc_clpp.R`) does not make this assumption and
  remains the estimator of record for *whether* a gene colocalizes at all.
* Only loci passing the GWAS window threshold enter, so this says nothing about
  sub-threshold signal.
* The rsID bridge (GTEx v8 WGS lookup) covers 98.6% of v11 variants. Dropped
  variants are dropped from both arms of a gene identically.

Regenerate: `05_genetic_anchoring/_h/19.coloc_modality_prep.sh`,
`20.coloc_modality_abf.sh` (array over 13 tissues), `21.coloc_modality_meta.sh`.
