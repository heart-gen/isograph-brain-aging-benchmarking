# Per-gene sQTL-vs-eQTL colocalization contrast

Does a disease association colocalize with a switch gene's **splicing** QTL more
than with the **same gene's expression** QTL, at the same locus, in the same GTEx
brain tissue? This puts the splicing-specificity claim on per-gene footing; the
`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.

Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics
(not credible sets), against the per-locus GWAS summary statistics already built
for the CLPP layer. Because the comparison is made within a gene, the module
membership that selected the gene cancels between the two modalities.

**Gene pool for this run (`--arm wgcna_multiplex`):** genes in the **matched WGCNA (multiplex) phenotype-significant modules**, built on the same switch+abundance features.

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
| aging__ad | ad | 2056 | 26 | 52 | 13 | 39 | 0.00041 | 0.00123 | 0.24 | 0.231 | 0.62 |
| aging__als | als | 972 | 20 | 25 | 13 | 18 | 0.473 | 0.473 | 0.286 | 0.276 | 0.59 |
| aging__lbd | lbd | 369 | 1 | 4 | 1 | 4 | 0.375 | 0.45 | 0.268 | 0.275 | 0.0842 |
| aging__pd | pd | 871 | 12 | 21 | 6 | 15 | 0.0784 | 0.157 | 0.24 | 0.241 | 0.217 |
| aging__scz | scz | 4399 | 78 | 140 | 46 | 108 | 6.33e-07 | 3.8e-06 | 0.262 | 0.247 | 0.000441 |
| brainseq-sczd__scz | scz | 976 | 21 | 30 | 16 | 25 | 0.211 | 0.317 | 0.28 | 0.271 | 0.488 |
| **POOLED** | ALL | 9643 | 158 | 272 | 95 | 209 | 5.46e-11 | - | 0.259 | 0.249 | 0.013 |

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
| primary | 1e-05 | 0.8 | 100 | 9643 | 95 | 209 | 5.46e-11 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 9643 | 68 | 149 | 3.94e-08 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 9643 | 19 | 43 | 0.00316 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 9643 | 253 | 448 | 1.65e-13 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 9629 | 95 | 209 | 5.46e-11 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 70 | 0 | 2 | 0.5 |
| aging__ad | NA | 1721 | 8 | 34 | 6.88e-05 |
| aging__ad | True | 265 | 5 | 3 | 0.727 |
| aging__als | False | 34 | 1 | 0 | 1 |
| aging__als | NA | 815 | 8 | 13 | 0.383 |
| aging__als | True | 123 | 4 | 5 | 1 |
| aging__lbd | False | 8 | 0 | 0 | n/a |
| aging__lbd | NA | 312 | 0 | 3 | 0.25 |
| aging__lbd | True | 49 | 1 | 1 | 1 |
| aging__pd | False | 27 | 0 | 1 | 1 |
| aging__pd | NA | 719 | 4 | 10 | 0.18 |
| aging__pd | True | 125 | 2 | 4 | 0.688 |
| aging__scz | False | 149 | 2 | 3 | 1 |
| aging__scz | NA | 3571 | 39 | 83 | 8.41e-05 |
| aging__scz | True | 679 | 5 | 22 | 0.00151 |
| brainseq-sczd__scz | False | 19 | 0 | 1 | 1 |
| brainseq-sczd__scz | NA | 917 | 16 | 21 | 0.511 |
| brainseq-sczd__scz | True | 40 | 0 | 3 | 0.25 |

## Scope and limits

* 740,511 `coloc.abf` fits; 114,026 paired (gene, locus, tissue) cells;
  9,643 (gene, locus) pairs after collapsing tissues by maximum.
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
