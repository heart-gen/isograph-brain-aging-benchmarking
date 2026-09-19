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
| aging__ad | ad | 2889 | 34 | 68 | 17 | 51 | 4.45e-05 | 8.91e-05 | 0.228 | 0.221 | 0.585 |
| aging__als | als | 1357 | 24 | 29 | 15 | 20 | 0.5 | 0.5 | 0.271 | 0.263 | 0.0773 |
| aging__lbd | lbd | 575 | 1 | 4 | 1 | 4 | 0.375 | 0.45 | 0.251 | 0.255 | 0.254 |
| aging__pd | pd | 1240 | 13 | 29 | 7 | 23 | 0.00522 | 0.00783 | 0.232 | 0.234 | 0.169 |
| aging__scz | scz | 6282 | 92 | 173 | 52 | 133 | 2.27e-09 | 1.36e-08 | 0.246 | 0.228 | 3.97e-05 |
| brainseq-sczd__scz | scz | 4195 | 74 | 137 | 46 | 109 | 4.47e-07 | 1.34e-06 | 0.269 | 0.247 | 0.00011 |
| **POOLED** | ALL | 16538 | 238 | 440 | 138 | 340 | 9.93e-21 | - | 0.249 | 0.237 | 7.16e-07 |

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
| primary | 1e-05 | 0.8 | 100 | 16538 | 138 | 340 | 9.93e-21 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 16538 | 89 | 234 | 3.54e-16 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 16538 | 26 | 60 | 0.000317 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 16538 | 386 | 738 | 4.72e-26 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 16525 | 138 | 340 | 9.93e-21 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 835 | 8 | 21 | 0.0241 |
| aging__ad | NA | 1644 | 6 | 19 | 0.0146 |
| aging__ad | True | 410 | 3 | 11 | 0.0574 |
| aging__als | False | 381 | 6 | 9 | 0.607 |
| aging__als | NA | 786 | 5 | 7 | 0.774 |
| aging__als | True | 190 | 4 | 4 | 1 |
| aging__lbd | False | 137 | 0 | 2 | 0.5 |
| aging__lbd | NA | 363 | 0 | 2 | 0.5 |
| aging__lbd | True | 75 | 1 | 0 | 1 |
| aging__pd | False | 336 | 3 | 10 | 0.0923 |
| aging__pd | NA | 715 | 4 | 10 | 0.18 |
| aging__pd | True | 189 | 0 | 3 | 0.25 |
| aging__scz | False | 2033 | 19 | 57 | 1.48e-05 |
| aging__scz | NA | 3255 | 18 | 50 | 0.000131 |
| aging__scz | True | 994 | 15 | 26 | 0.117 |
| brainseq-sczd__scz | False | 111 | 1 | 0 | 1 |
| brainseq-sczd__scz | NA | 3638 | 41 | 94 | 5.87e-06 |
| brainseq-sczd__scz | True | 446 | 4 | 15 | 0.0192 |

## Scope and limits

* 1,269,876 `coloc.abf` fits; 194,435 paired (gene, locus, tissue) cells;
  16,538 (gene, locus) pairs after collapsing tissues by maximum.
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
