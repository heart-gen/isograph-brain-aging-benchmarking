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
| aging__ad | ad | 1890 | 29 | 50 | 17 | 38 | 0.00646 | 0.0129 | 0.236 | 0.234 | 0.921 |
| aging__als | als | 927 | 20 | 22 | 12 | 14 | 0.845 | 0.845 | 0.284 | 0.273 | 0.245 |
| aging__lbd | lbd | 345 | 1 | 3 | 1 | 3 | 0.625 | 0.75 | 0.266 | 0.268 | 0.259 |
| aging__pd | pd | 789 | 11 | 22 | 6 | 17 | 0.0347 | 0.052 | 0.256 | 0.255 | 0.502 |
| aging__scz | scz | 4254 | 75 | 129 | 42 | 96 | 4.94e-06 | 2.96e-05 | 0.258 | 0.238 | 6.83e-06 |
| brainseq-sczd__scz | scz | 1694 | 41 | 69 | 23 | 51 | 0.00152 | 0.00455 | 0.281 | 0.26 | 0.054 |
| **POOLED** | ALL | 9899 | 177 | 295 | 101 | 219 | 3.71e-11 | - | 0.259 | 0.247 | 7.56e-05 |

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
| primary | 1e-05 | 0.8 | 100 | 9899 | 101 | 219 | 3.71e-11 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 9899 | 72 | 157 | 1.97e-08 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 9899 | 21 | 50 | 0.000767 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 9899 | 272 | 493 | 1.15e-15 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 9886 | 101 | 219 | 3.71e-11 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 227 | 2 | 5 | 0.453 |
| aging__ad | NA | 1506 | 14 | 27 | 0.0596 |
| aging__ad | True | 157 | 1 | 6 | 0.125 |
| aging__als | False | 94 | 2 | 2 | 1 |
| aging__als | NA | 746 | 8 | 11 | 0.648 |
| aging__als | True | 87 | 2 | 1 | 1 |
| aging__lbd | False | 37 | 0 | 0 | n/a |
| aging__lbd | NA | 273 | 1 | 2 | 1 |
| aging__lbd | True | 35 | 0 | 1 | 1 |
| aging__pd | False | 96 | 0 | 4 | 0.125 |
| aging__pd | NA | 632 | 6 | 9 | 0.607 |
| aging__pd | True | 61 | 0 | 4 | 0.125 |
| aging__scz | False | 561 | 6 | 13 | 0.167 |
| aging__scz | NA | 3270 | 34 | 71 | 0.000391 |
| aging__scz | True | 423 | 2 | 12 | 0.0129 |
| brainseq-sczd__scz | False | 13 | 0 | 0 | n/a |
| brainseq-sczd__scz | NA | 1570 | 22 | 47 | 0.00354 |
| brainseq-sczd__scz | True | 111 | 1 | 4 | 0.375 |

## Scope and limits

* 760,128 `coloc.abf` fits; 117,491 paired (gene, locus, tissue) cells;
  9,899 (gene, locus) pairs after collapsing tissues by maximum.
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
