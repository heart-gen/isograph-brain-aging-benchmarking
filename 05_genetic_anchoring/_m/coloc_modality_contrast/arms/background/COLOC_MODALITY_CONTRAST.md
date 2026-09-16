# Per-gene sQTL-vs-eQTL colocalization contrast

Does a disease association colocalize with a switch gene's **splicing** QTL more
than with the **same gene's expression** QTL, at the same locus, in the same GTEx
brain tissue? This puts the splicing-specificity claim on per-gene footing; the
`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.

Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics
(not credible sets), against the per-locus GWAS summary statistics already built
for the CLPP layer. Because the comparison is made within a gene, the module
membership that selected the gene cancels between the two modalities.

**Gene pool for this run (`--arm background`):** **every testable gene** in the same locus windows, regardless of module membership -- the arm that says whether the switch-gene result is specific to switch genes or is just what any gene at a brain GWAS locus does in GTEx.

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
| aging__ad | ad | 2614 | 34 | 60 | 19 | 45 | 0.00156 | 0.00313 | 0.234 | 0.224 | 0.0521 |
| aging__als | als | 1323 | 24 | 29 | 14 | 19 | 0.487 | 0.487 | 0.28 | 0.268 | 0.154 |
| aging__lbd | lbd | 486 | 1 | 5 | 1 | 5 | 0.219 | 0.263 | 0.263 | 0.256 | 0.696 |
| aging__pd | pd | 1065 | 12 | 29 | 6 | 23 | 0.00232 | 0.00347 | 0.248 | 0.246 | 0.556 |
| aging__scz | scz | 5824 | 90 | 170 | 52 | 132 | 3.28e-09 | 1.97e-08 | 0.257 | 0.231 | 5.72e-09 |
| brainseq-sczd__scz | scz | 2312 | 50 | 87 | 30 | 67 | 0.000219 | 0.000657 | 0.279 | 0.255 | 0.014 |
| **POOLED** | ALL | 13624 | 211 | 380 | 122 | 291 | 5.01e-17 | - | 0.257 | 0.24 | 1.1e-09 |

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
| primary | 1e-05 | 0.8 | 100 | 13624 | 122 | 291 | 5.01e-17 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 13624 | 85 | 202 | 3.8e-12 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 13624 | 26 | 61 | 0.000224 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 13624 | 326 | 662 | 5.12e-27 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 13607 | 122 | 291 | 5.01e-17 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 233 | 2 | 5 | 0.453 |
| aging__ad | NA | 2210 | 16 | 34 | 0.0153 |
| aging__ad | True | 171 | 1 | 6 | 0.125 |
| aging__als | False | 96 | 2 | 2 | 1 |
| aging__als | NA | 1134 | 10 | 16 | 0.327 |
| aging__als | True | 93 | 2 | 1 | 1 |
| aging__lbd | False | 40 | 0 | 0 | n/a |
| aging__lbd | NA | 407 | 1 | 4 | 0.375 |
| aging__lbd | True | 39 | 0 | 1 | 1 |
| aging__pd | False | 101 | 0 | 4 | 0.125 |
| aging__pd | NA | 897 | 6 | 15 | 0.0784 |
| aging__pd | True | 67 | 0 | 4 | 0.125 |
| aging__scz | False | 584 | 6 | 13 | 0.167 |
| aging__scz | NA | 4779 | 44 | 107 | 3.09e-07 |
| aging__scz | True | 461 | 2 | 12 | 0.0129 |
| brainseq-sczd__scz | False | 14 | 0 | 0 | n/a |
| brainseq-sczd__scz | NA | 2159 | 29 | 63 | 0.000509 |
| brainseq-sczd__scz | True | 139 | 1 | 4 | 0.375 |

## Scope and limits

* 1,058,898 `coloc.abf` fits; 156,775 paired (gene, locus, tissue) cells;
  13,624 (gene, locus) pairs after collapsing tissues by maximum.
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
