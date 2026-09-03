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
| aging__ad | ad | 2550 | 30 | 56 | 15 | 41 | 0.000686 | 0.00206 | 0.237 | 0.222 | 0.0387 |
| aging__als | als | 1242 | 24 | 28 | 15 | 19 | 0.608 | 0.608 | 0.282 | 0.269 | 0.333 |
| aging__lbd | lbd | 483 | 1 | 5 | 1 | 5 | 0.219 | 0.328 | 0.267 | 0.257 | 0.457 |
| aging__pd | pd | 1093 | 12 | 25 | 6 | 19 | 0.0146 | 0.0293 | 0.239 | 0.236 | 0.311 |
| aging__scz | scz | 5489 | 88 | 161 | 52 | 125 | 4.05e-08 | 2.43e-07 | 0.261 | 0.238 | 8.94e-09 |
| brainseq-sczd__scz | scz | 1212 | 24 | 30 | 19 | 25 | 0.451 | 0.542 | 0.278 | 0.262 | 0.0164 |
| **POOLED** | ALL | 12069 | 179 | 305 | 108 | 234 | 8.11e-12 | - | 0.258 | 0.241 | 1.37e-08 |

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
| primary | 1e-05 | 0.8 | 100 | 12069 | 108 | 234 | 8.11e-12 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 12069 | 72 | 166 | 1.01e-09 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 12069 | 22 | 48 | 0.00255 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 12069 | 287 | 537 | 2.3e-18 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 12051 | 108 | 234 | 8.11e-12 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 72 | 0 | 2 | 0.5 |
| aging__ad | NA | 2202 | 10 | 36 | 0.000156 |
| aging__ad | True | 276 | 5 | 3 | 0.727 |
| aging__als | False | 34 | 1 | 0 | 1 |
| aging__als | NA | 1079 | 10 | 14 | 0.541 |
| aging__als | True | 129 | 4 | 5 | 1 |
| aging__lbd | False | 8 | 0 | 0 | n/a |
| aging__lbd | NA | 425 | 0 | 4 | 0.125 |
| aging__lbd | True | 50 | 1 | 1 | 1 |
| aging__pd | False | 27 | 0 | 1 | 1 |
| aging__pd | NA | 937 | 4 | 14 | 0.0309 |
| aging__pd | True | 129 | 2 | 4 | 0.688 |
| aging__scz | False | 150 | 2 | 3 | 1 |
| aging__scz | NA | 4634 | 45 | 99 | 7.95e-06 |
| aging__scz | True | 705 | 5 | 23 | 0.000912 |
| brainseq-sczd__scz | False | 21 | 0 | 1 | 1 |
| brainseq-sczd__scz | NA | 1145 | 19 | 21 | 0.875 |
| brainseq-sczd__scz | True | 46 | 0 | 3 | 0.25 |

## Scope and limits

* 938,343 `coloc.abf` fits; 138,880 paired (gene, locus, tissue) cells;
  12,069 (gene, locus) pairs after collapsing tissues by maximum.
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
