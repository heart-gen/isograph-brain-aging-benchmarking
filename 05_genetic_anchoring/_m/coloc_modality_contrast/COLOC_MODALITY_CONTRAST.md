# Per-gene sQTL-vs-eQTL colocalization contrast

Does a disease association colocalize with a switch gene's **splicing** QTL more
than with the **same gene's expression** QTL, at the same locus, in the same GTEx
brain tissue? This puts the splicing-specificity claim on per-gene footing; the
`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.

Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics
(not credible sets), against the per-locus GWAS summary statistics already built
for the CLPP layer. Because the comparison is made within a gene, the IsoGraph
module membership that selected the gene cancels between the two arms.

## Primary result

Primary arm: `p12 = 1e-5` (coloc default), colocalized at `PP4 >= 0.80`, best of
the 13 GTEx brain tissues, GTEx's own grouped-permutation representative intron,
cells with `>= 100` shared SNPs in **both** modalities.

| analysis | trait | genes | sQTL coloc | eQTL coloc | splicing-only | expression-only | McNemar P | q | median cond PP4 sQTL | eQTL | Wilcoxon P |
|---|---|---|---|---|---|---|---|---|---|---|---|
| aging__ad | ad | 348 | 9 | 9 | 5 | 5 | 1 | 1 | 0.22 | 0.215 | 0.441 |
| aging__als | als | 163 | 9 | 9 | 5 | 5 | 1 | 1 | 0.276 | 0.266 | 0.55 |
| aging__lbd | lbd | 58 | 1 | 1 | 1 | 1 | 1 | 1 | 0.239 | 0.249 | 0.174 |
| aging__pd | pd | 156 | 5 | 8 | 2 | 5 | 0.453 | 0.906 | 0.235 | 0.234 | 0.562 |
| aging__scz | scz | 855 | 16 | 35 | 7 | 26 | 0.00132 | 0.00791 | 0.25 | 0.237 | 0.673 |
| brainseq-sczd__scz | scz | 67 | 0 | 4 | 0 | 4 | 0.125 | 0.375 | 0.235 | 0.278 | 0.264 |
| **POOLED** | ALL | 1647 | 40 | 66 | 20 | 46 | 0.00186 | - | 0.241 | 0.238 | 0.388 |

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
| primary | 1e-05 | 0.8 | 100 | 1647 | 20 | 46 | 0.00186 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 1647 | 14 | 40 | 0.000535 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 1647 | 3 | 11 | 0.0574 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 1647 | 60 | 90 | 0.0176 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 1643 | 20 | 46 | 0.00186 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 72 | 0 | 2 | 0.5 |
| aging__ad | True | 276 | 5 | 3 | 0.727 |
| aging__als | False | 34 | 1 | 0 | 1 |
| aging__als | True | 129 | 4 | 5 | 1 |
| aging__lbd | False | 8 | 0 | 0 | n/a |
| aging__lbd | True | 50 | 1 | 1 | 1 |
| aging__pd | False | 27 | 0 | 1 | 1 |
| aging__pd | True | 129 | 2 | 4 | 0.688 |
| aging__scz | False | 150 | 2 | 3 | 1 |
| aging__scz | True | 705 | 5 | 23 | 0.000912 |
| brainseq-sczd__scz | False | 21 | 0 | 1 | 1 |
| brainseq-sczd__scz | True | 46 | 0 | 3 | 0.25 |

## Scope and limits

* 126,390 `coloc.abf` fits; 19,727 paired (gene, locus, tissue) cells;
  1,647 (gene, locus) pairs after collapsing tissues by maximum.
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
