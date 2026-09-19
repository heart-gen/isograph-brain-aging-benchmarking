# Per-gene sQTL-vs-eQTL colocalization contrast

Does a disease association colocalize with a switch gene's **splicing** QTL more
than with the **same gene's expression** QTL, at the same locus, in the same GTEx
brain tissue? This puts the splicing-specificity claim on per-gene footing; the
`qtl_anchoring_meta` contrast is a set-level ratio of two odds ratios.

Estimator: `coloc::coloc.abf` on GTEx v11 cis **all-pairs** nominal statistics
(not credible sets), against the per-locus GWAS summary statistics already built
for the CLPP layer. Because the comparison is made within a gene, the module
membership that selected the gene cancels between the two modalities.

**Gene pool for this run (`--arm switch`):** IsoGraph phenotype-significant **switch genes**.

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
| aging__ad | ad | 1288 | 19 | 40 | 11 | 32 | 0.00191 | 0.00574 | 0.219 | 0.206 | 0.549 |
| aging__als | als | 599 | 17 | 20 | 10 | 13 | 0.678 | 0.813 | 0.254 | 0.248 | 0.806 |
| aging__lbd | lbd | 220 | 1 | 2 | 1 | 2 | 1 | 1 | 0.233 | 0.246 | 0.0927 |
| aging__pd | pd | 540 | 7 | 17 | 3 | 13 | 0.0213 | 0.0399 | 0.219 | 0.218 | 0.484 |
| aging__scz | scz | 3143 | 65 | 112 | 36 | 83 | 1.96e-05 | 0.000118 | 0.242 | 0.221 | 0.000847 |
| brainseq-sczd__scz | scz | 587 | 12 | 23 | 5 | 16 | 0.0266 | 0.0399 | 0.235 | 0.207 | 0.143 |
| **POOLED** | ALL | 6377 | 121 | 214 | 66 | 159 | 4.91e-10 | - | 0.234 | 0.219 | 0.00641 |

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
| primary | 1e-05 | 0.8 | 100 | 6377 | 66 | 159 | 4.91e-10 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 6377 | 49 | 120 | 4.66e-08 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 6377 | 13 | 33 | 0.00453 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 6377 | 195 | 334 | 1.6e-09 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 6365 | 66 | 159 | 4.91e-10 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 861 | 8 | 21 | 0.0241 |
| aging__ad | True | 427 | 3 | 11 | 0.0574 |
| aging__als | False | 402 | 6 | 9 | 0.607 |
| aging__als | True | 197 | 4 | 4 | 1 |
| aging__lbd | False | 142 | 0 | 2 | 0.5 |
| aging__lbd | True | 78 | 1 | 0 | 1 |
| aging__pd | False | 345 | 3 | 10 | 0.0923 |
| aging__pd | True | 195 | 0 | 3 | 0.25 |
| aging__scz | False | 2107 | 20 | 57 | 2.93e-05 |
| aging__scz | True | 1036 | 16 | 26 | 0.164 |
| brainseq-sczd__scz | False | 116 | 1 | 0 | 1 |
| brainseq-sczd__scz | True | 471 | 4 | 16 | 0.0118 |

## Scope and limits

* 488,712 `coloc.abf` fits; 76,899 paired (gene, locus, tissue) cells;
  6,377 (gene, locus) pairs after collapsing tissues by maximum.
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
