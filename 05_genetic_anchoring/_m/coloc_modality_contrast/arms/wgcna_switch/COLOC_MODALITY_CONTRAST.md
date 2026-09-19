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
| aging__ad | ad | 1699 | 24 | 41 | 14 | 31 | 0.0161 | 0.0322 | 0.231 | 0.226 | 0.262 |
| aging__als | als | 809 | 18 | 17 | 12 | 11 | 1 | 1 | 0.276 | 0.269 | 0.786 |
| aging__lbd | lbd | 341 | 1 | 3 | 1 | 3 | 0.625 | 0.75 | 0.255 | 0.265 | 0.0948 |
| aging__pd | pd | 760 | 9 | 16 | 5 | 12 | 0.143 | 0.215 | 0.232 | 0.241 | 0.113 |
| aging__scz | scz | 3761 | 62 | 113 | 33 | 84 | 2.67e-06 | 1.6e-05 | 0.247 | 0.235 | 0.0404 |
| brainseq-sczd__scz | scz | 2552 | 53 | 92 | 32 | 71 | 0.000153 | 0.00046 | 0.27 | 0.254 | 0.0226 |
| **POOLED** | ALL | 9922 | 167 | 282 | 97 | 212 | 5.35e-11 | - | 0.251 | 0.242 | 0.144 |

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
| primary | 1e-05 | 0.8 | 100 | 9922 | 97 | 212 | 5.35e-11 |
| p12_5e-6 | 5e-06 | 0.8 | 100 | 9922 | 67 | 165 | 1.01e-10 |
| p12_1e-6 | 1e-06 | 0.8 | 100 | 9922 | 22 | 38 | 0.0519 |
| pp4_call_0.5 | 1e-05 | 0.5 | 100 | 9922 | 268 | 443 | 5.45e-11 |
| min_shared_500 | 1e-05 | 0.8 | 500 | 9911 | 97 | 212 | 5.35e-11 |

## GO-invisible split

Whether the per-gene contrast localises to the GO-invisible modules. The
set-level `qtl_anchoring` result does **not** localise there after the
2026-08-29 refresh, so this is a check, not a confirmation.

| analysis | GO-invisible | genes | splicing-only | expression-only | McNemar P |
|---|---|---|---|---|---|
| aging__ad | False | 539 | 6 | 14 | 0.115 |
| aging__ad | NA | 802 | 5 | 9 | 0.424 |
| aging__ad | True | 358 | 3 | 8 | 0.227 |
| aging__als | False | 261 | 4 | 7 | 0.549 |
| aging__als | NA | 384 | 4 | 0 | 0.125 |
| aging__als | True | 164 | 4 | 4 | 1 |
| aging__lbd | False | 94 | 0 | 2 | 0.5 |
| aging__lbd | NA | 184 | 0 | 1 | 1 |
| aging__lbd | True | 63 | 1 | 0 | 1 |
| aging__pd | False | 230 | 3 | 6 | 0.508 |
| aging__pd | NA | 367 | 2 | 4 | 0.688 |
| aging__pd | True | 163 | 0 | 2 | 0.5 |
| aging__scz | False | 1343 | 16 | 38 | 0.00384 |
| aging__scz | NA | 1560 | 7 | 22 | 0.00813 |
| aging__scz | True | 858 | 10 | 24 | 0.0243 |
| brainseq-sczd__scz | False | 60 | 1 | 0 | 1 |
| brainseq-sczd__scz | NA | 2149 | 27 | 56 | 0.00193 |
| brainseq-sczd__scz | True | 343 | 4 | 15 | 0.0192 |

## Scope and limits

* 762,006 `coloc.abf` fits; 122,292 paired (gene, locus, tissue) cells;
  9,922 (gene, locus) pairs after collapsing tissues by maximum.
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
