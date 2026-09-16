# Locus event audit

Nomination source: `susie` at `PP4_sQTL >= 0.8`.

A colocalization posterior names a gene, not an event. This audit walks GWAS signal -> sQTL signal -> intron phenotype -> driver transcript, and asks whether the colocalizing intron matches a curated literature event.

## Curated events

| gene | trait | status | evidence class | interval | can promote |
|---|---|---|---|---|---|
| UNC13A | als | **reviewed** | `functional_validation` | chr19:17,641,557-17,642,844 | yes |
| PICALM | ad | **reviewed** | `twas_association` | chr11:86,026,368-86,031,468 | **no** |

Promotion to `known_mechanism_recovered` needs BOTH gates. `status: reviewed` means the coordinates were independently verified here against GENCODE v47. `evidence_class` must be `functional_validation` or `coloc_association`: a TWAS association does **not** establish a shared causal variant -- LD-driven co-regulation yields the same signal -- so a TWAS-only anchor caps a locus at `disease_locus_splice_linked` however well its coordinates match. Where that cap binds, `curated_interval_matched` still records whether the interval lined up.

## Was the curated event on trial?

sQTL arm: `representative`.

| gene | trait | curated interval | tested in GTEx | GTEx credible set |
|---|---|---|---|---|
| UNC13A | als | `chr19:17641557-17642844` | **yes** (2/13 tissues) | **no** |
| PICALM | ad | `chr11:86026368-86031468` | **yes** (10/13 tissues) | **no** |

`tested` is read from the GTEx sQTL all-pairs release, not from this pipeline's output: our output only contains introns that produced a signal pair, so inferring testability from it would call every event that failed to colocalize "untested". It is what licenses reading a non-match as a result rather than as a gap, and it is the fact that separates `context_distinct_splice_colocalization` from `disease_locus_splice_linked`.

`GTEx credible set` is reported and never gates. Its absence means there is no fine-mappable sQTL at that event in healthy brain -- which for an event that only appears under pathology is the expected biology, not a failed test.

## Tier counts

- `known_mechanism_recovered`: 0
- `context_distinct_splice_colocalization`: 0
- `disease_locus_splice_linked`: 1
- `novel_splice_led_candidate`: 39
- `not_resolved`: 0

`context_distinct_splice_colocalization` is unreachable in the `representative` arm and its count above is 0 by construction, not by result: GTEx's grouped-permutation winner is the only intron tested, so a curated non-representative event is never on trial. Only `--sqtl-arm all` can assign it.

## Evidence strength per nomination

Two descriptors travel with every nominated PP4, read at the headline tissue the gene-level number was taken from. Neither gates a tier.

- **Prior robustness**: the smallest `p12` in the sweep at which PP4 still clears 0.8. `robust` holds at 1e-6, `intermediate` from 5e-6, `primary_prior` only from the pre-specified 1e-5 upward. PP4 is monotone in p12, so this is the range the call survives rather than one arbitrary prior.
- **Why coloc.abf**, for a headline `coloc.susie` did not score: the first place the cell left the signal-level pipeline. `gwas_locus_over_max_snps` was never tested at signal level. `gwas_no_credible_set` was tested and the GWAS did not fine-map, which weakens any colocalization claimed there. `qtl_cs_not_matching_gtex` is the reference-LD artefact the agreement filter exists to remove.

| headline estimator | `robust` | `intermediate` | `primary_prior` | total |
|---|---|---|---|---|
| abf | 7 | 12 | 11 | 30 |
| susie | 0 | 9 | 1 | 10 |

| why coloc.abf | nominations |
|---|---|
| `gwas_no_credible_set` | 11 |
| `gwas_locus_over_max_snps` | 11 |
| `no_qtl_credible_set` | 4 |
| `qtl_cs_not_matching_gtex` | 4 |

| gene | trait | PP4 sQTL | headline tissue | estimator | why abf | PP4 at p12=1e-6 | calls from p12 | prior | robust tissues |
|---|---|---|---|---|---|---|---|---|---|
| SNAP91 | scz | 0.974 | Brain_Cerebellum | susie | — | 0.792 | 5e-06 | `intermediate` | 0/6 |
| UNC13A | als | 0.959 | Brain_Cerebellum | susie | — | 0.700 | 5e-06 | `intermediate` | 0/2 |
| DGKZ | scz | 0.953 | Brain_Frontal_Cortex_BA9 | susie | — | 0.671 | 5e-06 | `intermediate` | 0/9 |
| CDIP1 | scz | 0.951 | Brain_Caudate_basal_ganglia | susie | — | 0.661 | 5e-06 | `intermediate` | 0/6 |
| GGNBP2 | als | 0.912 | Brain_Cerebellum | susie | — | 0.509 | 5e-06 | `intermediate` | 0/1 |
| NDUFAF7 | scz | 0.909 | Brain_Putamen_basal_ganglia | susie | — | 0.501 | 5e-06 | `intermediate` | 0/2 |
| FGFR1 | scz | 0.908 | Brain_Caudate_basal_ganglia | susie | — | 0.497 | 5e-06 | `intermediate` | 0/2 |
| TTC19 | pd | 0.891 | Brain_Putamen_basal_ganglia | susie | — | 0.451 | 5e-06 | `intermediate` | 0/1 |
| KLC1 | scz | 0.890 | Brain_Cortex | susie | — | 0.448 | 5e-06 | `intermediate` | 0/1 |
| TPCN1 | ad | 0.829 | Brain_Cerebellum | susie | — | 0.327 | 1e-05 | `primary_prior` | 0/1 |
| ACTR1B | scz | 0.997 | Brain_Caudate_basal_ganglia | abf | `gwas_no_credible_set` | 0.967 | 1e-06 | `robust` | 9/11 |
| ACTR1B | scz | 0.997 | Brain_Caudate_basal_ganglia | abf | `gwas_no_credible_set` | 0.967 | 1e-06 | `robust` | 9/11 |
| MRPS33 | scz | 0.994 | Brain_Hypothalamus | abf | `gwas_no_credible_set` | 0.945 | 1e-06 | `robust` | 1/3 |
| IRF3 | scz | 0.994 | Brain_Cortex | abf | `gwas_locus_over_max_snps` | 0.939 | 1e-06 | `robust` | 9/11 |
| YPEL1 | scz | 0.985 | Brain_Cerebellar_Hemisphere | abf | `gwas_locus_over_max_snps` | 0.871 | 1e-06 | `robust` | 5/8 |
| POLG | scz | 0.982 | Brain_Cerebellar_Hemisphere | abf | `gwas_locus_over_max_snps` | 0.847 | 1e-06 | `robust` | 1/2 |
| PSMD6 | scz | 0.976 | Brain_Hippocampus | abf | `no_qtl_credible_set` | 0.803 | 1e-06 | `robust` | 1/1 |
| RAI1 | scz | 0.971 | Brain_Cerebellum | abf | `qtl_cs_not_matching_gtex` | 0.768 | 5e-06 | `intermediate` | 0/2 |
| TXNDC15 | als | 0.969 | Brain_Caudate_basal_ganglia | abf | `gwas_no_credible_set` | 0.759 | 5e-06 | `intermediate` | 0/1 |
| SIRPA | ad | 0.965 | Brain_Cerebellar_Hemisphere | abf | `gwas_locus_over_max_snps` | 0.736 | 5e-06 | `intermediate` | 0/13 |
| PRDM2 | als | 0.964 | Brain_Cortex | abf | `gwas_no_credible_set` | 0.729 | 5e-06 | `intermediate` | 0/6 |
| NUP50 | scz | 0.960 | Brain_Caudate_basal_ganglia | abf | `gwas_locus_over_max_snps` | 0.708 | 5e-06 | `intermediate` | 0/11 |
| RASA1 | scz | 0.946 | Brain_Cerebellum | abf | `gwas_no_credible_set` | 0.636 | 5e-06 | `intermediate` | 0/1 |
| PPIL2 | scz | 0.943 | Brain_Amygdala | abf | `gwas_locus_over_max_snps` | 0.624 | 5e-06 | `intermediate` | 0/8 |
| PTPRN | als | 0.918 | Brain_Nucleus_accumbens_basal_ganglia | abf | `qtl_cs_not_matching_gtex` | 0.529 | 5e-06 | `intermediate` | 0/5 |
| SYT5 | scz | 0.904 | Brain_Frontal_Cortex_BA9 | abf | `gwas_locus_over_max_snps` | 0.485 | 5e-06 | `intermediate` | 0/1 |
| ZNF232 | ad | 0.899 | Brain_Spinal_cord_cervical_c-1 | abf | `gwas_locus_over_max_snps` | 0.471 | 5e-06 | `intermediate` | 0/1 |
| WIPI2 | als | 0.893 | Brain_Nucleus_accumbens_basal_ganglia | abf | `gwas_no_credible_set` | 0.454 | 5e-06 | `intermediate` | 0/1 |
| YWHAB | scz | 0.891 | Brain_Frontal_Cortex_BA9 | abf | `gwas_no_credible_set` | 0.451 | 5e-06 | `intermediate` | 0/9 |
| DDRGK1 | pd | 0.884 | Brain_Anterior_cingulate_cortex_BA24 | abf | `gwas_no_credible_set` | 0.433 | 1e-05 | `primary_prior` | 0/1 |
| LPCAT4 | scz | 0.884 | Brain_Putamen_basal_ganglia | abf | `no_qtl_credible_set` | 0.432 | 1e-05 | `primary_prior` | 0/1 |
| NDUFS3 | ad | 0.864 | Brain_Cerebellar_Hemisphere | abf | `no_qtl_credible_set` | 0.388 | 1e-05 | `primary_prior` | 0/1 |
| RERE | scz | 0.856 | Brain_Nucleus_accumbens_basal_ganglia | abf | `no_qtl_credible_set` | 0.374 | 1e-05 | `primary_prior` | 0/2 |
| SH3GL2 | pd | 0.834 | Brain_Substantia_nigra | abf | `gwas_no_credible_set` | 0.335 | 1e-05 | `primary_prior` | 0/9 |
| COPA | scz | 0.830 | Brain_Cerebellum | abf | `gwas_locus_over_max_snps` | 0.327 | 1e-05 | `primary_prior` | 0/1 |
| TAOK2 | scz | 0.829 | Brain_Cerebellum | abf | `qtl_cs_not_matching_gtex` | 0.326 | 1e-05 | `primary_prior` | 0/1 |
| NEK4 | scz | 0.820 | Brain_Frontal_Cortex_BA9 | abf | `gwas_locus_over_max_snps` | 0.313 | 1e-05 | `primary_prior` | 0/1 |
| GLYCTK | scz | 0.812 | Brain_Anterior_cingulate_cortex_BA24 | abf | `gwas_locus_over_max_snps` | 0.301 | 1e-05 | `primary_prior` | 0/1 |
| MED19 | scz | 0.810 | Brain_Hippocampus | abf | `gwas_no_credible_set` | 0.299 | 1e-05 | `primary_prior` | 0/1 |
| NCOR1 | pd | 0.805 | Brain_Nucleus_accumbens_basal_ganglia | abf | `qtl_cs_not_matching_gtex` | 0.293 | 1e-05 | `primary_prior` | 0/1 |

## Loci with a curated event

| gene | trait | PP4 sQTL | tissues coloc | colocalizing introns | interval match | promotes | tier |
|---|---|---|---|---|---|---|---|
| UNC13A | als | 0.959 | 2/13 | chr19:17630750-17632782 | no | no | `disease_locus_splice_linked` |

## Top nominations by posterior

| gene | trait | PP4 sQTL | PP4 eQTL | introns named | driver tx | tier |
|---|---|---|---|---|---|---|
| ACTR1B | scz | 0.997 | 0.996 | 1 | 2 | `novel_splice_led_candidate` |
| ACTR1B | scz | 0.997 | 0.996 | 1 | 2 | `novel_splice_led_candidate` |
| MRPS33 | scz | 0.994 | 0.970 | 1 | 1 | `novel_splice_led_candidate` |
| IRF3 | scz | 0.994 | 0.977 | 1 | 3 | `novel_splice_led_candidate` |
| YPEL1 | scz | 0.985 | 0.864 | 1 | 3 | `novel_splice_led_candidate` |
| POLG | scz | 0.982 | 0.274 | 1 | 0 | `novel_splice_led_candidate` |
| PSMD6 | scz | 0.976 | 0.228 | 1 | 1 | `novel_splice_led_candidate` |
| SNAP91 | scz | 0.974 | 0.930 | 1 | 10 | `novel_splice_led_candidate` |
| RAI1 | scz | 0.971 | 0.717 | 1 | 0 | `novel_splice_led_candidate` |
| TXNDC15 | als | 0.969 | 0.379 | 1 | 1 | `novel_splice_led_candidate` |
| SIRPA | ad | 0.965 | 0.947 | 2 | 4 | `novel_splice_led_candidate` |
| PRDM2 | als | 0.964 | 0.495 | 1 | 5 | `novel_splice_led_candidate` |
| NUP50 | scz | 0.960 | 0.972 | 2 | 9 | `novel_splice_led_candidate` |
| UNC13A | als | 0.959 | 0.377 | 1 | 4 | `disease_locus_splice_linked` |
| DGKZ | scz | 0.953 | 0.858 | 1 | 2 | `novel_splice_led_candidate` |
| CDIP1 | scz | 0.951 | 0.848 | 2 | 7 | `novel_splice_led_candidate` |
| RASA1 | scz | 0.946 | 0.856 | 1 | 0 | `novel_splice_led_candidate` |
| PPIL2 | scz | 0.943 | 0.938 | 1 | 16 | `novel_splice_led_candidate` |
| PTPRN | als | 0.918 | 0.964 | 2 | 4 | `novel_splice_led_candidate` |
| GGNBP2 | als | 0.912 | 0.988 | 1 | 3 | `novel_splice_led_candidate` |
| NDUFAF7 | scz | 0.909 | 0.980 | 1 | 1 | `novel_splice_led_candidate` |
| FGFR1 | scz | 0.908 | 0.653 | 2 | 21 | `novel_splice_led_candidate` |
| SYT5 | scz | 0.904 | 0.911 | 1 | 1 | `novel_splice_led_candidate` |
| ZNF232 | ad | 0.899 | 0.750 | 1 | 6 | `novel_splice_led_candidate` |
| WIPI2 | als | 0.893 | 0.478 | 1 | 3 | `novel_splice_led_candidate` |
| TTC19 | pd | 0.891 | 0.794 | 1 | 2 | `novel_splice_led_candidate` |
| YWHAB | scz | 0.891 | 0.851 | 1 | 2 | `novel_splice_led_candidate` |
| KLC1 | scz | 0.890 | 0.033 | 1 | 4 | `novel_splice_led_candidate` |
| DDRGK1 | pd | 0.884 | 0.902 | 1 | 2 | `novel_splice_led_candidate` |
| LPCAT4 | scz | 0.884 | 0.999 | 1 | 1 | `novel_splice_led_candidate` |

## References for the curated events

Every interval above is traceable to a publication and to the local file it was read from. Nothing in the registry is written from recall.

**UNC13A** (als) — `reviewed`
- PMID [35197626](https://pubmed.ncbi.nlm.nih.gov/35197626/) · DOI [10.1038/s41586-022-04424-7](https://doi.org/10.1038/s41586-022-04424-7)
  - Ma et al., Nature 2022. "TDP-43 represses cryptic exon inclusion in the FTD-ALS gene UNC13A." Abstract: "The top variants associated with FTD or ALS risk in humans are located in the intron harbouring the cryptic exon."
- PMID [35197628](https://pubmed.ncbi.nlm.nih.gov/35197628/) · DOI [10.1038/s41586-022-04436-3](https://doi.org/10.1038/s41586-022-04436-3)
  - Brown et al., Nature 2022. "TDP-43 loss and ALS-risk SNPs drive mis-splicing and depletion of UNC13A." Full text: "Using two minigenes containing exon 20, intron 20 and exon 21, with and without the two ALS- and FTLD-linked variants, we determined that the risk variants enhanced CE upon TDP-43 loss."
- *Coordinate verification:* Derived locally, not quoted. rs12608932 is at chr19:17,641,880 in the project's own GRCh38 TOPMed panel (genotypes/qtl/all_samples/chr19.pvar). UNC13A (ENSG00000130477, minus strand) exon structure was read from the GENCODE v47 GTF cache; the longest protein-coding transcript, ENST00000551649.5 (45 exons), places that variant inside the 20th intron counted 5'->3', spanning chr19:17,641,557-17,642,844 (1,288 bp). That independently reproduces the papers' "intron 20, between exons 20 and 21" from annotation alone.
- *GTEx testability:* GTEx v11 DOES carry an sQTL intron phenotype at this exact intron, and its bounds (17,641,556-17,642,845) reproduce the GENCODE-derived interval to within the 1 bp that LeafCutter and GENCODE differ by at each end -- an independent confirmation of the coordinates above. Two things follow. First, the event is testable in GTEx, so a non-match is a real result rather than a missing phenotype. Second, it is NOT GTEx's grouped-permutation representative intron for UNC13A, so the representative arm never tested it: only the all-introns arm (COLOC_SIGNAL_SQTL=all) can. Re-derived across all 13 brain tissues by `curated_event_testability` on 2026-09-10: the phenotype exists in exactly 2 -- Cerebellar Hemisphere (clu_28403) and Cerebellum (clu_28681) -- and GTEx's own SuSiE calls NO credible set for it in either. Both facts are the expected signature of a TDP-43-loss- dependent cryptic exon in healthy donors: LeafCutter only forms a phenotype where the junction is detected at all, and a barely-spliced junction yields no fine-mappable sQTL. Decisively for the audit, those 2 tissues are exactly the 2 where the upstream UNC13A intron DOES colocalize with ALS, so the curated event was on trial in the very tissues carrying the signal -- which is what makes its silence informative rather than a coverage gap.
- *Scope:* The interval is the cryptic-exon-HARBOURING INTRON, not the cryptic exon itself. The CE is a sub-interval whose exact bounds are in the papers' supplementary data and were not verified here. Intron resolution is the right unit anyway: GTEx's sQTL phenotype is a LeafCutter intron-excision ratio, so the strongest available question is whether the colocalizing intron IS this intron.

**PICALM** (ad) — `reviewed`
- PMID [30297968](https://pubmed.ncbi.nlm.nih.gov/30297968/) · DOI [10.1038/s41588-018-0238-1](https://doi.org/10.1038/s41588-018-0238-1)
  - Raj et al., Nat Genet 2018. "Integrative transcriptome analyses of the aging brain implicate altered splicing in Alzheimer's disease susceptibility." Abstract: "We report that altered splicing is the mechanism for the effects of the PICALM, CLU and PTK2B susceptibility alleles."
  - Supplementary: Supplementary Table 11 (TWAS), file 41588_2018_238_MOESM10_ESM.xlsx, md5 c29d708ec9ece1c93eaaca7e7b3e5f13. PICALM row: Intron.Usage.Cluster "11:85737409:85742511:clu_7404", BEST.GWAS.ID rs3851179, BEST.GWAS.Z -7.91, TWAS.Z 3.52. Local copy: inputs/raw/literature/raj2018/.
- *Which event, and why:* The locus carries more than one PICALM splicing signal and they are NOT the same event, so picking the wrong one would have silently mis-stated the result. Supplementary Table 10 (the FDR<0.05 sQTL catalogue, MOESM9, md5 ff1b3ebf75e59345a0db0169ab4be949) lists exactly one PICALM sQTL -- "11:85689136:85692172:clu_7402", lead SNP rs540422 -- which lifts to chr11:85,978,093-85,981,129 (GRCh38). That is a real PICALM sQTL but it is NOT the event carrying the AD association. The TWAS table anchors PICALM's splicing association on rs3851179, the canonical PICALM AD GWAS lead (GWAS Z = -7.91), and on a DIFFERENT cluster, clu_7404. The reviewed interval above is clu_7404.
- *Coordinate verification:* Derived locally, not quoted. The paper imputed against HRC r1.1 and reports GRCh37 coordinates, so the interval was lifted rather than copied. 187 variants flanking the cluster carry an rsID in both the hg19 1000G EUR panel and this project's GRCh38 TOPMed panel; their offsets are +288,958 (107) and +288,957 (80), a 1 bp ambiguity. Both candidates resolve onto the SAME annotated feature: GENCODE v47 carries a PICALM intron at chr11:86,026,368-86,031,468 (5,101 bp), which fixes the interval independently of the offset. As a cross-check, rs3851179 itself lifts 85,868,640 -> 86,157,598.
- *GTEx testability:* GTEx v11 carries an sQTL phenotype at this exact intron, with LeafCutter bounds 86,026,367-86,031,469 -- the same 1 bp offset from GENCODE seen at UNC13A. So the event is testable and a non-match is a real result. It is NOT GTEx's grouped-permutation representative intron for PICALM, so the representative arm never tested it; only the all-introns arm can. Re-derived across all 13 brain tissues by `curated_event_testability` on 2026-09-10: present in 10 of 13, INCLUDING Cortex, the one tissue where PICALM's coloc.abf posterior reached 0.818. GTEx's own SuSiE calls no credible set for it in any of them; the only PICALM sQTL credible sets in brain are at chr11:85,974,812-85,981,129 (Cortex) and chr11:85,976,682-85,981,129 (Cerebellar Hemisphere), ~50 kb upstream and a different event. SIGNAL-LEVEL STATUS OF PICALM (updated 2026-09-11). On the uniform primary grid (MAX_SNPS = 12,000) PICALM still has no coloc.susie row in either sQTL arm: its AD locus (locus60_chr11, 15,713 SNPs) exceeds the guard, so its primary evidence is coloc.abf alone (Cortex PP4_sQTL 0.818) and the primary audit tiers it `disease_locus_splice_linked`. The prespecified AD sensitivity arm at MAX_SNPS = 30,000 (coloc_signal_susie/sensitivity/max_snps_30000/) fits the locus: coloc.susie gives Cortex PP4_sQTL 0.813 on the same chr11:85,974,812-85,981,129 intron, the curated clu_7404 intron is tested and does not colocalize, and the sensitivity audit tiers PICALM `context_distinct_splice_colocalization`. The locus LD audit (26.locus_ld_robustness) finds that fit LD-robust -- one 2-variant credible set (rs10792832/rs3851179, purity 0.994), stable under 100-500 kb boundary trims, kriging outlier removal and a residual-variance refit -- but prior-sensitive: PP4 0.813 at p12 = 1e-5, 0.685 at 5e-6, 0.303 at 1e-6. Quote the signal-level tier as the sensitivity arm's, never as the primary grid's.
- *Caution:* PICALM is an established AD locus with evidence that the susceptibility allele acts through splicing, but it is not an established causal-splicing positive control in the way UNC13A's cryptic exon is: the TWAS Z is 3.52, and some isoform-associated variants at this locus are not independently associated with AD after conditioning on the main signal. Decisively: the anchor is a **TWAS** association (TWAS.Z 3.52), not a colocalization. TWAS does not establish a shared causal variant -- LD-driven co-regulation at the locus yields the same signal -- so this record is capped at `disease_locus_splice_linked` by `evidence_class` regardless of how well the coordinates match. Promoting it would require a colocalization or a functional demonstration of the event, neither of which this source provides.

Article metadata retrieved from PubMed.

