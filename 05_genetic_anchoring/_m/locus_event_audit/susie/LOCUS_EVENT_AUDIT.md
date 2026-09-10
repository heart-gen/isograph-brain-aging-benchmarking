# Locus event audit

Nomination source: `susie` at `PP4_sQTL >= 0.8`.

A colocalization posterior names a gene, not an event. This audit walks GWAS signal -> sQTL signal -> intron phenotype -> driver transcript, and asks whether the colocalizing intron matches a curated literature event.

## Curated events

| gene | trait | status | evidence class | interval | can promote |
|---|---|---|---|---|---|
| UNC13A | als | **reviewed** | `functional_validation` | chr19:17,641,557-17,642,844 | yes |
| PICALM | ad | **reviewed** | `twas_association` | chr11:86,026,368-86,031,468 | **no** |

Promotion to `known_mechanism_recovered` needs BOTH gates. `status: reviewed` means the coordinates were independently verified here against GENCODE v47. `evidence_class` must be `functional_validation` or `coloc_association`: a TWAS association does **not** establish a shared causal variant -- LD-driven co-regulation yields the same signal -- so a TWAS-only anchor caps a locus at `disease_locus_splice_linked` however well its coordinates match. Where that cap binds, `curated_interval_matched` still records whether the interval lined up.

## Tier counts

- `known_mechanism_recovered`: 0
- `disease_locus_splice_linked`: 2
- `novel_splice_led_candidate`: 41
- `not_resolved`: 0

## Loci with a curated event

| gene | trait | PP4 sQTL | tissues coloc | colocalizing introns | interval match | promotes | tier |
|---|---|---|---|---|---|---|---|
| UNC13A | als | 0.961 | 2/13 | chr19:17630750-17632782 | no | no | `disease_locus_splice_linked` |
| PICALM | ad | 0.818 | 1/13 | chr11:85974812-85981129 | no | no | `disease_locus_splice_linked` |

## Top nominations by posterior

| gene | trait | PP4 sQTL | PP4 eQTL | introns named | driver tx | tier |
|---|---|---|---|---|---|---|
| TPP1 | als | 0.996 | 0.774 | 1 | 0 | `novel_splice_led_candidate` |
| GPM6A | scz | 0.991 | 0.990 | 1 | 3 | `novel_splice_led_candidate` |
| CTSH | ad | 0.989 | 0.993 | 2 | 44 | `novel_splice_led_candidate` |
| PITPNM2 | pd | 0.989 | 0.614 | 1 | 2 | `novel_splice_led_candidate` |
| RNASEH2C | scz | 0.981 | 0.975 | 1 | 1 | `novel_splice_led_candidate` |
| PGS1 | als | 0.976 | 0.165 | 1 | 5 | `novel_splice_led_candidate` |
| PLCB2 | scz | 0.975 | 0.157 | 1 | 1 | `novel_splice_led_candidate` |
| SNCA | lbd | 0.974 | 0.113 | 1 | 2 | `novel_splice_led_candidate` |
| SNCA | pd | 0.970 | 0.382 | 1 | 2 | `novel_splice_led_candidate` |
| TXNDC15 | als | 0.969 | 0.379 | 1 | 1 | `novel_splice_led_candidate` |
| DOC2A | scz | 0.966 | 0.755 | 2 | 6 | `novel_splice_led_candidate` |
| SIRPA | ad | 0.965 | 0.947 | 2 | 4 | `novel_splice_led_candidate` |
| TMED4 | scz | 0.965 | 0.954 | 2 | 2 | `novel_splice_led_candidate` |
| UNC13A | als | 0.961 | 0.377 | 1 | 4 | `disease_locus_splice_linked` |
| NT5C2 | scz | 0.958 | 0.963 | 1 | 6 | `novel_splice_led_candidate` |
| CDIP1 | scz | 0.949 | 0.854 | 2 | 7 | `novel_splice_led_candidate` |
| PPIL2 | scz | 0.943 | 0.938 | 1 | 16 | `novel_splice_led_candidate` |
| EDEM3 | scz | 0.926 | 0.996 | 1 | 11 | `novel_splice_led_candidate` |
| PPIP5K1 | scz | 0.926 | 0.761 | 2 | 6 | `novel_splice_led_candidate` |
| SPI1 | ad | 0.926 | 0.881 | 1 | 1 | `novel_splice_led_candidate` |
| PTPRN | als | 0.918 | 0.964 | 2 | 4 | `novel_splice_led_candidate` |
| GGNBP2 | als | 0.911 | 0.989 | 1 | 3 | `novel_splice_led_candidate` |
| SYT5 | scz | 0.904 | 0.911 | 1 | 1 | `novel_splice_led_candidate` |
| ZNF232 | ad | 0.899 | 0.750 | 1 | 6 | `novel_splice_led_candidate` |
| WIPI2 | als | 0.893 | 0.478 | 1 | 3 | `novel_splice_led_candidate` |
| TTC19 | pd | 0.891 | 0.794 | 1 | 2 | `novel_splice_led_candidate` |
| KLC1 | scz | 0.882 | 0.033 | 1 | 4 | `novel_splice_led_candidate` |
| DOC2A | ad | 0.881 | 0.986 | 1 | 6 | `novel_splice_led_candidate` |
| PAK6 | scz | 0.866 | 0.996 | 1 | 7 | `novel_splice_led_candidate` |
| NDUFS3 | ad | 0.864 | 0.099 | 1 | 1 | `novel_splice_led_candidate` |

## References for the curated events

Every interval above is traceable to a publication and to the local file it was read from. Nothing in the registry is written from recall.

**UNC13A** (als) — `reviewed`
- PMID [35197626](https://pubmed.ncbi.nlm.nih.gov/35197626/) · DOI [10.1038/s41586-022-04424-7](https://doi.org/10.1038/s41586-022-04424-7)
  - Ma et al., Nature 2022. "TDP-43 represses cryptic exon inclusion in the FTD-ALS gene UNC13A." Abstract: "The top variants associated with FTD or ALS risk in humans are located in the intron harbouring the cryptic exon."
- PMID [35197628](https://pubmed.ncbi.nlm.nih.gov/35197628/) · DOI [10.1038/s41586-022-04436-3](https://doi.org/10.1038/s41586-022-04436-3)
  - Brown et al., Nature 2022. "TDP-43 loss and ALS-risk SNPs drive mis-splicing and depletion of UNC13A." Full text: "Using two minigenes containing exon 20, intron 20 and exon 21, with and without the two ALS- and FTLD-linked variants, we determined that the risk variants enhanced CE upon TDP-43 loss."
- *Coordinate verification:* Derived locally, not quoted. rs12608932 is at chr19:17,641,880 in the project's own GRCh38 TOPMed panel (genotypes/qtl/all_samples/chr19.pvar). UNC13A (ENSG00000130477, minus strand) exon structure was read from the GENCODE v47 GTF cache; the longest protein-coding transcript, ENST00000551649.5 (45 exons), places that variant inside the 20th intron counted 5'->3', spanning chr19:17,641,557-17,642,844 (1,288 bp). That independently reproduces the papers' "intron 20, between exons 20 and 21" from annotation alone.
- *GTEx testability:* GTEx v11 DOES carry an sQTL intron phenotype at this exact intron, and its bounds (17,641,556-17,642,845) reproduce the GENCODE-derived interval to within the 1 bp that LeafCutter and GENCODE differ by at each end -- an independent confirmation of the coordinates above. Two things follow. First, the event is testable in GTEx, so a non-match is a real result rather than a missing phenotype. Second, it is NOT GTEx's grouped-permutation representative intron for UNC13A, so the representative arm never tested it: only the all-introns arm (COLOC_SIGNAL_SQTL=all) can. Of the brain tissues checked (Cerebellar Hemisphere, Cortex, Frontal Cortex BA9, Spinal cord) the phenotype is present ONLY in Cerebellar Hemisphere, which is consistent with a TDP-43-loss-dependent cryptic exon being close to absent in normal tissue.
- *Scope:* The interval is the cryptic-exon-HARBOURING INTRON, not the cryptic exon itself. The CE is a sub-interval whose exact bounds are in the papers' supplementary data and were not verified here. Intron resolution is the right unit anyway: GTEx's sQTL phenotype is a LeafCutter intron-excision ratio, so the strongest available question is whether the colocalizing intron IS this intron.

**PICALM** (ad) — `reviewed`
- PMID [30297968](https://pubmed.ncbi.nlm.nih.gov/30297968/) · DOI [10.1038/s41588-018-0238-1](https://doi.org/10.1038/s41588-018-0238-1)
  - Raj et al., Nat Genet 2018. "Integrative transcriptome analyses of the aging brain implicate altered splicing in Alzheimer's disease susceptibility." Abstract: "We report that altered splicing is the mechanism for the effects of the PICALM, CLU and PTK2B susceptibility alleles."
  - Supplementary: Supplementary Table 11 (TWAS), file 41588_2018_238_MOESM10_ESM.xlsx, md5 c29d708ec9ece1c93eaaca7e7b3e5f13. PICALM row: Intron.Usage.Cluster "11:85737409:85742511:clu_7404", BEST.GWAS.ID rs3851179, BEST.GWAS.Z -7.91, TWAS.Z 3.52. Local copy: inputs/raw/literature/raj2018/.
- *Which event, and why:* The locus carries more than one PICALM splicing signal and they are NOT the same event, so picking the wrong one would have silently mis-stated the result. Supplementary Table 10 (the FDR<0.05 sQTL catalogue, MOESM9, md5 ff1b3ebf75e59345a0db0169ab4be949) lists exactly one PICALM sQTL -- "11:85689136:85692172:clu_7402", lead SNP rs540422 -- which lifts to chr11:85,978,093-85,981,129 (GRCh38). That is a real PICALM sQTL but it is NOT the event carrying the AD association. The TWAS table anchors PICALM's splicing association on rs3851179, the canonical PICALM AD GWAS lead (GWAS Z = -7.91), and on a DIFFERENT cluster, clu_7404. The reviewed interval above is clu_7404.
- *Coordinate verification:* Derived locally, not quoted. The paper imputed against HRC r1.1 and reports GRCh37 coordinates, so the interval was lifted rather than copied. 187 variants flanking the cluster carry an rsID in both the hg19 1000G EUR panel and this project's GRCh38 TOPMed panel; their offsets are +288,958 (107) and +288,957 (80), a 1 bp ambiguity. Both candidates resolve onto the SAME annotated feature: GENCODE v47 carries a PICALM intron at chr11:86,026,368-86,031,468 (5,101 bp), which fixes the interval independently of the offset. As a cross-check, rs3851179 itself lifts 85,868,640 -> 86,157,598.
- *GTEx testability:* GTEx v11 carries an sQTL phenotype at this exact intron in 5 of the 6 brain tissues checked (all but Anterior cingulate BA24), with LeafCutter bounds 86,026,367-86,031,469 -- the same 1 bp offset from GENCODE seen at UNC13A. So the event is testable and a non-match is a real result. It is NOT GTEx's grouped-permutation representative intron for PICALM, so the representative arm never tested it; only the all-introns arm can.
- *Caution:* PICALM is an established AD locus with evidence that the susceptibility allele acts through splicing, but it is not an established causal-splicing positive control in the way UNC13A's cryptic exon is: the TWAS Z is 3.52, and some isoform-associated variants at this locus are not independently associated with AD after conditioning on the main signal. Decisively: the anchor is a **TWAS** association (TWAS.Z 3.52), not a colocalization. TWAS does not establish a shared causal variant -- LD-driven co-regulation at the locus yields the same signal -- so this record is capped at `disease_locus_splice_linked` by `evidence_class` regardless of how well the coordinates match. Promoting it would require a colocalization or a functional demonstration of the event, neither of which this source provides.

Article metadata retrieved from PubMed.

