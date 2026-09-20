# Genetics and orthogonal validation for the anchored switch genes

**30** genes (33 gene x trait rows) whose colocalizing sQTL junction resolves to the tissue-matched IsoGraph switch pair.

Sorted on **coloc**, which is the primary genetic evidence. eCAVIAR **CLPP is secondary** and is carried as a trailing column only: it is noisy per locus and, for the two largest values in this set, contradicted by every other layer. It does not enter `genetics_call`.

Generated 2026-09-20T04:34:26+00:00 from commit `35424460` **with a dirty working tree -- not reproducible from that commit alone**, pandas 2.3.3, numpy 2.4.4, over 12 input files. Digests, sizes and mtimes are in `provenance.json`; re-run `python -m isograph_benchmark.real_data.anchored_gene_summary` to regenerate. Do not edit any number here by hand.

## Genetics

`PP4` are coloc.abf posteriors over GTEx brain; `tis` counts tissues reaching a coloc call, which matters because GTEx eQTL power exceeds sQTL power everywhere. `SMR` is supported/tested instrumented probes, `!n` marks HEIDI rejections. A gene with sQTL support and no eQTL instrument at all is the cleanest splicing-specific case.

| gene | trait | tissue | PP4 sQTL | PP4 eQTL | tis s/e | SMR sQTL | SMR eQTL | GO-inv | call | CLPP |
| --- | --- | --- | ---: | ---: | :---: | :---: | :---: | :---: | --- | ---: |
| DOC2A | scz | Cerebellar_Hemisphere | 0.966 | 0.755 | 9/0 | 5/8 !3 | -- | yes | splicing-specific | 0.048 |
| PRDM2 | als | Cortex | 0.964 | 0.494 | 6/0 | 7/7 | -- | no | splicing-specific | 0.056 |
| THAP3 | scz | Cerebellar_Hemisphere | 0.964 | 0.219 | 2/0 | 1/1 | -- | no | splicing-specific | 0.023 |
| IFNAR2 | ad | Cerebellum | 0.870 | 0.829 | 7/1 | 3/15 !2 | -- | no | splicing-specific | 0.029 |
| TPP1 | als | Cerebellum | 0.996 | 0.774 | 1/0 | 2/2 | 1/1 | yes | splicing-led | 0.480 |
| SNAP91 | scz | Cerebellum | 0.974 | 0.930 | 6/2 | 6/10 | 0/1 | yes | splicing-led | 0.047 |
| NT5C2 | scz | Cerebellum | 0.959 | 0.963 | 7/2 | 4/4 | 0/3 | no | splicing-led | 0.017 |
| CRELD2 | scz | Cerebellum | 0.949 | 0.852 | 13/2 | 48/51 | 5/5 | yes | splicing-led | 0.062 |
| PTPRN | als | Cerebellum | 0.918 | 0.964 | 5/1 | 7/7 | 1/1 | yes | splicing-led | 0.021 |
| DNAJA3 | scz | Cerebellar_Hemisphere | 0.920 | 0.282 | 3/0 | 0/2 | -- | yes | splicing-led (coloc only) | 0.025 |
| INO80E | scz | Cerebellar_Hemisphere | 0.947 | 0.970 | 4/9 | 3/4 | 1/3 !2 | yes | colocalizes, not splicing-specific | 0.034 |
| DLG1 | scz | Cerebellar_Hemisphere | 0.275 | 0.134 | 0/0 | -- | -- | no | weak coloc | 0.026 |
| VAMP2 | als | Cortex | 0.275 | 0.283 | 0/0 | -- | -- | no | weak coloc | 0.013 |
| TBC1D15 | pd | Frontal_Cortex_BA9 | 0.237 | 0.043 | 0/0 | -- | -- | no | weak coloc | 0.017 |
| CTSB | pd | Amygdala | 0.215 | 0.940 | 0/3 | 2/3 | 1/1 | no | weak coloc | 0.011 |
| GPR135 | scz | Cerebellum | 0.121 | 0.117 | 0/0 | 0/7 | 0/1 | yes | weak coloc | 0.010 |
| FLCN | ad | Cerebellum | 0.096 | 0.014 | 0/0 | -- | -- | yes | weak coloc | 0.010 |
| BAIAP3 | als | Hypothalamus | 0.086 | 0.065 | 0/0 | -- | -- | no | weak coloc | 0.013 |
| CTC1 | als | Cerebellum | 0.083 | 0.758 | 0/0 | -- | -- | no | weak coloc | 0.053 |
| PIGQ | als | Putamen_basal_ganglia | 0.067 | 0.031 | 0/0 | -- | -- | yes | weak coloc | 0.015 |
| RPAIN | scz | Cortex | 0.066 | 0.509 | 0/0 | -- | -- | no | weak coloc | 0.020 |
| FLCN | scz | Cerebellar_Hemisphere | 0.060 | 0.014 | 0/0 | -- | -- | yes | weak coloc | 0.014 |
| SPG7 | scz | Cerebellum | 0.048 | 0.026 | 0/0 | -- | -- | yes | weak coloc | 0.037 |
| PCGF3 | pd | Cerebellum | 0.036 | 0.010 | 0/0 | -- | -- | no | weak coloc | 0.011 |
| NADSYN1 | scz | Cerebellar_Hemisphere | 0.031 | 0.235 | 0/0 | -- | -- | yes | weak coloc | 0.016 |
| TARBP1 | scz | Cerebellum | 0.012 | 0.062 | 0/0 | -- | -- | yes | weak coloc | 0.010 |
| B3GAT1 | scz | Cortex | 0.004 | 0.936 | 0/1 | -- | -- | yes | weak coloc | 0.011 |
| PPP6R2 | scz | Frontal_Cortex_BA9 | 0.001 | 0.092 | 0/0 | -- | -- | no | weak coloc | 0.058 |
| TMEM175 | als | Cerebellar_Hemisphere | 0.000 | 0.139 | 0/0 | 6/36 | 1/3 !1 | yes | weak coloc | 0.483 |
| RBFA | scz | Cerebellum | 0.000 | 0.001 | 0/0 | -- | -- | yes | weak coloc | 0.013 |
| TMEM175 | ad | Cerebellum | 0.000 | 0.137 | 0/0 | -- | -- | yes | weak coloc | 0.035 |
| TMEM175 | pd | Cortex | 0.000 | 0.031 | 0/0 | -- | -- | yes | weak coloc | 0.051 |
| GSTO2 | scz | Frontal_Cortex_BA9 | -- | -- | 0/0 | -- | -- | yes | weak coloc | 0.021 |

## Orthogonal validation

Three assays with different failure modes. **Long-read**: ONT DLPFC BA9/46, n = 12, Bambu -- do the switching transcripts exist and trade off on an independent platform? **Junction recount**: allele-aware junction counts from BrainSEQ BAMs at ~n = 500 -- is the junction itself observed, and in how many donors? Restricted to the switch pairs that actually CARRY the colocalizing junction, so the columns are about that junction rather than about the gene switching somewhere. The junction usually sits in several such pairs and the usage ratio depends on which partner it is measured against, so a RANGE over those pairs is reported rather than one pair standing for the rest. **PSI**: the LIBD event catalogue, the narrowest of the three -- `junction not in catalogue` is a gap in that resource, not a negative.

**Minor-form usage is a descriptor, not a criterion.** BrainSEQ and GTEx are neurotypical tissue, so an isoform that matters in disease is often a minority form here precisely because disease is what raises it. `usage` is reported for interpretation; the call keys on whether both forms are observed and in how many donors (`both` >= 30). For reference only, 0.05 is the pre-registered threshold the separate PSI arm applies.

| gene | trait | LR conf | LR pairs det/switch | ASE pairs test/carrying | ASE donors both | ASE usage range | PSI | switch replicates | call |
| --- | --- | :---: | :---: | :---: | ---: | :---: | --- | :---: | --- |
| DOC2A | scz | no | 4/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | tested, not confirmed |
| PRDM2 | als | yes | 10/4 | 11/12 | 498/498 | 0.034-0.457 | junction not in catalogue | yes | confirmed (2 assays) |
| THAP3 | scz | no | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| IFNAR2 | ad | no | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| TPP1 | als | yes | 1/1 | 0/0 | 0/0 | -- | no BrainSEQ region | no | confirmed (1 assay) |
| SNAP91 | scz | yes | 22/6 | 14/18 | 486/486 | 0.084-0.403 | no BrainSEQ region | no | confirmed (2 assays) |
| NT5C2 | scz | yes | 2/0 | 5/5 | 494/495 | 0.221-0.475 | no BrainSEQ region | no | confirmed (1 assay) |
| CRELD2 | scz | yes | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | tested, not confirmed |
| PTPRN | als | no | 6/0 | 5/14 | 451/451 | 0.046-0.460 | no BrainSEQ region | no | confirmed (1 assay) |
| DNAJA3 | scz | yes | 6/3 | 3/9 | 486/486 | 0.050-0.283 | no BrainSEQ region | no | confirmed (2 assays) |
| INO80E | scz | no | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| DLG1 | scz | yes | 6/0 | 24/34 | 498/498 | 0.142-0.496 | no BrainSEQ region | no | confirmed (1 assay) |
| VAMP2 | als | no | 2/0 | 4/9 | 451/451 | 0.273-0.491 | not validated | no | confirmed (1 assay) |
| TBC1D15 | pd | yes | 2/2 | 3/15 | 486/486 | 0.111-0.409 | junction not in catalogue | no | confirmed (2 assays) |
| CTSB | pd | yes | 4/0 | 14/15 | 451/451 | 0.000-0.373 | no BrainSEQ region | no | confirmed (1 assay) |
| GPR135 | scz | no | 0/0 | 0/1 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| FLCN | ad | yes | 0/0 | 2/13 | 486/486 | 0.226-0.413 | no BrainSEQ region | no | confirmed (1 assay) |
| BAIAP3 | als | no | 2/0 | 8/12 | 480/485 | 0.000-0.479 | no BrainSEQ region | no | confirmed (1 assay) |
| CTC1 | als | no | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| PIGQ | als | no | 0/0 | 4/15 | 486/486 | 0.085-0.455 | validated | no | confirmed (2 assays) |
| RPAIN | scz | yes | 0/0 | 7/12 | 451/451 | 0.088-0.455 | validated | no | confirmed (2 assays) |
| FLCN | scz | yes | 0/0 | 2/13 | 486/486 | 0.226-0.413 | no BrainSEQ region | no | confirmed (1 assay) |
| SPG7 | scz | yes | 0/0 | 7/15 | 498/498 | 0.052-0.455 | no BrainSEQ region | no | confirmed (1 assay) |
| PCGF3 | pd | no | 10/0 | 11/20 | 498/498 | 0.000-0.384 | no BrainSEQ region | no | confirmed (1 assay) |
| NADSYN1 | scz | yes | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | tested, not confirmed |
| TARBP1 | scz | no | 0/0 | 14/15 | 486/486 | 0.000-0.488 | no BrainSEQ region | no | confirmed (1 assay) |
| B3GAT1 | scz | no | 2/0 | 4/10 | 450/451 | 0.002-0.064 | validated | no | confirmed (2 assays) |
| PPP6R2 | scz | yes | 2/2 | 9/15 | 451/451 | 0.038-0.468 | junction not in catalogue | no | confirmed (2 assays) |
| TMEM175 | als | no | 10/2 | 0/0 | 0/0 | -- | no BrainSEQ region | no | tested, not confirmed |
| RBFA | scz | no | 0/0 | 0/0 | 0/0 | -- | no BrainSEQ region | no | not measurable |
| TMEM175 | ad | no | 10/2 | 0/0 | 0/0 | -- | no BrainSEQ region | no | tested, not confirmed |
| TMEM175 | pd | no | 10/2 | 0/0 | 0/0 | -- | validated | no | confirmed (1 assay) |
| GSTO2 | scz | yes | 18/6 | 0/0 | 0/0 | -- | validated | no | confirmed (2 assays) |

## Counts

| genetics call | rows |
| --- | ---: |
| weak coloc | 22 |
| splicing-led | 5 |
| splicing-specific | 4 |
| splicing-led (coloc only) | 1 |
| colocalizes, not splicing-specific | 1 |

| orthogonal call | rows |
| --- | ---: |
| confirmed (1 assay) | 13 |
| confirmed (2 assays) | 9 |
| not measurable | 6 |
| tested, not confirmed | 5 |

## Reproducibility

| input | bytes | mtime (UTC) |
| --- | ---: | --- |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/coloc/coloc_isoform_events_combined.parquet` | 65777 | 2026-09-20T03:50:44+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/coloc_signal_susie/genes.parquet` | 888100 | 2026-09-19T15:21:31+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/05_genetic_anchoring/_m/smr_heidi/gtex/smr_results.parquet` | 288031 | 2026-09-19T15:21:33+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/longread_switch_confirm/gene_confirmation.parquet` | 89537 | 2026-09-19T15:21:33+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/longread_switch_confirm/pair_confirmation.parquet` | 1413387 | 2026-09-19T15:21:33+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/junction_coloc_confirm/junction_confirm.parquet` | 27877 | 2026-09-20T04:10:54+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/caudate/pair_feasibility.parquet` | 2769570 | 2026-09-19T20:49:23+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/caudate/junction_allelic_counts.parquet` | 402553512 | 2026-09-19T20:49:45+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/dlpfc/pair_feasibility.parquet` | 894995 | 2026-09-19T20:49:23+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/dlpfc/junction_allelic_counts.parquet` | 59993539 | 2026-09-19T20:49:27+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/hippocampus/pair_feasibility.parquet` | 3099714 | 2026-09-19T20:49:22+00:00 |
| `/ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking/06_switch_mechanism/_m/ase_junction_switch/hippocampus/junction_allelic_counts.parquet` | 415206482 | 2026-09-19T20:49:52+00:00 |

Thresholds: `PP4_STRONG` = 0.8, `MIN_DONORS_BOTH` = 30, `MIN_PAIR_FRAGS` = 10, `MIN_USAGE_REFERENCE` = 0.05.