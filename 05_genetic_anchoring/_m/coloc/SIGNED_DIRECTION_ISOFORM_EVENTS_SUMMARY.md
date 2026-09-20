## Signed colocalization direction and resolved isoform events

### Purpose
Colocalization (eCAVIAR CLPP) establishes that a GWAS credible set and a GTEx brain QTL
credible set share a causal variant for an IsoGraph switch gene, but is unsigned and does
not name the isoform change. This analysis resolves, per prioritized (colocalized) switch
gene: (i) the **direction** in which the trait-increasing (risk) allele shifts splicing or
expression, and (ii) the **actual isoform event** — the LeafCutter intron junction and the
IsoGraph switching transcript pair it maps onto — with an independent-cohort replication
check. It covers a disease case (schizophrenia) and an aging case (Alzheimer's disease,
Parkinson's disease, Lewy body dementia, amyotrophic lateral sclerosis).

### Inputs
- Colocalized genes per trait: `05_genetic_anchoring/_m/coloc/<layer>__<trait>/coloc/coloc_colocalized_genes.tsv`
  (from the eCAVIAR CLPP capstone), for `brainseq-sczd__scz`, `aging__ad`, `aging__pd`,
  `aging__lbd`, `aging__als`.
- GWAS summary statistics (signed effect + effect allele), via the trait registry
  `isograph_benchmark/real_data/gwas_traits.py`: PGC3 SCZ (EUR) [@doi:10.1038/s41586-022-04434-5],
  Bellenguez AD [@doi:10.1038/s41588-022-01024-z], Nalls PD [@doi:10.1016/S1474-4422(19)30320-5],
  Chia LBD [@doi:10.1038/s41588-021-00785-3], van Rheenen ALS [@doi:10.1038/s41588-021-00973-1].
- GTEx v11 brain `sQTLs`/`eQTLs` `signif_pairs.parquet` (signed cis-QTL `slope`, per ALT
  allele) and SuSiE credible sets [citation needed: GTEx v11 data release].
- IsoGraph structural interpretation per brain region (`interpret_modules.py`): tissue-matched
  GTEx `structure_switch_pairs.parquet`, `transcript_polarity_table.parquet`,
  `structure_annotations.parquet`; independent BrainSeq cohorts (caudate, hippocampus, DLPFC,
  caudate SCZ-diagnosis) for replication.
- GENCODE v47 primary-assembly annotation (per-transcript exon coordinates → intron set)
  [citation needed: GENCODE].

### Methods Text
For each colocalized switch gene we defined the trait-increasing (risk) allele as the GWAS
effect allele when its signed effect was positive and the alternate allele otherwise, and
aligned it to the GTEx effect (ALT) allele of the colocalizing variant to obtain a signed
QTL effect (`coloc_direction.py`). For splicing QTLs, the GTEx phenotype identifier was
parsed into its LeafCutter intron junction [citation needed: LeafCutter]. Each junction was
mapped to the transcript(s) that splice it out by deriving per-transcript intron coordinates
from the GENCODE v47 exon annotation (intron = [exon_i end + 1, exon_{i+1} start − 1]; a
junction was assigned to a transcript when both boundaries matched within 2 bp), restricted
to the gene's own transcripts (`coloc_isoform_events.py`). A junction–switch concordance flag
was set when a junction-containing transcript was a member of the IsoGraph switching
transcript pair reported for that gene **in the same GTEx tissue**. Independent replication
was assessed by testing whether the junction transcript was also an IsoGraph switch-pair
isoform in the matched BrainSeq cohort. The structural nature of each switch was summarized
as the union of exon/CDS/UTR/biotype/coding-status changes across the switching pair isoforms
from `structure_annotations`. All GWAS effect-allele lookups used name-resolved streaming of
the registry sumstats; no per-SNP recomputation was performed.

### Results Text
*(Re-quoted 2026-09-19 against `coloc_direction_combined.parquet` and
`coloc_isoform_events_{combined,meta}.parquet` from the switching-filter re-run at Leiden
2.0. This file is hand-written and does not regenerate with the analysis.)*

Across the six anchoring analyses (five traits; SCZ enters twice, from the aging and the
SCZD switch layers), **all 467** colocalized gene–phenotype events were signed for
risk-allele direction — 379 splicing-QTL and 88 expression-QTL. Of the 379 sQTL events,
**267 mapped to at least one GENCODE transcript junction** and **76 junctions in 30 genes
were concordant with the tissue-matched IsoGraph switching transcript pair**; **49 of the
76 lie in GO-invisible switch modules**. Independent BrainSEQ replication of the switch is
recorded for 21 events over 7 genes (CAMK1, CDIP1, FANCL, GPM6A, NMRK1, PRDM2, TMEM107),
but **only 1 of those 21 is also concordant** (PRDM2, ALS, cortex); the other 20 replicate a
switch whose junction did not map into the tissue-matched switch pair. Concordance and
independent replication are therefore very nearly disjoint, which limits how strongly any
single concordant event can be presented as validated.

Concordance by trait: SCZ 51 (across the aging and SCZD layers), ALS 13 (BAIAP3, CTC1,
PIGQ, PRDM2, PTPRN, TMEM175, TPP1, VAMP2), AD 7 (FLCN, IFNAR2, TMEM175), PD 5 (CTSB,
PCGF3, TBC1D15, TMEM175), LBD 0. TMEM175 appears under AD, ALS and PD, but not as one
junction resolved three times: ALS and PD share chr4:932540-947709(+) (different tissues,
different leads rs873786 / rs77060135), while the AD events are two other junctions
(chr4:951717-952367, chr4:952450-953190) at a third lead, rs11552301. Its eCAVIAR CLPP of
0.483 is also not corroborated by the other genetic layers — coloc.abf PP4_sQTL is ~0 in
every trait (4.3e-12 PD, 1.7e-5 AD, 1.7e-4 ALS), no tissue reaches the sQTL or eQTL coloc
call, and SMR supports only 6 of 36 instrumented sQTL probes with its single eQTL signal
HEIDI-rejected. TMEM175 must not be used as a worked example on the strength of CLPP alone.

**SNCA no longer resolves, and is kept as a falsification example.** The LBD risk allele A
still increases usage of junction chr4:89,835,692–89,836,127 in cortex (risk QTL effect
+0.73), but the junction's transcripts are **not** in the IsoGraph switch pair for that
tissue, so the event is scored non-concordant. On the legacy expression filter it was
reported as one of nine concordant events and was the vignette in Fig 4 panel A. That panel
must not show SNCA as resolved; the PI decision of 2026-09-18 keeps SNCA in the manuscript
as an example of a nomination the method declines to support, not as evidence for it.

This is a secondary/interpretive analysis built on the primary colocalization result;
concordance reflects mapping of a genetically anchored splice change onto an IsoGraph
switch, not a formal statistical test, and the count scales with how many events the
colocalization layer nominates.

### Figure and Table Notes
- Potential supplementary table: `05_genetic_anchoring/_m/coloc/coloc_isoform_events_combined.parquet`
  (+ `coloc_direction_combined.parquet`).
  - Rationale: full per-event resolution (risk allele, signed effect, junction, mapped
    transcript, switch-pair concordance, BrainSeq replication, structural consequence).
  - Key columns: `trait, gene_name, kind, tissue, junction, risk_allele, risk_qtl_effect,
    junction_transcripts, concordant, brainseq_replicates_switch, structural_consequence,
    resolved_event`.
- Potential supplementary table: `05_genetic_anchoring/_m/coloc/coloc_isoform_events_meta.parquet`
  (rendered `COLOC_ISOFORM_EVENTS_META.md`) — per-analysis rollup of concordant / GO-invisible
  / BrainSeq-replicated counts.
- Candidate main-figure vignette: **not SNCA** (it no longer resolves; see above). On the
  re-run the best-supported worked example is *PRDM2* (ALS, cortex): sQTL PP4 0.964 vs eQTL
  0.494, colocalizing in 6 tissues on splicing and 0 on expression, SMR 7/7 supported with no
  HEIDI rejection and no eQTL instrument, and the only concordant event that replicates the
  switch in an independent BrainSEQ cohort. Its weaknesses are a modest CLPP (0.056), a
  GO-visible module, and a junction/switch polarity r of only 0.11.

### Reproducibility Information
- Analysis directory: `05_genetic_anchoring/_m/coloc/`
- Primary scripts: `isograph_benchmark/real_data/coloc_direction.py`,
  `isograph_benchmark/real_data/coloc_isoform_events.py` (both **uncommitted** at time of writing).
- Input files: as listed under Inputs (coloc outputs, GTEx v11 signif_pairs, GENCODE v47
  gtf cache `04_module_characterization/_m/tmp/gencode.v47.primary_assembly.annotation.gtf_cache.parquet`).
- Output files: `coloc_direction.parquet` + `COLOC_DIRECTION.md` per analysis;
  `coloc_isoform_events.parquet` + `COLOC_ISOFORM_EVENTS.md` per analysis; combined parquets
  and `COLOC_ISOFORM_EVENTS_META.md` at the coloc root.
- Execution command: `python -m isograph_benchmark.real_data.coloc_direction` then
  `python -m isograph_benchmark.real_data.coloc_isoform_events`.
- Execution date: 2026-07-18.
- Git commit: repository HEAD 86a3ef4; the two analysis scripts are untracked (not yet committed).
- Compute environment: conda env `/ocean/projects/bio260021p/shared/opt/envs/isograph`;
  Python 3.12.13, pandas 2.3.3, numpy 2.4.4. Run on the interactive node (light join/stream;
  GWAS effect-allele streaming is I/O-bound over the full sumstats).
- Random seed: not applicable (deterministic joins).
- Missing reproducibility information: GTEx v11 and GENCODE v47 release citekeys are
  placeholders; LeafCutter citekey is a placeholder; the upstream coloc `susieR` version is
  recorded in the coloc CLPP job logs, not re-captured here.

### Limitations and Integration Notes
Concordance requires the GTEx sQTL junction to fall within 2 bp of an annotated GENCODE
intron of a named IsoGraph switch-pair isoform; genuine but unannotated or novel junctions,
and switches whose top pair excludes the junction-containing transcript, are counted as
non-concordant, so the 76 concordant events are a conservative floor. Conversely the test is
permissive in a way that cuts the other direction: `n_switch_pairs` is at its cap of 15 for
the median concordant event, so "concordant" means the junction transcript appears in one of
up to 15 reported switch pairs, and the junction-to-switch polarity r is below 0.15 for 24 of
the 30 concordant genes -- membership in the pair list is not evidence that the junction
drives the switch axis. Separately, `structural_consequence` is `no annotated structural
change` for all 76 concordant events on this re-run and empty for the other 391, where the
legacy run carried UTR/CDS remodeling; this looks like a join defect in
`coloc_isoform_events.py` rather than a biological result, and it should be checked before
any figure panel relies on the structural layer. The IsoGraph `::switch`
feature is a gene-level composite over the gene's transcripts; transcript-pair identity is
supplied by the structural-annotation layer, not the VAE. This analysis should integrate with
(i) the eCAVIAR colocalization summary (which it signs and resolves), (ii) the S-LDSC
partitioned-heritability result (which establishes size-robust genetic anchoring), and (iii)
the GO-invisible gate (every concordant aging event is GO-invisible, reinforcing that the
anchored biology is missed by abundance/pathway enrichment). The aging-switch-layer → SCZ
projection (`aging__scz`) is pending and will add one row to the rollup.
