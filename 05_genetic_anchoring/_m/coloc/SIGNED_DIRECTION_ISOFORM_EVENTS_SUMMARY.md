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
Across the five anchoring analyses, 72 of 73 colocalized gene–phenotype events were signed
for risk-allele direction. Of 55 splicing-QTL events, 42 mapped to at least one GENCODE
transcript junction, and **9 splicing junctions were concordant with the tissue-matched
IsoGraph switching transcript pair; all 9 lay in GO-invisible switch modules**. Concordant
aging events included *SNCA* (α-synuclein), where the risk allele increased usage of junction
chr4:89,835,692–89,836,127 in both LBD (risk allele A) and PD (risk allele C) in cortex /
frontal cortex; *CTSH* and *TPCN1* (AD); and *PGS1*, *PPP6R2*, and *GGNBP2* (ALS), with
*PPP6R2* showing two competing junctions moving in opposite directions — a canonical splice
switch. Six events replicated the switch in an independent BrainSeq cohort, including *MYO18A*
(schizophrenia; five junctions in the SCZ-diagnosis caudate cohort) and *MRPS10* (AD, caudate).
Schizophrenia produced no tissue-matched GTEx concordance but five BrainSeq-replicated events,
consistent with its switch layer being defined natively in the BrainSeq caudate cohort. This
is a secondary/interpretive analysis built on the primary colocalization result; concordance
reflects mapping of a genetically anchored splice change onto an IsoGraph switch, not a formal
statistical test.

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
- Candidate main-figure vignette: *SNCA* signed splice event concordant in both LBD and PD —
  the strongest single illustration that the switch layer is genetically anchored to a known
  neurodegeneration gene; pair with the coloc CLPP panel.

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
non-concordant, so the 9 concordant events are a conservative floor. The IsoGraph `::switch`
feature is a gene-level composite over the gene's transcripts; transcript-pair identity is
supplied by the structural-annotation layer, not the VAE. This analysis should integrate with
(i) the eCAVIAR colocalization summary (which it signs and resolves), (ii) the S-LDSC
partitioned-heritability result (which establishes size-robust genetic anchoring), and (iii)
the GO-invisible gate (every concordant aging event is GO-invisible, reinforcing that the
anchored biology is missed by abundance/pathway enrichment). The aging-switch-layer → SCZ
projection (`aging__scz`) is pending and will add one row to the rollup.
