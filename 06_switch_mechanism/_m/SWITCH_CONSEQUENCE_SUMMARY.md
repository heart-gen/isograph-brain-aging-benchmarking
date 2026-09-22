## Switch coding-consequence and NMD enrichment

### Purpose
Characterize *what kind* of transcript remodeling the IsoGraph switch axis encodes. For every gene with a phenotype-significant isoform switch, we ask whether the switched isoform pair differs by specific structural/coding consequences (UTR change, CDS change, first/last/internal-exon change, biotype switch, coding-status change, NMD status) more often than expected under a within-gene null. This distinguishes a *productive* remodeling layer (UTR/CDS tuning of a still-coding transcript) from a *degradative* one (routing to nonsense-mediated decay or loss of coding potential). It also tests whether the GO-invisible switch modules — the ones missed by gene-abundance/pathway analysis — carry the same consequence signature as GO-visible ones.

### Inputs
- Per-region IsoGraph switch structure for the phenotype-significant switch genes: `02_module_discovery/{brainseq,gtex}/<region>/_m/isograph_vae/module_interpret/` (`structure_switch_pairs.parquet`, `structure_annotations.parquet`, `transcript_polarity_table.parquet`), with GO-invisible / GO-visible module tags from the module-trust layer.
- GENCODE v47 transcript models for the pairwise structural comparison (`isograph.explain.structure.compare_pair`), including CDS/exon coordinates and strand for the strand-aware 50-nt NMD rule.
- Ten of the 17 array regions had ≥1 phenotype-significant switch gene (FDR < 0.05) and were analyzed: **nine aging** analyses (BrainSeq caudate, hippocampus; GTEx amygdala, anterior_cingulate_cortex_ba24, cerebellum, cortex, frontal_cortex_ba9, hippocampus, hypothalamus) and the **SCZD** caudate analysis. The remaining 7 (BrainSeq dlpfc; GTEx caudate_basal_ganglia, cerebellar_hemisphere, nucleus_accumbens_basal_ganglia, putamen_basal_ganglia, spinal_cord_cervical_c_1, substantia_nigra) were skipped with the reason "no phenotype-significant switch genes" and contribute no output. (Region list verified 2026-09-21 against the per-region `consequence_enrichment.parquet` files; the previous list, which named substantia_nigra instead of BrainSeq hippocampus, was from a superseded run.)
- Stratum coverage among the nine aging analyses: GO-invisible k = 6 (BrainSeq hippocampus; GTEx cerebellum, cortex, frontal_cortex_ba9, hippocampus, hypothalamus); GO-visible k = 8 (BrainSeq caudate; GTEx amygdala, anterior_cingulate_cortex_ba24, cerebellum, cortex, frontal_cortex_ba9, hippocampus, hypothalamus).

### Methods Text
For each region, every phenotype-significant switch gene was decomposed into its IsoGraph switch-pair isoforms, and each pair was annotated for nine mutually informative consequence classes using a GENCODE v47 pairwise comparison: first-exon change, last-exon change, internal-exon difference, UTR change, CDS change, biotype switch, coding-status change, a combined coding consequence (CDS or coding-status change), and an NMD-status switch. NMD status was assigned per isoform with the strand-aware 50-nt rule (a premature termination codon more than 50 nt upstream of the final exon–exon junction), and `nmd_switch` was set when the two isoforms of a pair differed in NMD status. To remove the transcript-count and transcript-length confounds that differ across genes, enrichment was assessed against a *within-gene* permutation null: for each gene we drew the same number of random transcript pairs from that gene's own transcript pool and recomputed each consequence rate, repeating 1,000 times (seed 13). Per region and consequence class we report the observed rate, the null mean, their ratio (enrichment), a one-sided empirical permutation p-value and a gene-level block-bootstrap standard error of the log enrichment (2,000 resamples), stratified into all switch genes, GO-invisible modules, and GO-visible modules. Aging and disease (SCZD) analyses were pooled separately. Within each class, per-analysis log enrichments were combined per consequence class and stratum by a DerSimonian–Laird random-effects model, reporting the pooled enrichment ratio, its 95% CI, and I². The number of regions enriched or depleted at empirical p < 0.05, the median enrichment, and a Fisher combination of the empirical p-values (floored at 1×10⁻⁶ to respect permutation resolution) were retained as secondary summaries.

### Results Text
*(Re-quoted 2026-09-19 against `switch_consequence_meta.parquet` from the switching-filter
re-run at Leiden 2.0. This file is hand-written and does not regenerate with the analysis;
the aging and disease arms are pooled **separately** here, which the legacy version did not
do — it pooled 10 regions across both and is superseded.)*

Across the **9 aging** analyses the switch axis is consistently a **productive
UTR/CDS-remodeling** layer, not a decay-routing one. UTR change is enriched in **9/9**
regions (pooled RE ratio **1.247**, 95% CI 1.22–1.28, p = 5.1e-81, I² = 0.83; observed in
50% of switches; Fisher-combined p = 6.3e-18), and CDS change — equivalently the combined
coding consequence — is enriched in **8/9** (pooled **1.040**, 1.04–1.04, p = 1.4e-126,
I² = 0.00; observed rate 78%; Fisher p = 3.6e-16). Every consequence associated with loss of
a productive coding transcript is **never** enriched in any region: NMD-status switch
(pooled 0.933, observed in only 17% of switches, 0/9 enriched), coding-status change (0.884,
0/9) and biotype switch (0.896, 0/9), all with Fisher-combined p = 1.00. First-exon,
last-exon and internal-exon differences are near-neutral (pooled 0.958–0.994, 0/9 enriched).

The **disease** (SCZD caudate) arm gives the same signature on its single analysis: UTR
change 1.186 (p = 2.6e-31, observed 39%) and CDS change 1.040 (p = 1.2e-12, observed 71%)
enriched; NMD switch 0.895, coding-status change 0.958 and biotype switch 0.941 all
depleted. Because it is one analysis, it is reported beside the aging pool, never merged
into it.

Critically, the GO-invisible modules carry the *same* signature as GO-visible ones — UTR
change pooled 1.266 vs 1.242, CDS 1.033 vs 1.043, NMD depleted in both (0.915 vs 0.940) —
showing that the switch biology missed by pathway analysis is structurally identical
productive remodeling rather than an artifact class. The result is
exploratory-confirmatory: a structural characterization of the switch layer, internally
validated by the within-gene null and reproducible across nine aging analyses, but it does
not by itself establish downstream functional impact.

**A note on the pooled p-values.** The fixed-effect intervals here are extremely tight
because the per-analysis denominators are large; the effect sizes (1.04 for CDS, 1.25 for
UTR) are modest and are what should be quoted. Read I² alongside them: the UTR result is
heterogeneous (0.83) while the CDS result is homogeneous (0.00).

### Figure and Table Notes
- Potential supplementary table: `06_switch_mechanism/_m/switch_consequence_meta.parquet` (rendered `06_switch_mechanism/_m/SWITCH_CONSEQUENCE_META.md`)
  - Rationale: full per-class × stratum × consequence rollup (random-effects pooled ratio and CI, I², regions enriched/depleted, median enrichment, Fisher p) supporting the "productive remodeling, not NMD" claim.
  - Key columns: `analysis_class`, `stratum`, `consequence`, `n_regions`, pooled ratio / 95% CI / p / I² (from `stats.meta_analysis.meta`), `n_enriched_p05`, `n_depleted_p05`, `median_enrichment`, `median_obs_rate`, `fisher_p`.
- Potential main-figure panel: a diverging bar / dot panel of median enrichment per consequence class (UTR and CDS above 1; NMD, coding-status, biotype below 1), faceted or overlaid GO-invisible vs GO-visible, to make the "switch = productive UTR/CDS remodeling" message a single visual. Underlying per-region points come from the per-region `switch_consequence/consequence_enrichment.parquet` files. No such figure has been generated yet.

### Reproducibility Information
- Analysis directory: `06_switch_mechanism/_m/` (rollup: `switch_consequence_meta.parquet`, `SWITCH_CONSEQUENCE_META.md`); per-region outputs under `02_module_discovery/{brainseq,gtex}/<region>/_m/isograph_vae/switch_consequence/`.
- Primary scripts: `isograph_benchmark/real_data/switch_consequence.py` (per-region permutation test and bootstrap SE), `isograph_benchmark/real_data/switch_consequence_meta.py` (cross-region random-effects rollup, aging and disease pooled separately; Fisher p secondary).
- Execution command: `sbatch 06_switch_mechanism/_h/01a.switch_consequence.sh` (SLURM array 0–16, one region per task, tree-qualified region names), then `06_switch_mechanism/_h/02a.switch_consequence_meta.sh`.
- Parameters: `--n-perm 1000`, `--n-boot 2000`, `--seed 13`, switch-gene FDR default 0.05; NMD 50-nt rule, ±2 bp tolerance in structural matching.
- Execution date: 2026-07-18 (array logs `02_module_discovery/brainseq/_m/logs/switch-consequence-42399857_*.log`, all 17 tasks reached "Complete"; meta refreshed 12:48).
- Git commit: baseline `86a3ef4`; analysis code committed in the same session (uncommitted at run time).
- Compute environment: PSC Bridges-2 RM-shared, account bio260021p; conda env `/ocean/projects/bio260021p/shared/opt/envs/isograph`.
- Key package versions (from env at analysis time): Python 3.12.13, numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6.
- Random seed: 13.
- Missing reproducibility information: per-run package versions are not echoed into the SLURM logs; they are recorded here from the runtime conda env rather than a per-run `sessionInfo`-equivalent.

### Limitations and Integration Notes
The within-gene null controls transcript count and length but not sequence conservation or expression level, so enrichment reflects relative consequence composition, not absolute functional burden. "UTR change" and "CDS change" are structural annotations, not measured effects on stability or translation — the RBP-regulon analysis (`07_rbp_regulation/_m/rbp/`) is the mechanistic follow-up that asks *which* trans-factors could read out the enriched 3′UTR remodeling, and the clinical-consequence analysis (gnomAD constraint / ClinVar density over switched exons) is the outstanding orthogonal test of functional impact. This summary should sit alongside the signed-direction colocalization / isoform-event summary (`05_genetic_anchoring/_m/coloc/SIGNED_DIRECTION_ISOFORM_EVENTS_SUMMARY.md`): together they argue that IsoGraph's genetically-anchored, GO-invisible switches are productive UTR/CDS-remodeling events rather than decay artifacts. [citation needed: GENCODE v47]; NMD 50-nt rule [citation needed: NMD 50-nt boundary rule].
