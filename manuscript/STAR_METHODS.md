# STAR★Methods

*Drafted 2026-09-20 to the Cell Press STAR Methods schema recorded in `MANUSCRIPT_PLAN.md` §7
(indexed 2026-07-20; cell.com returns HTTP 403 to direct fetch, so the live pages must be
eyeballed at final formatting). Section order is the required one. Every Key Resources Table
entry is referenced in Method Details. Versions were read from the execution environments on
2026-09-20, not from documentation.*

---

## KEY RESOURCES TABLE

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
|---|---|---|
| **Deposited data** | | |
| BrainSEQ Phase II/III postmortem brain RNA-seq and genotypes (caudate, DLPFC, hippocampus; SCZD caudate) | Lieber Institute for Brain Development | dbGaP: *accession to be confirmed by the lead contact before submission* |
| GTEx v11 brain RNA-seq, cis-eQTL and cis-sQTL summary statistics (13 brain regions) | GTEx Consortium | GTEx Portal v11; RRID:SCR_013042 |
| ONT long-read RNA-seq, aged DLPFC BA9/46 (n = 12; 6 AD, 6 control), Bambu quantifications | Aguzzoli-Heberle et al., Nat Biotechnol 2024 | SRA: PRJNA1008058 (SRP456327, SRR25723935–SRR25723946); Synapse: syn52047893; Zenodo: 8180677 |
| Schizophrenia GWAS (PGC3, European) | Psychiatric Genomics Consortium wave 3 | PGC3_SCZ_wave3.european.autosome.public.v3 |
| Alzheimer's disease GWAS | Bellenguez et al. 2022 | GWAS Catalog: GCST90027158 (PMID 35379992) |
| Parkinson's disease GWAS | Nalls et al. | GWAS Catalog: GCST009325 |
| Dementia with Lewy bodies GWAS | Chia et al. 2021 | GWAS Catalog: GCST90001390 (PMID 33589841) |
| Amyotrophic lateral sclerosis GWAS | van Rheenen et al. 2021 | GWAS Catalog: GCST90027164 (PMID 34873335) |
| GENCODE v47 primary-assembly annotation and transcript FASTA | GENCODE | GENCODE release 47; RRID:SCR_014966 |
| gnomAD v2.1.1 gene constraint (LOEUF, missense o/e) | Karczewski et al. 2020 | https://doi.org/10.1038/s41586-020-2308-7; RRID:SCR_014964 |
| ClinVar variant classifications | NCBI ClinVar | RRID:SCR_006169 |
| ENCODE eCLIP peaks (HepG2, K562) | ENCODE Consortium | RRID:SCR_015482; per-experiment accessions in `07_rbp_regulation/_m/neuronal_clip_manifests/` |
| Mouse neuronal CLIP (NOVA2 CTAG, NOVA2 cKO, replicate signal) | GEO | GEO: GSE103315, GSE103314, GSE206650 |
| ATtRACT RBP motif database | Giudice et al. | RRID:SCR_016984 |
| 1000 Genomes Phase 3 EUR LD reference (baselineLD, LDSC) | Price lab / LDSC resources | baselineLD v2.2; https://doi.org/10.1038/ng.3954 |
| **Software and algorithms** | | |
| IsoGraph (isoform-switch network inference) | This paper | v0.1.5; Zenodo DOI *to be minted before acceptance* |
| Analysis code for this study | This paper | GitHub repository + Zenodo DOI *to be minted before acceptance* |
| Python | Python Software Foundation | 3.12.13 |
| numpy / pandas / scipy / scikit-learn | — | 2.4.4 / 2.3.3 / 1.17.1 / 1.8.0 |
| PyTorch (VAE fitting) | Paszke et al. | 2.11.0+cu126 |
| python-igraph / leidenalg (module detection) | Csárdi & Nepusz; Traag et al. | 1.0.0 / 0.11.0 |
| statsmodels / pyarrow / matplotlib | — | 0.14.6 / 24.0.0 / 3.10.9 |
| R | R Core Team | 4.4.3 (QTL/coloc environment); 4.5.3 (figures) |
| coloc | Giambartolomei et al. 2014 | 5.2.3; https://doi.org/10.1371/journal.pgen.1004383 |
| susieR | Wang et al. 2020 | 0.14.2; https://doi.org/10.1111/rssb.12388 |
| eCAVIAR | Hormozdiari et al. 2016 | https://doi.org/10.1016/j.ajhg.2016.10.003 |
| stratified LD score regression (S-LDSC) | Finucane et al. 2015 | https://doi.org/10.1038/ng.3404 |
| WGCNA (matched-feature and classical baselines) | Langfelder & Horvath | 1.73; RRID:SCR_003302 |
| satuRn (held-out DTU evidence) | Gilis et al. | 1.18.0 |
| SMR / HEIDI | Zhu et al. 2016 | 1.4.2; RRID:SCR_018627 |
| MAGMA | de Leeuw et al. | v1.10; RRID:SCR_005757 |
| PLINK 2 | Chang et al. | v2.0.0-a.6.9LM (29 Jan 2025); RRID:SCR_001757 |
| bcftools / samtools | Danecek et al. | bcftools 1.23; RRID:SCR_005227 |
| phASER (haplotypic expression) | Castel et al. | outputs mirrored per sample; RRID:SCR_017038 |
| Bambu (long-read quantification) | Chen et al. | as released with the source dataset |
| ggplot2 / patchwork / arrow (R) | — | 4.0.3 / 1.3.2 / 24.0.0 |
| **Other** | | |
| Analysis inventory (one row per analysis, with run status and provenance) | This paper | `reports/pi/_evidence/inventory.parquet` |
| Per-analysis provenance blocks (git commit, input digests, thresholds) | This paper | e.g. `08_integration/_m/anchored_gene_summary/provenance.json` |

---

## RESOURCE AVAILABILITY

### Lead contact

Further information and requests for resources should be directed to and will be fulfilled by
the lead contact, **[NAME, EMAIL — to be completed before submission]**.

### Materials availability

This study did not generate new unique reagents. It is computational; all primary data are
previously published or externally held.

### Data and code availability

- **Previously published data.** All primary molecular data analysed here are previously
  published and are listed in the Key Resources Table with their accessions. BrainSEQ
  genotypes and TOPMed-imputed data are **controlled access** and are available through the
  relevant data-access committee; the authors cannot redistribute them. GTEx v11 summary
  statistics, GWAS summary statistics, GENCODE, gnomAD, ClinVar, ENCODE and ATtRACT are
  publicly available at the identifiers listed.
- **Derived data.** All derived results needed to reproduce the figures and tables are
  deposited: lean per-analysis result tables are in the code repository under git-LFS, and
  four heavy artifacts (~7.8 GB; per-region edge and reconstruction tables) are deposited at
  Zenodo — file list in `zenodo/MANIFEST.tsv`. **DOI to be minted before acceptance.**
- **Code.** All original code — the IsoGraph package (v0.1.5) and the analysis repository,
  including every SLURM wrapper — is public and archived with a **Zenodo DOI to be minted
  before acceptance**, and will be public by the publication date.
- **Provenance.** Each analysis writes a provenance record (git commit, dirty state, package
  versions, thresholds, and per-input size/mtime/SHA-256) so a reader can verify that a table
  matches the inputs it claims. Filesystem mtimes under `_m/` are a checkout artifact and are
  not run dates; SLURM logs carry the run times.
- Any additional information required to reanalyse the data reported here is available from
  the lead contact on request.

---

## EXPERIMENTAL MODEL AND SUBJECT DETAILS

**Human postmortem brain — BrainSEQ.** Bulk RNA-seq and genotypes from the Lieber Institute
for Brain Development. Four artifact stores are used: caudate, DLPFC and hippocampus (aging
analyses) and a schizophrenia case/control caudate cohort (`caudate_sczd`). Donor counts by
region in the allele-aware arm, after QC: caudate 486, DLPFC 498, hippocampus 451 samples with
a BAM, its index and a phASER VCF. The cohort is approximately half African-American (EA
194–248 and AA 213–226 donors per region in the ancestry-split analyses), which is handled
explicitly wherever a European-ancestry GWAS is involved: the two ancestries are fitted
separately and never pooled for orientation.

**Human postmortem brain — GTEx v11.** Thirteen brain regions, used both for module discovery
and for cis-eQTL/cis-sQTL colocalization. GTEx eQTL power exceeds sQTL power in every tissue,
which is why per-gene modality comparisons are interpreted against tissue counts rather than
posteriors alone.

**Human postmortem DLPFC — ONT long read.** Twelve donors (6 AD, 6 control), BA9/46, used as
an orthogonal platform for switch confirmation. The chemistry is cDNA (SQK-PCS111) and
truncates at the 5′ end, which is recorded as an assay limitation wherever a 5′ event fails to
confirm.

**Ethics.** All human data are previously published and de-identified; each source study
carries its own consent and IRB approvals, cited in the Key Resources Table. No new human
subjects were recruited for this work.

---

## METHOD DETAILS

*This section is assembled from the per-analysis summaries under `<stage>/_m/*.md` and the
stage reports in `reports/pi/`. Each subsection below names the CLI that produced it; every CLI
is in the deposited repository and every wrapper is a committed SLURM script, so each number in
the paper is traceable to one command.*

### Isoform-switch feature construction and network inference (IsoGraph)

Per gene, transcript-level quantifications are reduced to a switch coordinate (the first
principal component of the gene's within-gene isoform proportions) and an abundance
coordinate. Genes enter on the switching transcript filter (production since 2026-09-14).
A variational autoencoder embeds the per-gene multiplex features; a gene–gene similarity graph
is built in the latent space and partitioned with Leiden at **resolution 2.0** (production
since 2026-09-16; resolution 5.0 is retained as a disclosed sensitivity arm). Module detection
is run identically for the matched WGCNA baselines, which receive **the same feature matrix**,
so that any difference is attributable to the inference and not to the representation.
CLIs: `real_data/run_models.py`, `real_data/run_matched_wgcna.py`, `real_data/sweep_leiden.py`.

### Module trust, reproducibility and replication

Split-half stability, per-module trust, within-cohort sign concordance, cross-cohort eigengene
projection and permutation nulls. Cross-cohort module-*pair* replication is reported as a
**negative**: BrainSEQ is quantified with Salmon and GTEx with RSEM, and changing the
quantifier destroys more per-gene switch–age signal than changing the brain region does, so
the pair test confounds quantification with biology. CLIs: `real_data/stability.py`,
`module_trust.py`, `replication.py`, `replication_permutation.py`, `eigengene_projection.py`.

### Module characterization and the DTU-without-DGE bound

GO:BP enrichment, cell-type composition (MuSiC deconvolution), incremental association with
phenotype conditional on abundance, and the GO-invisible gate. The scope bound is stated
explicitly: IsoGraph does not outperform WGCNA on per-module GO or phenotype rates; its
contribution is the abundance-independent switch layer. CLIs:
`real_data/module_enrichment.py`, `celltype_composition.py`, `incremental_association.py`,
`go_invisible_gate.py`, `abundance_structure_separation.py`.

### Genetic anchoring: QTL enrichment, colocalization and heritability

Module-level sQTL-versus-eQTL enrichment contrasts with matched baselines; colocalization at
three estimator tiers (`coloc.susie` > `coloc.abf` > eCAVIAR CLPP) against five GWAS; and
S-LDSC partitioned heritability with the baselineLD model. Colocalization uses a 12,000-SNP
guard per locus (30,000 as a sensitivity arm for AD). SMR/HEIDI is run beneath colocalization
on the same summary statistics as a consistency check rather than as independent evidence.
CLIs: `real_data/qtl_anchoring.py`, `coloc_{prep,summary,isoform_events,direction,meta}.py`,
`coloc_signal_susie.py`, `smr_heidi.py`, `ldsc_annot_prep.py`, `ldsc_summary.py`.

### Switch mechanism and orthogonal confirmation

Structural consequence annotation of switch pairs; ONT long-read confirmation against an
abundance-matched null; ISA/satuRn concordance; short-read junction confirmation over every
concordant gene using the LIBD PSI event catalogue; and preprocessing sensitivity, both with
the partition held fixed and with a full refit per setting. CLIs:
`real_data/switch_consequence.py`, `longread_switch_confirm.py`, `isa_concordance.py`,
`junction_coloc_confirm.py`, `switch_feature_sensitivity.py`, `switch_feature_refit.py`.

### Allele-specific isoform usage

An allele-aware junction recount on WASP-tagged BrainSEQ BAMs assigns fragments to phASER `PW`
haplotypes; a beta-binomial GLMM with a donor random intercept contrasts the two haplotypes
within donors heterozygous at the gene's switch-QTL lead, with donors homozygous at the lead
providing a built-in null. Two orientation arms tie the effect to disease: one transfers the
sign from the switch-QTL lead to the locus GWAS lead through signed LD in the in-sample panel
(gated at |r| ≥ 0.8), and one anchors the contrast **on the GWAS lead itself**, which makes
orientation exact at the cost of power. Anchor genotypes are read from phASER's own per-sample
VCFs, because the TOPMed genotype panel is phased in a different frame (measured agreement
0.51 against phASER, versus 0.9996–1.0000 for phASER's own VCFs); the run verifies this before
fitting. CLIs: `real_data/ase_junction_switch.py`, `ase_junction_allelic.py`,
`ase_risk_orientation.py`, `ase_gwas_lead_arm.py`.

### RBP regulon nomination

ATtRACT motif scanning over mature and intronic sequence, a binomial GLM adjusting for motif
*opportunity* (transcript length, GC, UTR/CDS composition, transcript count), and orthogonal
eCLIP binding evidence. The adjustment reorders rather than thins the nominations, and the
layer is reported as candidate regulons, not regulation. CLIs: `real_data/rbp_scan.py`,
`rbp_scan_intronic.py`, `rbp_regulon.py`, `rbp_binding.py`, `neuronal_clip_*.py`.

### Cross-layer integration

Per-gene joins over all layers (`gene_deep_dive.py`), the anchored-gene evidence table
(`anchored_gene_summary.py`), out-of-cohort projection of age-sensitive modules into the SCZ
case/control cohort (`scz_age_projection.py`), functional preservation of cross-cohort matched
modules against a size-matched null (`replication_functional.py`) and the perturbation target
panel (`rbp_target_panel.py`, `rbp_pair_assayability.py`).

---

## QUANTIFICATION AND STATISTICAL ANALYSIS

**Multiple testing.** Benjamini–Hochberg within each pre-specified family, with the family
stated at the point of use. Families are not pooled across stages. Where a correction family
was redefined (for example the EA-only allelic refits, which are a different family from the
pooled Quest fits), the new correction is applied and the change is recorded.

**Nulls.** Every module-similarity and module-preservation comparison uses a **size-matched**
permutation null, because all similarity measures grow with module size. Long-read switch
confirmation uses an abundance-matched null. Cross-cohort replication uses a Freedman–Lane
style permutation.

**Colocalization.** `coloc.abf` PP4 ≥ 0.8 is the nomination threshold, with the prior
p12 = 1×10⁻⁵ pre-specified and a sweep (1×10⁻⁶, 5×10⁻⁶) reported as robustness. eCAVIAR CLPP
is retained as a trailing, secondary column and does not enter any call: it is a per-variant
posterior product, noisy at a single locus, and in the anchored set its two largest values are
contradicted by every other layer.

**Allelic tests.** Beta-binomial GLMM with a donor random intercept; units are donor ×
haplotype. A gene must have ≥ 5 donors informative on both haplotypes to be fitted. Fits at
the optimiser bound are quasi-separation: their sign is informative and their magnitude is
not, so they are counted apart and never reported as effect sizes. A model-free robust score
test is reported beside the GLMM.

**Module projections.** An allelic effect is signed along a module's switch axis only when the
gene is in the module through switching (not abundance) **and** at least one transcript of the
pair correlates with the module score at q ≤ 0.05; otherwise the projection is withheld with a
recorded reason. This was added after a gene whose module membership was abundance-driven
produced four significant pairs whose projected signs disagreed while the underlying cis
effect was a single coherent one.

**Orthogonal confirmation.** One independent assay counts as confirmation; genes confirmed by
two or three are marked. The three assays fail for unrelated reasons (platform, tissue
coverage, catalogue coverage), so a gene reachable by only one is less-sampled rather than
weaker. Minor-form usage is reported as a descriptor and never used as a gate: these are
neurotypical tissues, and an isoform that matters in disease is expected to be a minority form
here. A usage range whose minimum is exactly zero is reported as **indeterminate** — read
depth and genuine absence are not separable — and is read in neither direction.

**Seeds and determinism.** Random seeds are recorded per analysis in the corresponding
`*_stats.json` or `params.json`. The per-gene joins and the anchored summary are deterministic.
Module partitions are re-assigned at every fit, so every module-id join is gated by a partition
fingerprint; a stale join is refused rather than silently mis-mapped.

**Software.** Versions as listed in the Key Resources Table, read from the execution
environments on 2026-09-20.
