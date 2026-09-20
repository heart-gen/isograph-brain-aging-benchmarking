## RBP-regulon enrichment of co-switch modules

### Purpose
Test the trans-regulatory hypothesis behind IsoGraph co-switch modules. The coding-consequence analysis showed the switch axis is dominated by 3′UTR/CDS remodeling — usage that is heavily controlled by sequence-specific RNA-binding proteins (RBPs). This motivates the coordination hypothesis: a co-switch module is coordinated because its member genes share a common trans-acting RBP, the one mechanism that cis-sQTL anchoring cannot explain. For each module and RBP we ask whether RBP binding-site *switching* (a motif gained or lost between the two switch-pair isoforms) is over-represented among the module's genes relative to the region's switch-gene background.

The central methodological risk is **motif opportunity**: long, GC-rich or UTR-heavy transcripts present more sequence for a degenerate PWM to match, so a module of such genes can look like a regulon without any shared regulation. Every module × RBP cell is therefore reported under two arms — an unadjusted hypergeometric and a covariate-adjusted binomial GLM — and the two are allowed to disagree in the written result rather than being reconciled in favour of the larger number.

### Inputs
- Switch-pair isoform sequences: GENCODE v47 mature transcript FASTA (`inputs/raw/gencode_v47/gencode.v47.transcripts.fa.gz`), restricted to transcripts appearing in any region's `structure_switch_pairs.parquet`.
- Genome FASTA + GTF (GENCODE v47, shared PSC staging) for the intronic splice-site-flank scope.
- RBP motif models: ATtRACT human position-weight matrices (`inputs/rbp_motifs/ATtRACT_db.txt` + `pwm.txt`, *Homo sapiens* only), mapped Matrix_id → Gene_name; 160 RBPs carry at least one matrix.
- Motif-similarity families (`rbp_motif_families.parquet`, 136 families over 1,194 matrices) for the deduplicated `--unit family_id` arm.
- Per-gene opportunity covariates (`rbp_scan_opportunity.parquet`, `transcript_regions.parquet`): transcript length, GC, 5′UTR/CDS/3′UTR base composition, transcript count.
- Module / switch-gene assignments and GO-invisible tags per region: 14 regions (BrainSEQ + GTEx).
- Intermediates: `rbp_counts.parquet` (mature), `rbp_counts_intronic.parquet` (intronic flanks), `rbp_family_counts.parquet`.

### Methods Text
The analysis is two-stage. Stage 1 (dedicated `motif` conda env) collected every transcript in any region's IsoGraph switch pairs, extracted its mature sequence from the GENCODE v47 transcript FASTA, and scanned it against the ATtRACT human RBP PWMs with MOODS. Because ATtRACT motifs are RNA and single-stranded, matrices were treated with the U column as T and only the sense strand was scanned; PWMs were converted to log-odds with a 0.1 pseudocount and a per-position match p-value threshold of 1×10⁻⁴ defined a motif hit. Rather than a flat 0.25 background, which systematically inflates AU-rich binders (ELAVL, CPEB, hnRNPD), each transcript was assigned to one of 20 GC × 2 purine composition bins and one threshold set was solved per bin against a background estimated from that bin's sequence; flat-background counts were retained in parallel under `bg_mode='flat'` (preserved as the `*_flatbg` tables). Hits were tallied per RBP, per motif-similarity family, and partitioned into 5′UTR/CDS/3′UTR/noncoding classes by projecting GTF CDS spans into transcript coordinates. The same pipeline was run over intronic splice-site flanks derived from the GTF (pre-mRNA sense, minus strand reverse-complemented) to cover the binding niche of splicing-regulatory RBPs invisible to a mature-transcript scan.

Stage 2 (isograph env) made a within-pair site-switch call: for each gene and RBP the switch was scored as altering that RBP's binding when the motif was present in exactly one of the two switch-pair isoforms (a presence gain/loss, which controls for transcript length by comparing two isoforms of the same gene). Then, per region and module, two tests were applied to every module × RBP cell with at least 3 module genes. First, a **hypergeometric** test asking whether the number of the module's genes with an RBP site-switch exceeded expectation given the region-wide switch-gene pool, Benjamini–Hochberg corrected across the full cell set. Second, a **covariate-adjusted binomial GLM** regressing the per-gene site-switch indicator on module membership with log transcript length, GC, 5′UTR/CDS/3′UTR fractions and log transcript count as covariates, reported as an adjusted odds ratio with a 95% Wald interval. Cells whose (in-module × switched) 2×2 has an empty margin admit no finite maximum-likelihood estimate; these were detected before fitting by a zero-cell criterion and after fitting by pinned fitted probabilities or degenerate coefficients, excluded from the model, and reported as an explicit estimability breakdown, so BH is applied over exactly the set of tests actually performed. Each module's GO-invisible status was carried through for stratification.

### Results Text

*Numbers below are the switching-filter re-run at canonical Leiden resolution 2.0 (2026-09-19).
The motif-family arm is the exception and is flagged where it appears: its table has not been
regenerated since 2026-08-29 and still describes the legacy production.*

Across 16 regions, 11,790 switch genes were scored against 160 ATtRACT human RBP models, yielding 11,040 module × RBP cells. On the unadjusted hypergeometric, **384 module–RBP pairs were significant at BH q < 0.05 (61 of them in GO-invisible modules), spanning 125 distinct RBPs and 36 distinct region × module regulons**. Of the 11,040 cells, 9,955 (90.2%) were estimable under the adjusted GLM; the remainder were dropped for complete separation (510 zero-cell), no variation to model (530), or quasi-separation (45). Resolution 2.0 yields far fewer, larger modules than the legacy res-5.0 production, so the cell count fell roughly threefold and every count below is on that smaller test set.

**The opportunity adjustment removes most of this signal, and the survivors are largely not the same cells.** Only **89 of the 384** hypergeometric hits also reach adjusted q < 0.05 (66 RBPs, 15 region × module regulons, 15 GO-invisible); 216 have an adjusted CI entirely above 1, and **1 has an adjusted CI entirely below 1** — enriched on raw counts yet significantly *depleted* once length, GC and UTR composition are accounted for. The adjusted arm is not simply a subset: it calls 416 cells, most of which the raw test does not. Raw enrichment and adjusted odds ratio are only weakly concordant among the significant cells (Spearman ρ = 0.35, n = 383), so the two arms rank regulons differently rather than agreeing at different thresholds.

The reordering is sharpest at the top of the adjusted arm, and module identifiers cannot be
carried over from the legacy run: Leiden re-assigns ids at every fit, so the previously quoted
regulons (GTEx frontal cortex "M008" and its ELAVL4/KHDRBS1 hits) name different gene sets now
and have been dropped rather than re-quoted. On the re-run the strongest adjusted results are
**depletions** in GTEx frontal cortex BA9 M005, a GO-invisible module: CELF5 (OR 0.54, 95% CI
0.45–0.65, q = 8.8×10⁻⁷), HNRNPLL (0.54, 0.45–0.65) and CELF4 (0.57, 0.47–0.69) — fewer
site-switches than length, GC and UTR composition predict. The strongest adjusted *enrichments*
are in putamen M001: ACO1 (1.60, 1.37–1.85), SRSF10 (1.61, 1.37–1.90) and PABPC1 (1.63,
1.38–1.92). Cross-region recurrence also differs between the arms: PPRC1 (8 of 16 regions) and
IGF2BP3 (7) lead the raw test, while RBFOX2 (6) and IGF2BP3 (5) lead the adjusted one.

**The scope arms behave the same way.** Extending to intronic splice-site flanks — the binding niche of splicing regulators, and the scope the neuronal-CLIP arm builds on — gives **488** raw hits (77 GO-invisible) and **188** adjusted-significant. The motif-family arm (136 motif-similarity families instead of individual RBPs) has **not been regenerated on the switching filter**: `rbp_regulon_family.parquet` still dates from 2026-08-29, so its legacy figures (555 raw, 36 adjusted, 97.2% estimable) are not comparable with the counts above and are not quoted as current. Re-run it before it appears in the manuscript.

Orthogonal ENCODE eCLIP evidence is consistent with binding *capacity* but does not rescue factor specificity. It also comes from the wrong tissue: ENCODE profiles its ~168 RBPs in **HepG2 and K562, not brain**, so every eCLIP statement below is about whether these factors *can* occupy these intervals in a cell line, not whether they *do* in postmortem cortex or caudate. Read it as a capacity floor, never as brain occupancy: 21 of 35 testable regulon RBPs are binding-supported (preferential switched-interval binding at BH q ≤ 0.05; 31 of 35 show the preference before correction), representing 47 of the 60 motif families among those testable RBPs, with a median switched−constitutive bound-rate gap of only 0.020 — a small, near-universal alternative-exon skew rather than selective occupancy by the predicted factors.

The defensible claim is therefore narrower than the unadjusted counts suggest: co-switch modules show widespread RBP binding-site switching, but after opportunity adjustment only a small minority of module–RBP pairs are supported, and those are led by factors (NONO, RBM8A, CPEB4, ZC3H10, NOVA2) that raw enrichment ranking would have missed. This remains an **exploratory, hypothesis-generating** result: motif-presence is a computational prediction of altered RBP binding, not measured binding, and the adjusted-significant set is the appropriate basis for experimental follow-up.

### Figure and Table Notes
- Potential supplementary table: `07_rbp_regulation/_m/rbp/rbp_regulon.parquet` (rendered `07_rbp_regulation/_m/rbp/RBP_REGULON.md`)
  - Rationale: full per-region module × RBP enrichment under both arms, supporting the candidate-regulon claim and its qualification.
  - Key columns: `region`, `module_id`, `rbp`, `module_size`, `n_switched`, `pool_switched`, `enrichment`, `q` (hypergeometric); `odds_ratio`, `ci_low`, `ci_high`, `pvalue_glm`, `qvalue_glm`, `model_status` (adjusted); `go_invisible`.
- Parallel scopes/units: `rbp_regulon_family.parquet`, `rbp_regulon_intronic.parquet`, `rbp_regulon_combined.parquet` (+ matching `RBP_REGULON_*.md`); `*_flatbg` variants preserve the superseded flat-background results.
- Potential supplementary table: `07_rbp_regulation/_m/rbp/rbp_switch_calls.parquet` (per gene × RBP gain/loss calls).
- Potential main/supplementary figure: raw enrichment versus adjusted odds ratio across the 384 significant cells, with the single contradicted cell marked — this makes the opportunity confound visible in one panel and motivates the adjusted arm. A recurrence panel (RBP × number of regions) should be drawn from the adjusted arm, not the raw one. Neither is yet generated.

### Reproducibility Information
- Analysis directory: `07_rbp_regulation/_m/rbp/`.
- Primary scripts: `isograph_benchmark/real_data/rbp_scan.py` and `rbp_scan_intronic.py` (stage 1, motif env), `isograph_benchmark/real_data/rbp_regulon.py` (stage 2, isograph env), `rbp_motif_families.py` (family definitions).
- Execution commands:
  - `sbatch 07_rbp_regulation/_h/02.rbp_regulon.sh --stage regulon` (mature, per-RBP; job 44484238)
  - `sbatch 07_rbp_regulation/_h/02.rbp_regulon.sh --stage regulon --unit family_id` (job 44484239)
  - `sbatch 07_rbp_regulation/_h/03.rbp_regulon_intronic.sh` (intronic scan + intronic and combined scopes; job 44484240)
  - `--stage regulon` re-tests the frozen stage-1 count tables without repeating the MOODS scan, which is the expensive frozen input.
- Parameters: motif hit p-threshold 1×10⁻⁴, pseudocount 0.1, **composition-matched background over 20 GC × 2 purine bins** (flat background retained in parallel), sense-strand only, U→T; hypergeometric over-representation with BH across module × RBP cells; adjusted binomial GLM with covariates `log_length`, `gc`, `frac_5utr`, `frac_cds`, `frac_3utr`, `log_n_transcripts`, BH across estimable cells only; estimability guards `_SEP_TOL` 1×10⁻⁸, `_MAX_SE` 100, `_MAX_ABS_COEF` 30; minimum module size 3 switch genes; switch-gene FDR default 0.05.
- Output files (2026-08-26): `rbp_counts_intronic.parquet` (stage 1, 00:40); `rbp_switch_calls*.parquet`, `rbp_regulon*.parquet`, `RBP_REGULON*.md` (stage 2/3, 00:54–01:38).
- Execution date: 2026-08-26 (supersedes the 2026-07-18 flat-background run preserved as `*_flatbg`).
- Git commit: `012b1a3`; the estimability gate and two-arm report were uncommitted at run time.
- Compute environment: PSC Bridges-2 RM-shared, account bio260021p, 8 cpus-per-task (16 GB).
  - Stage 1 env `/ocean/projects/bio260021p/shared/opt/envs/motif`: Python 3.11.15, MOODS-python 1.9.4.1, pandas 3.0.3, pyfaidx.
  - Stage 2 env `/ocean/projects/bio260021p/shared/opt/envs/isograph`: Python 3.12.13, numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6.
- Determinism: the mature per-RBP run was executed twice (jobs 44452849 and 44484238) and reproduced every reported count exactly. No sampling is involved; no seed applies.
- Missing reproducibility information: per-run package versions are not echoed into SLURM logs; recorded here from the runtime conda envs.

### Limitations and Integration Notes
The adjusted arm is the load-bearing result and it is deliberately conservative: the covariates (length, GC, UTR composition, transcript count) are themselves correlated with the biology, so an RBP that genuinely acts on long 3′UTRs will be partly adjusted away. The 43 surviving pairs are best read as a lower bound on real regulons, and the 19 contradicted pairs as a firm exclusion; the 344 cells with an adjusted CI above 1 but not surviving multiplicity are the ambiguous middle and should not be counted as findings. Motif presence is a sequence prediction, not measured binding, and ATtRACT PWMs are short and degenerate, so individual calls are noisy — cross-region recurrence on the adjusted arm is more trustworthy than any single module–RBP pair. The eCLIP layer establishes binding capacity at alternative exons generally, not factor-specific occupancy, and it is measured in HepG2/K562 rather than brain — a cell-type mismatch that bounds every binding claim here, since RBP expression, competing factors and isoform availability all differ between a transformed cell line and postmortem neural tissue. It therefore cannot substitute for perturbation evidence in a neural context; `rbp_target_panel.py` derives the ranked perturbation candidate set from these tables.

This analysis is the trans-regulatory complement to the two cis analyses: coding-consequence (`06_switch_mechanism/_m/SWITCH_CONSEQUENCE_SUMMARY.md`) established the 3′UTR-remodeling substrate, and sQTL/eQTL colocalization (`05_genetic_anchoring/_m/coloc/`) established cis genetic anchoring — the shared-RBP result offers a candidate trans mechanism for *why* the module's genes co-switch. RBP motifs: ATtRACT database [@doi:10.1093/database/baw035]; scanning with MOODS [@doi:10.1093/bioinformatics/btp554]. [citation needed: GENCODE v47]. [citation needed: ENCODE eCLIP].
