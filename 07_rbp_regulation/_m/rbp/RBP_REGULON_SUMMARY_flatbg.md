## RBP-regulon enrichment of co-switch modules

### Purpose
Test the trans-regulatory hypothesis behind IsoGraph co-switch modules. The coding-consequence analysis showed the switch axis is dominated by 3′UTR/CDS remodeling — usage that is heavily controlled by sequence-specific RNA-binding proteins (RBPs). This motivates the coordination hypothesis: a co-switch module is coordinated because its member genes share a common trans-acting RBP, the one mechanism that cis-sQTL anchoring cannot explain. For each module and RBP we ask whether RBP binding-site *switching* (a motif gained or lost between the two switch-pair isoforms) is over-represented among the module's genes relative to the region's switch-gene background — a significant RBP marks the module as a candidate regulon. Results are stratified GO-invisible vs GO-visible.

### Inputs
- Switch-pair isoform sequences: GENCODE v47 mature transcript FASTA (`inputs/raw/gencode_v47/gencode.v47.transcripts.fa.gz`), restricted to transcripts appearing in any region's `structure_switch_pairs.parquet`.
- RBP motif models: ATtRACT human position-weight matrices (`inputs/rbp_motifs/ATtRACT_db.txt` + `pwm.txt`, *Homo sapiens* only), mapped Matrix_id → Gene_name.
- Module / switch-gene assignments and GO-invisible tags per region (same 10 regions with phenotype-significant switch genes as the coding-consequence analysis).
- Intermediate: per-(transcript, RBP) motif hit counts `07_rbp_regulation/_m/rbp/rbp_counts.parquet`.

### Methods Text
The analysis is two-stage. Stage 1 (dedicated `motif` conda env) collected every transcript in any region's IsoGraph switch pairs, extracted its mature sequence from the GENCODE v47 transcript FASTA, and scanned it against the ATtRACT human RBP PWMs with MOODS. Because ATtRACT motifs are RNA and single-stranded, matrices were treated with the U column as T and only the sense strand was scanned; PWMs were converted to log-odds against a flat background with a 0.1 pseudocount, and a per-position match p-value threshold of 1×10⁻⁴ defined a motif hit. Matrix hits were summed per RBP (an RBP may have several matrices), giving per-(transcript, RBP) counts. Stage 2 (isograph env) made a within-pair site-switch call: for each gene and RBP the switch was scored as altering that RBP's binding when the motif was present in exactly one of the two switch-pair isoforms (a presence gain/loss, which controls for transcript length by comparing two isoforms of the same gene). Then, per region and module, a hypergeometric test asked whether the number of the module's genes with an RBP site-switch exceeded expectation given the region-wide switch-gene pool, over all module × RBP pairs with non-zero counts; p-values were Benjamini–Hochberg corrected across the full module × RBP test set. Modules with fewer than 3 switch genes were dropped. Each module's GO-invisible status was carried through for stratification.

### Results Text
Across the 10 regions, 9,430 switch genes were scored against 157–160 ATtRACT human RBP models, yielding 16,973 module × RBP hypergeometric tests. **829 module–RBP pairs were significant at BH q < 0.05 (245 of them in GO-invisible modules), spanning 129 distinct RBPs** — co-switch modules are broadly enriched for shared RBP binding-site switching rather than being driven by one or two factors. The most reproducible regulators recurred across many regions: KHDRBS1 (SAM68) was a significant regulon in 8/10 regions, and A1CF, KHDRBS3, RBMS3, PPIE, RNASEL, and U2AF2 in 7/10 each, with CPEB2/CPEB4, KHDRBS2, NUDT21 (the 3′-end/APA factor), and ZCRB1 in 6/10. The strongest single-region regulons were neuronal splicing/UTR factors: in GTEx frontal cortex module M008, PPIE (enrichment 3.22×, q = 1.4×10⁻⁴⁰), the neuronal ELAV proteins ELAVL4/ELAVL3/ELAVL2 (3.0×, 1.7×, 3.2×), KHDRBS1 (2.0×), HNRNPD (2.1×), SYNCRIP (2.0×), and IGF2BP3 (2.2×) were all over-represented; module M007 was enriched for RBMS1/RBMS3, CPEB2, A1CF, and ADAR. The recurrence of 3′UTR/APA regulators (NUDT21, CPEB2/4, ELAV, RBMS) directly mirrors the 3′UTR remodeling seen in the coding-consequence analysis. This is an **exploratory, hypothesis-generating** result: the motif-presence call is a computational prediction of altered RBP binding, not measured binding, so significant RBPs are *candidate* regulons for experimental follow-up.

### Figure and Table Notes
- Potential supplementary table: `07_rbp_regulation/_m/rbp/rbp_regulon.parquet` (rendered `07_rbp_regulation/_m/rbp/RBP_REGULON.md`)
  - Rationale: full per-region module × RBP enrichment supporting the candidate-regulon claim.
  - Key columns: `region`, `module_id`, `rbp`, `module_size`, `n_switched`, `pool_switched`, `enrichment`, `q`, `go_invisible`.
- Potential supplementary table: `07_rbp_regulation/_m/rbp/rbp_switch_calls.parquet` (per gene × RBP gain/loss calls).
- Potential main/supplementary figure: a recurrence panel of the most reproducible RBP regulons (RBP × number of regions with a significant module), optionally paired with a frontal-cortex M008 highlight showing the neuronal ELAV/PPIE/KHDRBS1 cluster. Not yet generated.

### Reproducibility Information
- Analysis directory: `07_rbp_regulation/_m/rbp/`.
- Primary scripts: `isograph_benchmark/real_data/rbp_scan.py` (stage 1, motif env), `isograph_benchmark/real_data/rbp_regulon.py` (stage 2, isograph env).
- Execution command: `sbatch 07_rbp_regulation/_h/02.rbp_regulon.sh` (stage 1 in `motif` env, stage 2 in `isograph` env).
- Parameters: motif hit p-threshold 1×10⁻⁴, pseudocount 0.1, flat background, sense-strand only, U→T; hypergeometric over-representation, BH across module × RBP tests, minimum module size 3, switch-gene FDR default 0.05.
- Output files: `rbp_counts.parquet` (stage 1, 13:09), `rbp_switch_calls.parquet` + `rbp_regulon.parquet` + `RBP_REGULON.md` (stage 2, 13:27).
- Execution date: 2026-07-18.
- Git commit: baseline `86a3ef4`; analysis code committed in the same session (uncommitted at run time).
- Compute environment: PSC Bridges-2 RM-shared, account bio260021p.
  - Stage 1 env `/ocean/projects/bio260021p/shared/opt/envs/motif`: Python 3.11.15, MOODS-python 1.9.4.1, pandas 3.0.3, pyfaidx.
  - Stage 2 env `/ocean/projects/bio260021p/shared/opt/envs/isograph`: Python 3.12.13, numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6.
- Random seed: not applicable (deterministic scan + hypergeometric; no sampling).
- Missing reproducibility information: per-run package versions are not echoed into SLURM logs; recorded here from the runtime conda envs.

### Limitations and Integration Notes
Scope is deliberately light: only the mature transcript sequence is scanned (3′UTR/exonic single-stranded sites), not intronic splice-site flanks, so splicing-regulatory RBPs acting on introns are undercounted — the full genomic-intronic scope (hg38 FASTA + splice-site flanks) is the outstanding heavier extension. Motif presence is a sequence prediction, not measured binding or CLIP evidence, and ATtRACT PWMs are short and degenerate, so individual calls are noisy; the reproducible cross-region recurrence (e.g. KHDRBS1 in 8/10 regions) is the more trustworthy signal than any single module–RBP pair. This analysis is the trans-regulatory complement to the two cis analyses: coding-consequence (`06_switch_mechanism/_m/SWITCH_CONSEQUENCE_SUMMARY.md`) established the 3′UTR-remodeling substrate, and sQTL/eQTL colocalization (`05_genetic_anchoring/_m/coloc/`) established cis genetic anchoring — the shared-RBP result offers a candidate trans mechanism for *why* the module's genes co-switch. RBP motifs: ATtRACT database [@doi:10.1093/database/baw035]; scanning with MOODS [@doi:10.1093/bioinformatics/btp554]. [citation needed: GENCODE v47].
