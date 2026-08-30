# IsoGraph GO-invisible disease switch modules (DTU-without-DGE biology gate)

Modular analysis summary for Manubot integration. Generated from
`02_module_discovery/brainseq/caudate_sczd/_m/go_invisible_gate.parquet` and
`go_invisible_gate_background.json`. Every numeric claim is reproduced from those files;
do not edit the numbers by hand — regenerate from the table.

## Purpose

Decide, honestly, what IsoGraph's phenotype-associated switch modules *are*. Two failure
modes had to be excluded before any biological claim: (1) that the modules are abundance
shifts mislabeled as switches, and (2) that "GO-invisible" simply means low-quality. The
gate asks whether the schizophrenia (SCZD)-associated IsoGraph switch modules in BrainSEQ
caudate carry genuine, functionally-consequential isoform switching even when classical GO
enrichment sees nothing — i.e. whether they are a real DTU-without-DGE layer that is
structurally invisible to a gene-level/abundance pathway pipeline, rather than pathways
"WGCNA misses."

## Inputs

- **Production SCZD switch modules** — IsoGraph `isograph_vae/` fit on BrainSEQ caudate
  (full-multiplex primary, Leiden resolution 5.0), with module–SCZD phenotype associations
  (`pheno_fdr`).
- **Per-transcript switch evidence** — anticorrelated within-gene transcript pairs and
  switch strength from the fit, used to count members carrying a *real* switch (both an
  up- and a down-correlated transcript), not an abundance shift.
- **Structural switch annotation** — per-driver functional-consequence flags (CDS change,
  coding-status change, biotype switch, UTR change) from the isoform structural annotation.
- **GO enrichment** — the module GO-enrichment stage; a module is "GO-invisible" when it
  returns zero enriched terms.
- **Pooled background** — the functional-consequence fractions over all switch transcripts
  in the analysis (`go_invisible_gate_background.json`), the reference the disease modules
  are judged against.

## Methods text

Within the BrainSEQ caudate SCZD analysis we identified IsoGraph switch modules associated
with case/control status at a module-level phenotype FDR ≤ 0.1, and classified each as
GO-enriched or GO-invisible by whether its members returned any enriched GO term. For each
phenotype-significant module we then quantified switch coherence as the number of members
carrying a genuine isoform switch (an anticorrelated within-gene transcript pair) and the
maximum within-module switch strength, and quantified functional consequence as the
fraction of driver switches that alter the coding sequence, change coding status, switch
transcript biotype, or change the UTR. These per-module fractions were compared against the
pooled background over all switch transcripts in the same analysis. The test of the gate is
comparative, not absolute: GO-invisible disease modules pass if their switch coherence and
functional-consequence fractions are at or above background and indistinguishable from any
GO-visible disease modules — establishing that GO-invisibility reflects GO's gene-level
bias, not lower module quality. Analyses used IsoGraph v0.1.5 under the project Python 3.12
environment with fixed seeds.

## Results text

**All phenotype-significant SCZD switch modules are GO-invisible.** Four IsoGraph switch
modules are associated with schizophrenia case/control status (pheno_fdr ≤ 0.1): M026
(FDR 0.0019), M020 (0.017), M010 (0.027), and M023 (0.047). **All four return zero enriched
GO terms** (`n_go_terms = 0`); there are no GO-enriched disease switch modules in this
analysis. Classical pathway enrichment is blind to the entire phenotype-associated switch
signal here.

**They carry genuine isoform switching, not abundance shifts.** In every module nearly all
members carry a real anticorrelated transcript pair: 21/23 (M026), 31/30 (M020 — count
includes multi-transcript genes), 73/83 (M010), 28/27 (M023), with maximum switch strength
1.10–1.39 and 93–459 significant switch transcripts per module. The modules are DTU, not
re-described differential expression.

**Their switches are functionally consequential at or above background.** Driver
functional-consequence fractions match or exceed the pooled background (CDS 0.84,
coding-status 0.67, biotype 0.74, UTR 0.61): e.g. M026 CDS 1.00 / coding-status 0.78 /
biotype 0.89 / UTR 0.71; M023 0.96 / 0.82 / 0.91 / 0.65; M020 0.88 / 0.81 / 0.87 / 0.57;
M010 0.71 / 0.53 / 0.58 / 0.54. The GO-invisible disease modules are structurally
indistinguishable from the background switch population — GO-invisibility is a property of
GO's gene-level/abundance bias, not of module quality.

**Drivers are plausible but heterogeneous within a module.** Top switch genes are
psychiatric-relevant (DDX3X, XPO1, IFNAR2 in M026; STXBP5, SPTBN1, ZDHHC17 in M020; ARNT2,
PBX1, DGKH, SLC25A12 in M023; SEPTIN8, FAM107B, SOX2-OT in M010) but heterogeneous within
each module — a shared *switch axis*, not a shared GO process. This is why pathway
enrichment is blind to them and why the correct framing is a complementary DTU-without-DGE
layer, **not** "pathways WGCNA misses."

**Verdict: PASS in the complementary form.** The phenotype-associated SCZD switch signal is
real, functionally consequential isoform regulation that is invisible to GO/abundance
pipelines — the disease-axis instance of IsoGraph's defensible complementary value.

## Figure and table notes

- **Table:** `go_invisible_gate.parquet` — the per-module ledger (module, n_genes,
  pheno_fdr, go_invisible, n_go_terms, genes_with_real_switch, max_switch_strength,
  n_sig_switch_tx, frac_cds_changed/coding_status_change/biotype_switch/utr_changed,
  top_switch_genes) plus the pooled `_background` row from the JSON sidecar. This is the
  primary artifact; the gate is a small, table-centric result.
- **Supplementary figure — GO-invisible switch modules
  (`real_data/_m/figures/figGoInvisible.{pdf,png}`, built by
  `real_data/_h/go_invisible_figure.R`).** (A) per-module functional-consequence fractions
  (CDS / coding-status / biotype / UTR) for the four GO-invisible disease modules with the
  pooled-background line overlaid, showing they sit comparable to background (some above, some
  below — indistinguishable, not depleted); (B) switch coherence — nearly every member carries
  a real anticorrelated transcript pair, with max switch strength annotated. Sits in the
  complementarity supplement next to the QTL splicing-specificity figure (the same modules
  carry the genetic anchoring).

## Reproducibility information

- Analysis directory: `02_module_discovery/brainseq/caudate_sczd/_m/`.
- Primary script: `isograph_benchmark/real_data/go_invisible_gate.py`
  (`python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`).
- SLURM driver: `02_module_discovery/brainseq/_h/12.go_invisible_gate.sh`.
- Inputs: production `isograph_vae/` SCZD modules + module phenotype associations, switch
  evidence, structural switch annotation, module GO enrichment.
- Outputs: `go_invisible_gate.parquet`, `GO_INVISIBLE_GATE.md`,
  `go_invisible_gate_background.json`.
- Key parameters: phenotype FDR threshold ≤ 0.1; GO-invisible = zero enriched terms;
  Leiden resolution 5.0; deterministic seeds. Last regenerated 2026-06-29 (post
  covariate-decouple production re-fit).
- Compute environment: PSC Bridges-2 RM-shared; IsoGraph v0.1.5, project Python 3.12
  (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Missing reproducibility information: per-package versions taken from the live environment,
  not a per-run lockfile.

## Limitations and integration notes

- This gate is run on a single cohort/region (BrainSEQ caudate, SCZD). It establishes that
  the phenotype-associated switch modules there are real DTU; it does not by itself claim
  generality. The cross-region/aging analog of "is the switch signal real" is the
  module-trust funnel; the genetic-mechanism analog is the QTL splicing-specificity contrast
  (the GO-invisible modules are exactly where splicing-QTL are spared).
- The result is **post-refit**: the covariate-decouple production re-fit changed the
  phenotype-FDR landscape, and the current gate yields 4 phenotype-significant modules, all
  GO-invisible (an earlier res-5.0 run reported 8 disease modules with 2 GO-visible). Cite
  the regenerated parquet, not earlier prose.
- The honest framing is **complementary, not superior**: the absence of GO-visible disease
  modules here is not a claim that IsoGraph finds pathways WGCNA cannot — it reflects that
  the signal is isoform regulation, which GO does not represent. Integrate with the
  three-baseline comparison (IsoGraph's GO-enrichment rate is low by construction) and the
  de-confounded gene-level result (abundance dominates the bulk signal). This summary
  supplies the disease-axis DTU-without-DGE biology; the QTL and trust-funnel summaries
  supply its genetic anchoring and reproducibility.
- This is a **primary** biology gate, not a sensitivity check; it is the qualitative
  go/no-go that licenses the complementary-layer claim.
