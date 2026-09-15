# 07 — RBP regulation

**Question:** what *trans* factors could drive the co-switching — do modules share
candidate RNA-binding-protein regulators?

This is the mechanism arm. It is **motif prediction plus public CLIP**, not new binding
data, and the language in the paper should say so: these are candidate regulators.

## Order

Run the whole stage with `bash 07_rbp_regulation/_h/run_stage.sh` (add `--dry-run` to print the
plan). The leading number of a wrapper is its tier; steps in one tier run in parallel.

| Step | Wrapper | Waits on | Produces |
|---|---|---|---|
| 01a | `rbp_motif_families` | — | ATtRACT PWM family collapse |
| 02a | `rbp_regulon` | 01a, stage 03 | Motif scan + per-module × RBP over-representation in switched exons (mature scope), GC-binned background |
| 03a | `rbp_regulon_intronic` | 02a | Intronic and combined scopes |
| 03b | `rbp_binding_fetch` | 02a (**login node**) | ENCODE eCLIP peaks for every nominated RBP; the runner holds 04a until it has run |
| 04a | `rbp_binding` | 01a, 02a, 03b | Public CLIP binding evidence joined to the regulon calls |
| 04b | `neuronal_clip_freeze` | 03a; `inputs/_h/download_neuronal_clip.sh` | Frozen module-RBP candidate set + dataset/context QC manifests |
| 05a | `neuronal_clip_windows` | 04b | Pair-specific splice-flank windows + matched within-gene controls |
| 06a | `neuronal_clip_motif_qc` | 05a | Intronic-opportunity audit; NOVA-family coordinates |
| 07a | `neuronal_clip_overlap` | 06a | Assay-callable window calls in human TDP-43 / PTBP2 contexts |
| 07b | `nova_family_renomination` | 06a | NOVA-family nominations with opportunity-adjusted module enrichment |
| 08a | `nova2_ctag_clip` | 07b | NOVA2 cTag-CLIP overlap |
| 09a | `nova2_perturbation` | 08a | Nova2-cKO splice events localized within the frozen windows |

The freeze (04b) guards the candidate count and identity hash recorded in
`configs/neuronal_clip.yaml`, as do the NOVA steps. A re-fit that changes the nominations fails
those guards by design; update the expected values deliberately, never to get past them.

The wet-lab target panel and the pair assayability read the per-gene deep dive, so they moved to
`08_integration/_h/02a` and `03a` on 2026-09-15.

## Key results

**714** significant module × RBP tests (hypergeometric q<0.05) across **138** RBPs of
**160** tested, spanning 14 regions; **149** of the 714 fall in GO-invisible modules.
Under the opportunity-adjusted GLM (transcript length, GC, UTR composition) the count is
**310** hits across **107** RBPs in 12 regions (50 GO-invisible). **The two arms are not
nested: only 43 cells are significant in both**, and 267 of the 310 adjusted hits are not
hypergeometric hits at all. Quote the GLM number for a regulon claim and the
hypergeometric one for screening breadth, but never present 310 as "the survivors" of
714 — the adjustment reshuffles the ranking rather than thinning it.

Recurrent neuronal 3′UTR/splicing factors, by the number of the 14 regions in which the
RBP has a hypergeometric hit: **KHDRBS1, ELAVL4 and PPRC1 at 8/14**; CPEB4, RBMS3,
RNASEL, ELAVL3, RBM14 and ADAR at 7/14 — i.e. the ELAV and CPEB families recur, as does
KHDRBS1. On the covariate-adjusted GLM the most recurrent are DHX58, RBFOX2 and SART3
at 6/14. The GLM regulon test applies an estimability gate that excludes separated cells
(2,529 of 34,240) before BH correction.

*(Counts recomputed 2026-09-03 from `_m/rbp/rbp_regulon.parquet`. A previous version of
this section read "829 ... across 129 RBPs ... ~245 GO-invisible" with a 7/10 recurrence
list; none of those reproduced, and the denominator is 14 regions, not 10.)*

`_m/neuronal_clip/` is gitignored (large downloaded CLIP tracks); the manifests are
tracked under `_m/neuronal_clip_manifests/`.

**CLIs:** `isograph_benchmark/real_data/{rbp_motif_families,rbp_scan,rbp_scan_intronic,rbp_regulon,rbp_binding,rbp_pair_assayability,rbp_target_panel,neuronal_clip_*,nova*}.py`.

## Display items

S-real-6 `figRbpRegulon`; supplementary table S10 (per-gene RBP regulators).
Analysis specs: `docs/RBP_NEURONAL_CLIP_ANALYSIS_SPEC.md`,
`docs/NOVA2_PERTURBATION_ANALYSIS_SPEC.md`.
