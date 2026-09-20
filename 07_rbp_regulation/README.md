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

**384** significant module × RBP tests (hypergeometric q<0.05) across **125** RBPs of
**160** tested, spanning 12 of 16 regions; **61** of the 384 fall in GO-invisible modules.
Under the opportunity-adjusted GLM (transcript length, GC, UTR composition) the count is
**416** hits across **124** RBPs in 13 regions (159 GO-invisible). **The two arms are not
nested: only 89 cells are significant in both**, and 327 of the 416 adjusted hits are not
hypergeometric hits at all. Quote the GLM number for a regulon claim and the
hypergeometric one for screening breadth, but never present the adjusted count as "the
survivors" of the raw one — the adjustment reshuffles the ranking rather than thinning it.

Recurrent factors, by the number of the 16 regions in which the RBP has a hypergeometric
hit: **PPRC1 at 8/16** and **IGF2BP3 at 7/16**; DHX9, G3BP2, KHDRBS1, RBM41, PABPC5 and
NOVA2 at 6/16. On the covariate-adjusted GLM the most recurrent are RBFOX2 at 6/16, then
IGF2BP3, PABPC5, TARDBP and SART3 at 5/16. The GLM regulon test applies an estimability
gate that excludes separated cells (510 zero-cell of 11,040) before BH correction.

*(Counts recomputed 2026-09-19 from `_m/rbp/rbp_regulon.parquet` on the switching-filter
re-run at Leiden resolution 2.0. The previous section quoted the legacy res-5.0 production
— 714 hypergeometric / 310 adjusted over 34,240 cells in 14 regions — which does not apply
to the current partitions. Resolution 2.0 gives fewer, larger modules, so the cell count
fell roughly threefold. The ELAV/CPEB recurrence story did not survive: the legacy list
led with KHDRBS1, ELAVL4 and PPRC1 at 8/14.)*

`_m/neuronal_clip/` is gitignored (large downloaded CLIP tracks); the manifests are
tracked under `_m/neuronal_clip_manifests/`.

**CLIs:** `isograph_benchmark/real_data/{rbp_motif_families,rbp_scan,rbp_scan_intronic,rbp_regulon,rbp_binding,rbp_pair_assayability,rbp_target_panel,neuronal_clip_*,nova*}.py`.

## Display items

S-real-6 `figRbpRegulon`; supplementary table S10 (per-gene RBP regulators).
Analysis specs: `docs/RBP_NEURONAL_CLIP_ANALYSIS_SPEC.md`,
`docs/NOVA2_PERTURBATION_ANALYSIS_SPEC.md`.
