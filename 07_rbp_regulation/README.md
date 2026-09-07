# 07 — RBP regulation

**Question:** what *trans* factors could drive the co-switching — do modules share
candidate RNA-binding-protein regulators?

This is the mechanism arm. It is **motif prediction plus public CLIP**, not new binding
data, and the language in the paper should say so: these are candidate regulators.

## Order

| Step | Wrapper | Produces |
|---|---|---|
| 01 | `rbp_motif_families` | ATtRACT PWM family collapse |
| 02–03 | `rbp_regulon`, `rbp_regulon_intronic` | Per-module × RBP over-representation in switched exons (mature / intronic / combined), GC-binned background |
| 04–05 | `rbp_binding_fetch`, `rbp_binding` | Public CLIP binding evidence |
| 06–07 | `rbp_pair_assayability`, `build_rbp_target_panel` | Assayability of switch pairs; wet-lab target panel |
| 08–11 | `neuronal_clip_{freeze,windows,motif_qc,overlap}` | Neuronal CLIP freeze, windows, motif QC, overlap |
| 12–14 | `nova_family_renomination`, `nova2_ctag_clip`, `nova2_perturbation` | NOVA family arm |

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
