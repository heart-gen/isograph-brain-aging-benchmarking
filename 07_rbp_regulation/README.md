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

829 significant module × RBP tests (q<0.05) across 129 RBPs; recurrent neuronal
3′UTR/splicing factors (KHDRBS1 8/10 regions; A1CF/KHDRBS3/RBMS3/PPIE/RNASEL/U2AF2 7/10;
ELAV/CPEB families), ~245 of them in GO-invisible modules. The GLM regulon test applies
an estimability gate that excludes separated cells before BH correction.

`_m/neuronal_clip/` is gitignored (large downloaded CLIP tracks); the manifests are
tracked under `_m/neuronal_clip_manifests/`.

**CLIs:** `isograph_benchmark/real_data/{rbp_motif_families,rbp_scan,rbp_scan_intronic,rbp_regulon,rbp_binding,rbp_pair_assayability,rbp_target_panel,neuronal_clip_*,nova*}.py`.

## Display items

S-real-6 `figRbpRegulon`; supplementary table S10 (per-gene RBP regulators).
Analysis specs: `docs/RBP_NEURONAL_CLIP_ANALYSIS_SPEC.md`,
`docs/NOVA2_PERTURBATION_ANALYSIS_SPEC.md`.
