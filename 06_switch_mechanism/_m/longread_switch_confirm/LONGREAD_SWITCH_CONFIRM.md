# Long-read confirmation of the switch layer — interpretation

Written 2026-09-12 from `summary.json` here (`longread_switch_confirm.py`, wrapper
`06_switch_mechanism/_h/07`) and the matched nulls in
`../switch_orthogonal_confirm/{GLOBAL_NULL.md,anchored_summary.json,signal_coloc/}`. Until
now these results were carried only in `TODO.md`.

## Question

Do IsoGraph's GTEx cortical aging switch genes and transcript pairs exist, and behave as
switches, on an independent platform?

## Data

ONT long-read DLPFC BA9/46, Aguzzoli-Heberle et al. 2024 (Zenodo 8180677, Bambu
quantification), **n = 12** (6 AD, 6 control). Independent platform, laboratory, cohort and
quantifier. Compared against GTEx frontal cortex BA9 and cortex switch sets (coding pairs only;
count ≥ 5 in ≥ 3 samples; isoform fraction ≥ 0.05).

## Result

| Quantity | Frontal cortex BA9 | Cortex | Pooled |
|---|---:|---:|---:|
| Switch genes expressed in long read | 6,869 / 7,004 | 1,654 / 1,676 | 8,523 |
| Expressed genes confirmed multi-isoform with both switch isoforms | 4,103 (59.7%) | 1,051 (63.5%) | **60.5%** |
| Switch transcript pairs detected | 16,394 / 61,302 (26.7%) | 4,512 / 15,878 (28.4%) | **27.1%** |

**Is detected-pair behaviour above a null?** Switch-like = anti-correlated usage across the 12
samples. Rates below use the **detected-pair** denominator and come with matched nulls:

| Pair set | switch-like | matched null | p |
|---|---:|---:|---:|
| All IsoGraph switch pairs vs abundance-matched non-switch pairs, same genes | 0.647 | 0.639 | 0.022 |
| Genetically anchored pairs vs matched switch pairs | **0.453** | 0.252 | 5e-4 |
| Anchored pairs, anchored isoform usably expressed | **0.600** | 0.304 | 5e-4 |
| Signal-level coloc nominations, anchored pairs (n = 124 detected) | **0.565** | 0.237 | 5e-4 |
| Signal-level, anchored isoform usably expressed (n = 61) | **0.672** | 0.299 | 5e-4 |

## Interpretation

- **Existence is well supported; switch behaviour of the layer as a whole is not.** About 60%
  of expressed switch genes show both isoforms in long read, but the all-pairs switch-like
  rate is only 0.0075 above its null. Within-gene isoform fractions sum to one, so any two
  isoforms of a gene are anti-correlated by construction — the null sits at 0.639 before any
  biology — which makes a bare negative correlation close to vacuous.
- **The informative comparisons are within the switch universe.** Genetically anchored pairs
  beat matched switch pairs by ~2× (0.453 vs 0.252; 0.600 vs 0.304 at usable abundance), and
  the signal-level nominations reproduce that (0.565 vs 0.237). That is the claim this arm can
  carry: the anchored subset behaves as switches in independent long read, as a set.
- **The retired 8.4%.** `pair_switch_like_rate` = 0.0844 in `summary.json` is switch-like pairs
  over *all* prespecified pairs, including the ~73% never detected at n = 12. It is a detection
  rate multiplied by a concordance rate, compares to nothing, and must not be quoted.
- **Single loci are not confirmed here.** SNCA's anchored isoform is 0.29% of the gene's
  long-read output and CTSH's 0.4%; both fail the usable-abundance qualification. For SNCA this
  is an assay limitation — the BrainSEQ short-read junction test measures the exact contrast at
  0.189 (DLPFC) and 0.234 (caudate) minor-form usage — so the 0.29% must never be cited as
  evidence against SNCA (`../junction_coloc_confirm/`).

## Limits

n = 12 bounds power for any per-pair statement; tissue is DLPFC against GTEx cortical regions,
not a matched region; half the donors have AD, so age-related switching is read against a
disease-mixed reference; Bambu and RSEM differ, which is the point of the arm but also means
detection depends on annotation agreement.

## Confidence

High for gene-level existence; moderate for the anchored-set switch behaviour (small n, but a
strict matched null and two independent nomination layers agree); not informative for
individual loci.
