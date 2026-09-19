# 08 — Integration

**Question:** what do the layers say together about specific genes, programs and regulators?

Every analysis here reads more than one of stages 04–07 at once, which is the only reason it has
its own stage: placed anywhere earlier, it would read a stage that had not run yet. Until
2026-09-15 these lived in the stage whose question they finish (the deep dive and the SCZ age
projection in 05, functional preservation in trust, the RBP panel in 07) and had to be run out of
order by hand.

## Order

Run the whole stage with `bash 08_integration/_h/run_stage.sh` (add `--dry-run` to print the plan).
The leading number of a wrapper is its tier; steps in one tier run in parallel.

| Step | Wrapper | Reads | Produces |
|---|---|---|---|
| 01a | `build_deep_dive` | 05 coloc events + direction, 06 clinical consequence, 07 RBP calls and regulons | Per-gene vignettes, `deep_dive_panel`, `deep_dive_{rbp,exon_clinical,literature}`, and the genetic-anchoring figure set. The per-event table is written earlier, by `05_genetic_anchoring/_h/05c.deep_dive_events.sh`, because stage 06 reads it |
| 01b | `scz_age_projection` | 02 fits, 05 coloc events + MAGMA, 07 combined regulons | Are age-sensitive switch programs disrupted in SCZ? (the *convergence* sub-test is RETRACTED 2026-08-30 — wrong background; see `05_genetic_anchoring` `module_coloc_convergence`) |
| 01c | `replication_functional` | 04 matched pairs, 03 enrichment + composition, 06 switch consequence | Functional preservation of cross-cohort matched modules against size-matched random pairs (array: {isograph, wgcna} × {linear, spline}) |
| 02a | `build_rbp_target_panel` | 01a panel, 07 regulon + binding | Wet-lab perturbation panel per RBP arm (KHDRBS1, NONO, ELAVL1) |
| 03a | `rbp_pair_assayability` | 02a | Bench assayability of each panel's switch pairs |

`_m/` holds `deep_dive/`, `scz_age_projection/`, `functional_preservation/` and
`rbp_target_panel/`. The SCZ projection's genotype extract (`scz_age_projection/genotypes/`) is
controlled-access and gitignored.

**CLIs:** `isograph_benchmark/real_data/{gene_deep_dive,scz_age_projection,replication_functional,rbp_target_panel,rbp_pair_assayability}.py`.

## Display items

**Table 2** and S8–S12 (deep dive), **Fig 4E** and Table S18 (SCZ age projection), and the
Fig 4 panels built by `01a`.
