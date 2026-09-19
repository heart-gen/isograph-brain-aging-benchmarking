# Biology gate — GO-invisible switch modules (brainseq-sczd)

Phenotype-significant IsoGraph switch modules (pheno_fdr <= 0.1): **6** — 1 GO-enriched, **5 GO-invisible** (M004, M007, M006, M008, M009).

Reproduce: `python -m isograph_benchmark.real_data.go_invisible_gate --analysis brainseq-sczd`

## Per-module switch coherence

| module | go_invisible | n_genes | pheno_fdr | genes_with_real_switch | max_switch_strength | n_sig_switch_tx | top_switch_genes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M005 | False | 382 | 0.0243 | 348 | 1.814 | 1798 | RPAIN, COX11, GET4, TAF9, MAEA, AACS, RBX1, HARS1 |
| M004 | True | 432 | 0.003 | 415 | 1.824 | 1652 | PIN1, CWC27, CADPS, FBXO9, ILKAP, PLEKHO1, FRA10AC1, COX7A2L |
| M007 | True | 266 | 0.0064 | 262 | 1.83 | 1417 | TBCD, PRKAG2, BBIP1, DDX56, CZIB, DRG2, NPAS2, TTC19 |
| M006 | True | 343 | 0.0179 | 343 | 1.795 | 1567 | MBTPS1, YEATS2, ERCC3, N4BP2L2, INTS11, FLCN, CBR4, UBE2Z |
| M008 | True | 230 | 0.0197 | 228 | 1.785 | 1276 | PTRH2, COPS5, OAZ2, SPIRE2, DDRGK1, MPV17, WDR11, DDX19A |
| M009 | True | 170 | 0.0243 | 169 | 1.839 | 1101 | FIG4, LGI1, MARK3, MACROH2A1, PRMT8, NUP93, RUVBL1, RAE1 |

Pooled background functional-consequence fractions: cds_changed 0.887, coding_status_change 0.775, biotype_switch 0.835, utr_changed 0.518.

## Reading

- `genes_with_real_switch` near `n_genes` = members carry genuine DTU (both up- and down-correlated transcripts), not abundance shifts.
- Functional-consequence fractions at/above background and indistinguishable from GO-visible disease modules => GO-invisible modules are not lower quality; GO-invisibility reflects GO's gene-level/abundance bias, blind to DTU.
- `top_switch_genes` heterogeneous within a module (shared switch axis, not a shared GO process) => frame as a complementary DTU-without-DGE layer, NOT pathway discovery WGCNA misses.
