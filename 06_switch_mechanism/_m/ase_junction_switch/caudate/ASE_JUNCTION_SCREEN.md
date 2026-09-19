# Allele-aware junction screen — caudate

Can reads that cross an isoform-specific splice junction **and** carry a phased heterozygous site support a within-donor allelic test of an isoform switch? This is the two-sided rescue of the exonic phASER screen (`ase_switch_direction.py`), counted from the WASP-tagged STAR BAMs on Quest.

## Pre-registered gate (unchanged from the exonic screen)

- per pair: at least **30** donors, each with at least **5** informative fragments, informative on **both** isoforms
- to proceed: at least **30** pairs pass in the gate family (pairs in genes with an all_samples switch-QTL, q < 0.05)

## Result

| family | pairs | genes | two-sided | any informative | **pass** | genes passing |
|---|---|---|---|---|---|---|
| gate_family | 8,480 | 661 | 6,449 | 8,178 | **3,896** | 531 |
| coloc_nominated | 624 | 45 | 489 | 606 | **258** | 37 |
| all_pairs | 48,412 | 3,932 | 37,500 | 45,566 | **18,266** | 2,618 |

**The gate passes: 3896 gate-family pairs clear it.** Stage 4 (allele orientation at the lead, then the beta-binomial contrast) may be built; its design is in the module docstring.

## Quality

- samples counted: 486 (missing: 0)
- haplotype-conflict rate among fragments: 0.0001
- WASP-fail rate among fragments: 0.0029; fragments with a het site but no WASP tag: 0
- het sites phased from GT because PW was absent: 0

## Caveats

- Cell composition cancels within a donor only if allelic effects do not differ by cell type.
- Some reference mapping bias survives WASP.
- The per-sample WASP VCFs carry `0|0`/`1|1` rows as well as hets (the `-g het` filter was applied before sample subsetting); heterozygous sites here are taken from phASER's `PW` genotype, so this does not reach the counts.

## Top gate-family pairs

| gene | T1 | T2 | donors (>= min) | both-isoform donors | het at lead | isoforms | pass |
|---|---|---|---|---|---|---|---|
| ENSG00000136205.18 | ENST00000311160.14 | ENST00000428457.1 | 472 | 0 | 22 | 1 | no |
| ENSG00000136205.18 | ENST00000705350.1 | ENST00000428457.1 | 472 | 0 | 22 | 1 | no |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000491228.1 | 468 | 0 | 157 | 1 | no |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000493501.1 | 467 | 0 | 157 | 1 | no |
| ENSG00000136205.18 | ENST00000457718.6 | ENST00000428457.1 | 464 | 0 | 22 | 1 | no |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000649153.1 | 459 | 5 | 156 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000650252.1 | 459 | 6 | 156 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000649153.1 | 459 | 339 | 157 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000650252.1 | 453 | 325 | 157 | 2 | **yes** |
| ENSG00000197077.14 | ENST00000406486.8 | ENST00000461374.1 | 452 | 0 | 141 | 1 | no |
| ENSG00000197077.14 | ENST00000358431.8 | ENST00000461374.1 | 447 | 0 | 136 | 1 | no |
| ENSG00000197077.14 | ENST00000406486.8 | ENST00000494730.2 | 447 | 0 | 138 | 1 | no |
| ENSG00000154358.23 | ENST00000680850.1 | ENST00000664957.1 | 440 | 294 | 221 | 2 | **yes** |
| ENSG00000154358.23 | ENST00000570156.7 | ENST00000664957.1 | 440 | 294 | 221 | 2 | **yes** |
| ENSG00000154358.23 | ENST00000680850.1 | ENST00000474237.2 | 439 | 35 | 221 | 2 | **yes** |
| ENSG00000154358.23 | ENST00000570156.7 | ENST00000474237.2 | 439 | 35 | 221 | 2 | **yes** |
| ENSG00000197077.14 | ENST00000358431.8 | ENST00000494730.2 | 438 | 0 | 129 | 1 | no |
| ENSG00000154358.23 | ENST00000662438.1 | ENST00000664957.1 | 436 | 294 | 221 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000443831.1 | 436 | 73 | 152 | 2 | **yes** |
| ENSG00000154310.17 | ENST00000436636.7 | ENST00000496492.5 | 436 | 0 | 91 | 1 | no |
| ENSG00000154310.17 | ENST00000341852.10 | ENST00000496492.5 | 436 | 0 | 91 | 1 | no |
| ENSG00000073910.23 | ENST00000642040.1 | ENST00000641614.1 | 436 | 14 | 141 | 2 | **yes** |
| ENSG00000073910.23 | ENST00000642040.1 | ENST00000490410.1 | 436 | 0 | 141 | 1 | no |
| ENSG00000073910.23 | ENST00000645780.1 | ENST00000641614.1 | 436 | 0 | 141 | 1 | no |
| ENSG00000073910.23 | ENST00000645780.1 | ENST00000490410.1 | 436 | 0 | 141 | 1 | no |

