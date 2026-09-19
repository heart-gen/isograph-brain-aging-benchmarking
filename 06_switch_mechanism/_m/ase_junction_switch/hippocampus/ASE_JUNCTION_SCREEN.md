# Allele-aware junction screen — hippocampus

Can reads that cross an isoform-specific splice junction **and** carry a phased heterozygous site support a within-donor allelic test of an isoform switch? This is the two-sided rescue of the exonic phASER screen (`ase_switch_direction.py`), counted from the WASP-tagged STAR BAMs on Quest.

## Pre-registered gate (unchanged from the exonic screen)

- per pair: at least **30** donors, each with at least **5** informative fragments, informative on **both** isoforms
- to proceed: at least **30** pairs pass in the gate family (pairs in genes with an all_samples switch-QTL, q < 0.05)

## Result

| family | pairs | genes | two-sided | any informative | **pass** | genes passing |
|---|---|---|---|---|---|---|
| gate_family | 5,529 | 429 | 4,119 | 5,305 | **2,219** | 322 |
| coloc_nominated | 621 | 46 | 486 | 605 | **221** | 34 |
| all_pairs | 55,508 | 4,450 | 43,369 | 52,052 | **17,206** | 2,673 |

**The gate passes: 2219 gate-family pairs clear it.** Stage 4 (allele orientation at the lead, then the beta-binomial contrast) may be built; its design is in the module docstring.

## Quality

- samples counted: 451 (missing: 0)
- haplotype-conflict rate among fragments: 0.0001
- WASP-fail rate among fragments: 0.0032; fragments with a het site but no WASP tag: 0
- het sites phased from GT because PW was absent: 0

## Caveats

- Cell composition cancels within a donor only if allelic effects do not differ by cell type.
- Some reference mapping bias survives WASP.
- The per-sample WASP VCFs carry `0|0`/`1|1` rows as well as hets (the `-g het` filter was applied before sample subsetting); heterozygous sites here are taken from phASER's `PW` genotype, so this does not reach the counts.

## Top gate-family pairs

| gene | T1 | T2 | donors (>= min) | both-isoform donors | het at lead | isoforms | pass |
|---|---|---|---|---|---|---|---|
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000491228.1 | 422 | 0 | 129 | 1 | no |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000493501.1 | 420 | 0 | 128 | 1 | no |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000650252.1 | 415 | 3 | 128 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000649153.1 | 414 | 3 | 128 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000649153.1 | 411 | 303 | 127 | 2 | **yes** |
| ENSG00000187240.17 | ENST00000375735.7 | ENST00000533027.1 | 408 | 0 | 55 | 1 | no |
| ENSG00000187240.17 | ENST00000650373.2 | ENST00000533027.1 | 408 | 0 | 55 | 1 | no |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000650252.1 | 405 | 276 | 126 | 2 | **yes** |
| ENSG00000187240.17 | ENST00000375735.7 | ENST00000648198.1 | 404 | 0 | 55 | 1 | no |
| ENSG00000187240.17 | ENST00000650373.2 | ENST00000648198.1 | 404 | 0 | 55 | 1 | no |
| ENSG00000155093.20 | ENST00000389418.9 | ENST00000648371.1 | 403 | 0 | 159 | 1 | no |
| ENSG00000155093.20 | ENST00000389416.8 | ENST00000648371.1 | 403 | 0 | 159 | 1 | no |
| ENSG00000155093.20 | ENST00000409483.5 | ENST00000648371.1 | 403 | 0 | 159 | 1 | no |
| ENSG00000187240.17 | ENST00000528670.5 | ENST00000649323.1 | 400 | 267 | 55 | 2 | **yes** |
| ENSG00000197077.14 | ENST00000406486.8 | ENST00000461374.1 | 392 | 0 | 126 | 1 | no |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000493501.1 | 391 | 216 | 125 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000443831.1 | ENST00000491228.1 | 390 | 21 | 123 | 2 | **yes** |
| ENSG00000081479.15 | ENST00000649046.1 | ENST00000443831.1 | 388 | 58 | 117 | 2 | **yes** |
| ENSG00000155093.20 | ENST00000389413.7 | ENST00000648371.1 | 386 | 0 | 150 | 1 | no |
| ENSG00000006071.16 | ENST00000683693.1 | ENST00000529967.6 | 385 | 219 | 138 | 2 | **yes** |
| ENSG00000197077.14 | ENST00000358431.8 | ENST00000461374.1 | 383 | 0 | 117 | 1 | no |
| ENSG00000182871.16 | ENST00000651438.1 | ENST00000459895.1 | 378 | 0 | 155 | 1 | no |
| ENSG00000182871.16 | ENST00000355480.10 | ENST00000459895.1 | 378 | 0 | 155 | 1 | no |
| ENSG00000182871.16 | ENST00000359759.8 | ENST00000459895.1 | 378 | 0 | 155 | 1 | no |
| ENSG00000006071.16 | ENST00000647086.1 | ENST00000529967.6 | 375 | 228 | 141 | 2 | **yes** |

