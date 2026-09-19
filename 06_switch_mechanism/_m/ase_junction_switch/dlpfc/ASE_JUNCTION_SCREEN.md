# Allele-aware junction screen — dlpfc

Can reads that cross an isoform-specific splice junction **and** carry a phased heterozygous site support a within-donor allelic test of an isoform switch? This is the two-sided rescue of the exonic phASER screen (`ase_switch_direction.py`), counted from the WASP-tagged STAR BAMs on Quest.

## Pre-registered gate (unchanged from the exonic screen)

- per pair: at least **30** donors, each with at least **5** informative fragments, informative on **both** isoforms
- to proceed: at least **30** pairs pass in the gate family (pairs in genes with an all_samples switch-QTL, q < 0.05)

## Result

| family | pairs | genes | two-sided | any informative | **pass** | genes passing |
|---|---|---|---|---|---|---|
| gate_family | 2,219 | 163 | 1,682 | 2,170 | **1,014** | 126 |
| coloc_nominated | 251 | 18 | 199 | 229 | **102** | 15 |
| all_pairs | 16,206 | 1,324 | 12,566 | 15,343 | **5,849** | 855 |

**The gate passes: 1014 gate-family pairs clear it.** Stage 4 (allele orientation at the lead, then the beta-binomial contrast) may be built; its design is in the module docstring.

## Quality

- samples counted: 498 (missing: 0)
- haplotype-conflict rate among fragments: 0.0001
- WASP-fail rate among fragments: 0.0030; fragments with a het site but no WASP tag: 0
- het sites phased from GT because PW was absent: 0

## Caveats

- Cell composition cancels within a donor only if allelic effects do not differ by cell type.
- Some reference mapping bias survives WASP.
- The per-sample WASP VCFs carry `0|0`/`1|1` rows as well as hets (the `-g het` filter was applied before sample subsetting); heterozygous sites here are taken from phASER's `PW` genotype, so this does not reach the counts.

## Top gate-family pairs

| gene | T1 | T2 | donors (>= min) | both-isoform donors | het at lead | isoforms | pass |
|---|---|---|---|---|---|---|---|
| ENSG00000133816.19 | ENST00000683283.1 | ENST00000525618.5 | 425 | 53 | 109 | 2 | **yes** |
| ENSG00000133816.19 | ENST00000256194.8 | ENST00000525618.5 | 424 | 53 | 108 | 2 | **yes** |
| ENSG00000059122.17 | ENST00000253928.14 | ENST00000570752.5 | 422 | 15 | 206 | 2 | **yes** |
| ENSG00000059122.17 | ENST00000253928.14 | ENST00000574985.1 | 422 | 0 | 206 | 1 | no |
| ENSG00000059122.17 | ENST00000416288.6 | ENST00000570752.5 | 422 | 15 | 206 | 2 | **yes** |
| ENSG00000059122.17 | ENST00000416288.6 | ENST00000574985.1 | 422 | 0 | 206 | 1 | no |
| ENSG00000133816.19 | ENST00000683283.1 | ENST00000530691.5 | 421 | 2 | 108 | 2 | **yes** |
| ENSG00000133816.19 | ENST00000256194.8 | ENST00000530691.5 | 420 | 2 | 107 | 2 | **yes** |
| ENSG00000059122.17 | ENST00000344592.9 | ENST00000253928.14 | 413 | 14 | 201 | 2 | **yes** |
| ENSG00000059122.17 | ENST00000344592.9 | ENST00000416288.6 | 413 | 14 | 201 | 2 | **yes** |
| ENSG00000115306.17 | ENST00000333896.5 | ENST00000467371.1 | 412 | 217 | 181 | 2 | **yes** |
| ENSG00000133816.19 | ENST00000675839.1 | ENST00000525618.5 | 411 | 51 | 99 | 2 | **yes** |
| ENSG00000133816.19 | ENST00000675839.1 | ENST00000530691.5 | 409 | 2 | 97 | 2 | **yes** |
| ENSG00000115306.17 | ENST00000356805.9 | ENST00000602898.1 | 408 | 1 | 181 | 2 | **yes** |
| ENSG00000115306.17 | ENST00000356805.9 | ENST00000496323.1 | 408 | 0 | 181 | 1 | no |
| ENSG00000133816.19 | ENST00000707072.1 | ENST00000525618.5 | 408 | 147 | 109 | 2 | **yes** |
| ENSG00000115306.17 | ENST00000333896.5 | ENST00000602898.1 | 407 | 1 | 181 | 2 | **yes** |
| ENSG00000115306.17 | ENST00000356805.9 | ENST00000467371.1 | 401 | 0 | 181 | 1 | no |
| ENSG00000166501.14 | ENST00000321728.12 | ENST00000498739.1 | 397 | 8 | 237 | 2 | **yes** |
| ENSG00000182871.16 | ENST00000651438.1 | ENST00000459895.1 | 396 | 0 | 174 | 1 | no |
| ENSG00000182871.16 | ENST00000355480.10 | ENST00000459895.1 | 396 | 0 | 174 | 1 | no |
| ENSG00000182871.16 | ENST00000359759.8 | ENST00000459895.1 | 396 | 0 | 174 | 1 | no |
| ENSG00000059122.17 | ENST00000253928.14 | ENST00000573564.5 | 395 | 0 | 198 | 1 | no |
| ENSG00000059122.17 | ENST00000416288.6 | ENST00000573564.5 | 395 | 0 | 198 | 1 | no |
| ENSG00000166501.14 | ENST00000643927.1 | ENST00000498739.1 | 395 | 8 | 238 | 2 | **yes** |

