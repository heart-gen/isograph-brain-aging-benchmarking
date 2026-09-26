# Donor and tissue-sample sharing across the 16 regions and 17 analyses

16 regions are analysed by 17 analyses. The 4,322 sample-analysis rows cover 4,084 distinct tissue samples from 844 donors: BrainSEQ caudate is one region used by two analyses, so its control samples are counted by both. Donor ids are cohort-local (`BrNum` for BrainSEQ, `SUBJID` for GTEx) and are never matched across cohorts; `sample_id` is the tissue-specific id (`RNum` for BrainSEQ) and is checked above to be region-specific.

## Per region

| region | cohort | n_donors | n_tissue_samples | n_analyses |
| --- | --- | --- | --- | --- |
| BrainSEQ DLPFC | BrainSEQ | 222 | 222 | 1 |
| BrainSEQ caudate | BrainSEQ | 390 | 390 | 2 |
| BrainSEQ hippocampus | BrainSEQ | 238 | 238 | 1 |
| GTEx ACC BA24 | GTEx | 233 | 233 | 1 |
| GTEx amygdala | GTEx | 181 | 181 | 1 |
| GTEx caudate | GTEx | 300 | 300 | 1 |
| GTEx cerebellar hem. | GTEx | 277 | 277 | 1 |
| GTEx cerebellum | GTEx | 266 | 266 | 1 |
| GTEx cortex | GTEx | 270 | 270 | 1 |
| GTEx frontal ctx BA9 | GTEx | 269 | 269 | 1 |
| GTEx hippocampus | GTEx | 255 | 255 | 1 |
| GTEx hypothalamus | GTEx | 257 | 257 | 1 |
| GTEx n. accumbens | GTEx | 285 | 285 | 1 |
| GTEx putamen | GTEx | 254 | 254 | 1 |
| GTEx spinal cord C1 | GTEx | 204 | 204 | 1 |
| GTEx substantia nigra | GTEx | 183 | 183 | 1 |

## Per analysis

| analysis | region | cohort | phenotype | n_samples | n_donors |
| --- | --- | --- | --- | --- | --- |
| BrainSEQ caudate (aging) | BrainSEQ caudate | BrainSEQ | aging | 238 | 238 |
| BrainSEQ caudate (SCZD) | BrainSEQ caudate | BrainSEQ | diagnosis | 390 | 390 |
| BrainSEQ hippocampus | BrainSEQ hippocampus | BrainSEQ | aging | 238 | 238 |
| BrainSEQ DLPFC | BrainSEQ DLPFC | BrainSEQ | aging | 222 | 222 |
| GTEx amygdala | GTEx amygdala | GTEx | aging | 181 | 181 |
| GTEx ACC BA24 | GTEx ACC BA24 | GTEx | aging | 233 | 233 |
| GTEx caudate | GTEx caudate | GTEx | aging | 300 | 300 |
| GTEx cerebellar hem. | GTEx cerebellar hem. | GTEx | aging | 277 | 277 |
| GTEx cerebellum | GTEx cerebellum | GTEx | aging | 266 | 266 |
| GTEx cortex | GTEx cortex | GTEx | aging | 270 | 270 |
| GTEx frontal ctx BA9 | GTEx frontal ctx BA9 | GTEx | aging | 269 | 269 |
| GTEx hippocampus | GTEx hippocampus | GTEx | aging | 255 | 255 |
| GTEx hypothalamus | GTEx hypothalamus | GTEx | aging | 257 | 257 |
| GTEx n. accumbens | GTEx n. accumbens | GTEx | aging | 285 | 285 |
| GTEx putamen | GTEx putamen | GTEx | aging | 254 | 254 |
| GTEx spinal cord C1 | GTEx spinal cord C1 | GTEx | aging | 204 | 204 |
| GTEx substantia nigra | GTEx substantia nigra | GTEx | aging | 183 | 183 |

## Donors per cohort

| cohort | n_donors | median_regions | max_regions | n_in_one_region_only |
| --- | --- | --- | --- | --- |
| BrainSEQ | 451 | 2 | 3 | 213 |
| GTEx | 393 | 9 | 13 | 20 |

## Most-overlapping region pairs

| region_a | region_b | n_shared_donors | jaccard |
| --- | --- | --- | --- |
| GTEx caudate | GTEx n. accumbens | 249 | 0.741 |
| GTEx cerebellar hem. | GTEx caudate | 238 | 0.702 |
| GTEx putamen | GTEx caudate | 233 | 0.726 |
| GTEx caudate | GTEx frontal ctx BA9 | 232 | 0.688 |
| GTEx n. accumbens | GTEx cerebellar hem. | 228 | 0.683 |
| GTEx frontal ctx BA9 | GTEx n. accumbens | 226 | 0.689 |
| GTEx caudate | GTEx hypothalamus | 225 | 0.678 |
| GTEx n. accumbens | GTEx putamen | 224 | 0.711 |
| GTEx cerebellar hem. | GTEx cerebellum | 222 | 0.692 |
| GTEx cortex | GTEx caudate | 222 | 0.638 |

## Identifier notes

- BrainSEQ caudate (aging) and BrainSEQ caudate (SCZD) share 238 tissue samples within BrainSEQ caudate

## Interpretation

- **These are 17 applications over 16 regions, not 17 independent cohorts.** The GTEx regions are largely re-measurements of one donor pool, the three BrainSEQ regions are largely the same brains, and the two BrainSEQ caudate analyses are the same region: the aging analysis is the 238 controls, and the schizophrenia analysis is those same 238 control samples plus 152 patients.
- Any statement that a result 'replicates across regions' within a cohort is therefore a statement about the same donors measured in different tissue, not about independent samples. The cross-cohort comparisons (BrainSEQ vs GTEx) are the ones that carry independent donors, and they are also the ones that cross a transcript-processing pipeline boundary.
- The split-half and permutation nulls elsewhere in the repository are computed within an analysis, so they are unaffected by this sharing; it bounds how much *between-analysis* agreement should be read as independent replication.