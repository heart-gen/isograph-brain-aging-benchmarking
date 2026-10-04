# Junction-pair corroboration

**Pair-discriminating short-read junction support for disease-prioritized transcript pairs.** Each small point is one exact nominated transcript pair; large diamonds show gene–region medians and horizontal lines span the observed pair range, not a confidence interval. The x axis is the fraction of adult control donors with at least one fragment assigned to each alternative, among donors with at least 10 total pair-discriminating fragments. Pairs require at least 30 covered donors. Right-hand labels give measured/targeted pairs and the range of covered donor counts. Of 55 signal-colocalization-prioritized genes, 26 have a primary BrainSEQ regional counterpart; 30 pair–region combinations in 7 genes are measurable with existing counts. Remaining regional nominations are retained in the coverage ledger rather than treated as negative. Primary mappings include BA9–DLPFC, hippocampus–hippocampus and generic cortex–DLPFC (an anatomical proxy); adjacent and secondary regions are separate sensitivities. Counts are pooled across haplotypes including unassigned fragments, without allelic-significance or heterozygosity selection. Support is for the pair-discriminating junction models, not necessarily direct observation of the particular colocalizing junction or full-length transcript. The measurements reuse BrainSEQ reads and are complementary corroboration, not independent phenotype replication.

```json
{
  "n_nominated_genes": 55,
  "n_nominating_events": 115,
  "n_primary_region_genes": 26,
  "n_primary_pair_regions": 260,
  "n_structurally_eligible_primary_pair_regions": 31,
  "n_measured_genes": 7,
  "n_measured_pair_regions": 30,
  "n_measured_distinct_pairs": 30,
  "both_form_donor_threshold": 30,
  "n_pairs_both_in_min_donors": 29,
  "fraction_both_min": 0.11312217194570136,
  "fraction_both_max": 1.0,
  "covered_donors_min": 220,
  "covered_donors_max": 237,
  "n_descriptive_concordances": 30,
  "median_descriptive_adjusted_rho": 0.3152587448955686,
  "sample_membership": [
    {
      "region": "dlpfc",
      "n_bundle_donors": 222,
      "n_count_completed": 221
    },
    {
      "region": "hippocampus",
      "n_bundle_donors": 238,
      "n_count_completed": 237
    },
    {
      "region": "caudate",
      "n_bundle_donors": 238,
      "n_count_completed": 237
    }
  ],
  "primary_status_counts": {
    "pair_not_counted": 224,
    "measured": 30,
    "not_two_sided": 5,
    "insufficient_coverage": 1
  },
  "nominated_genes_without_primary_region": [
    "ASB3",
    "C9orf72",
    "CCDC62",
    "CD46",
    "CDHR3",
    "DNAJA3",
    "DOC2A",
    "FGFR1",
    "GPM6A",
    "GPR135",
    "IFNAR2",
    "INO80E",
    "NCOR1",
    "NDUFAF7",
    "NSMAF",
    "NUCB2",
    "PAK6",
    "PCBP3",
    "PPIL2",
    "PPIP5K1",
    "RPS6KL1",
    "SCFD1",
    "SLC39A13",
    "SNCA",
    "SPI1",
    "TAOK2",
    "THAP3",
    "TPCN1",
    "ZDHHC12"
  ],
  "primary_region_genes_without_measured_pairs": [
    "ACTR1B",
    "CDIP1",
    "CRELD2",
    "FNBP1",
    "HMOX2",
    "IDH3B",
    "INTS8",
    "ITGB1BP1",
    "KLC1",
    "MAD1L1",
    "MED19",
    "NMRAL1",
    "NUP50",
    "PTPRN",
    "SIRPA",
    "SNAP91",
    "TMEM175",
    "YPEL1",
    "ZSWIM7"
  ]
}
```
