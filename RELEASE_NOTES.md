# v1.1.0 — IsoGraph brain analysis -- Cell Genomics submission

Release notes for the Cell Genomics submission of **IsoGraph separates
transcript-usage and abundance contributions to aging brain networks**, the
retitled manuscript. Changes since v1.0.0:

- A stricter long-read null for the colocalization-anchored and genome-wide switch
  pairs, matched jointly on the abundance quintile and the tertiles of detection
  rate, long-read isoform count and anchored usage SD (`longread_null.JOINT_KEYS`),
  reported beside the abundance-only null (06/02b). A pair anchored from both
  transcripts takes the mean of the two anchored usage SDs. Anchored distinct
  pairs: 58.1% switch-like against 47.0% under the joint null (P = 0.0015).
- `06_switch_mechanism/_h/02b.switch_orthogonal_confirm.local.sh` reruns those arms
  locally with the locked package versions; a stable sort makes the
  orientation-level draws independent of the numpy build.
- Long-read confirmation of module driver switches (06/02f).
- Fig. 6a redrawn with the joint null as the primary comparison; Data S21 gains a
  Set-level tests sheet with both nulls, and the supplementary builder writes
  extra workbook sheets.
- Supplementary tables and Data renumbered by first citation, with workbook titles
  and legends taken from the manuscript supplement.

Archive: all-versions DOI
[10.5281/zenodo.23149043](https://doi.org/10.5281/zenodo.23149043); v1.0.0 is
[10.5281/zenodo.23149044](https://doi.org/10.5281/zenodo.23149044). As for
v1.0.0, verify materialized Git LFS contents before publishing. Citation and
license terms are unchanged.

---

# v1.0.0 — IsoGraph brain analysis -- initial manuscript submission

Draft release notes for the initial manuscript submission of **IsoGraph resolves
coordinated transcript choice and gene abundance in the aging human brain**.

**Authors:** Alexis Bennett, Elisa Kain Johnson, and Kynon J. M. Benjamin.

[Manuscript source](https://github.com/heart-gen/isograph-brain-manuscript) ·
[HTML manuscript](https://heart-gen.github.io/isograph-brain-manuscript/) ·
[IsoGraph software](https://github.com/heart-gen/IsoGraph)

## Included in this release

- Synthetic benchmarking and BrainSEQ/GTEx adult brain analysis workflows.
- Module discovery, usage/abundance characterization, donor-subset stability,
  cross-cohort transfer and held-out differential transcript usage analyses.
- Transcript architecture, RBP regulation, QTL/allelic evidence, disease
  colocalization, and junction and independent long-read corroboration analyses.
- Figure/table builders, tracked result artifacts, supplementary exports,
  configurations and environment snapshots.
- Manuscript-aligned repository documentation, citation metadata and archive
  metadata, including all three authors and their ORCIDs.
- Repaired pipeline and stage-runner paths following the move to `00_scripts/`.

This release follows the manuscript's paired transcript-usage and abundance
framing. Disease analyses include the PRDM2 alternative terminal-exon event
prioritized at an amyotrophic lateral sclerosis locus. Benchmark performance
depends on the signal regime; candidate events do not establish causal mechanisms.

## Reproducibility and archive status

Use the `v1.0.0` tag with Git LFS to retrieve tracked result artifacts. Follow the
root README for environment snapshots, source-data preparation and SLURM runners.
Heavy compute requires the configured HPC resources and authorized source inputs.
The release preparation checks metadata and runner dry runs; it does not rerun
the scientific analyses.

The manuscript source used for metadata alignment is commit
`14a74b1fa6faf32c9b2c008ab467236f70a30136` of `heart-gen/isograph-brain-manuscript`.
Preprint, analysis-archive and data-bundle DOIs are pending. Large ignored artifacts
are intended for a separate deposit; the existing inventory still needs final
checksums and coverage verification. Verify materialized Git LFS contents before
publishing an archive. See [archiving instructions](zenodo/README.md).

## Citation and license

Cite this release using [CITATION.cff](CITATION.cff) and cite the companion
manuscript. The separate IsoGraph package archive cited by the manuscript is
[10.5281/zenodo.21707653](https://doi.org/10.5281/zenodo.21707653).

Original code: [Apache-2.0](LICENSE). Original derived data, figures, tables and
documentation: [CC-BY-4.0](LICENSE-DATA.md). Third-party resources retain their
original terms; controlled-access source data are not redistributed.
