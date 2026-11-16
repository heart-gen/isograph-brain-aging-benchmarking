# Archiving the analysis releases

The current analysis release is **v1.1.0**, accompanying the Cell Genomics
submission of
[IsoGraph separates transcript-usage and abundance contributions to aging brain networks](https://github.com/heart-gen/isograph-brain-manuscript).
Release v1.0.0 accompanied the initial submission, under the earlier title
*IsoGraph resolves coordinated transcript choice and gene abundance in the aging
human brain*, and is archived as
[10.5281/zenodo.23149044](https://doi.org/10.5281/zenodo.23149044). The
manuscript cites the all-versions DOI
[10.5281/zenodo.23149043](https://doi.org/10.5281/zenodo.23149043), which
resolves to the newest release. The data-bundle DOI is pending.

## Analysis repository archive

The repository-root [`.zenodo.json`](../.zenodo.json) describes the software
archive; [`CITATION.cff`](../CITATION.cff) records the authors and release version.
Tag annotations are `IsoGraph brain analysis -- initial manuscript submission`
(v1.0.0) and `IsoGraph brain analysis -- Cell Genomics submission` (v1.1.0).
[Release notes](../RELEASE_NOTES.md) provide the GitHub release body for each.

The GitHub–Zenodo integration is enabled. To publish a new version, push the
release commit and tag and publish the GitHub release; Zenodo deposits it as a new
version under the all-versions DOI, which the citation metadata and the companion
manuscript's availability statements cite and which does not change between
releases. Preserve earlier tags; subsequent metadata changes belong in a new
commit and, if archived, a new release.

Tracked lean results and figures use Git LFS. Before finalizing an archive,
verify that it contains artifact bytes rather than LFS pointer files. In a full
checkout, run `git lfs pull` and `git lfs fsck`; inspect the actual archive contents
and provide materialized artifacts as deposit files when needed. A source-only
archive does not include ignored inputs or heavy intermediates.

## Separate data-bundle deposit

[`zenodo/.zenodo.json`](.zenodo.json) describes the planned dataset deposit under
CC-BY-4.0. It is distinct from the root software-archive metadata. The staging
script collects these ignored heavy outputs from stages 02–07:

| Artifact | Contents |
|---|---|
| `feature_scores.parquet` | Per-sample switch/abundance feature matrix |
| `feature_reconstruction.parquet` | VAE feature reconstruction |
| `high_vs_low_table.parquet` | Per-module high-versus-low expression interpretation |
| `edges.parquet` | IsoGraph network edge lists |

Run from the repository root on the machine holding the final analysis artifacts:

```bash
bash 00_scripts/stage_zenodo_bundle.sh --checksums --tar
```

The committed `MANIFEST.tsv` currently contains paths and byte sizes, without
checksums. It is an inventory, not proof of a completed deposit. Regenerate it
against the final outputs before upload; compare its coverage with the intended
release contents. The script covers only the four artifact types above. Processed
dataset bundles and other planned public artifacts require separate inventory and
staging. Follow [`inputs/README.md`](../inputs/README.md) for source provenance;
do not include controlled-access source files in the public deposit.

The manuscript also cites the separate IsoGraph software archive
(`10.5281/zenodo.21707653`) and the independent long-read dataset
(`10.5281/zenodo.8180677`); neither DOI identifies this analysis release or its
data bundle. Add the preprint DOI to related identifiers when it is assigned.
