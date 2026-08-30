# EGA access request for neuronal eCLIP datasets

These controlled datasets are an optional extension. They are not required for the
public-data primary analysis and do not block completion of neuronal CLIP validation.

## Datasets

If pursuing the optional extension, request access separately for both datasets in study
[EGAS00001005880](https://ega-archive.org/studies/EGAS00001005880):

- [EGAD00001008426](https://ega-archive.org/datasets/EGAD00001008426):
  TDP-43 eCLIP from iPSC-derived motor neurons in two control lines; four samples.
- [EGAD00001008428](https://ega-archive.org/datasets/EGAD00001008428):
  NOVA1, NOVA2 and RBFOX2 eCLIP from iPSC-derived motor neurons in two control
  lines; twelve samples.

Both datasets are controlled by `EGAC00001002455`, the Data Access Committee of
the Department of Molecular Neurology, UKEr. The listed contact is Martin
Regensburger (`martin.regensburger@uk-erlangen.de`).

## Request procedure

1. [Register for an EGA account](https://www.ega-archive.org/register/) using an
   institutional email address. EGA states that account validation may take two
   working days.
2. After validation, log in to EGA and open each dataset page linked above.
3. Select **Request Access** on each page. Requests are dataset-specific, so submit
   one request for `EGAD00001008426` and another for `EGAD00001008428`.
4. Download the Data Access Agreement when prompted, complete and sign it, and
   upload it with the request.
5. Include the project title, scientific purpose, principal investigator,
   institution, authorized personnel, storage and access-control plan, intended
   outputs, and publication plans.
6. Retain the EGA request identifiers. If the requests remain unanswered, contact
   the DAC and cite both dataset accessions and the request identifiers.

EGA accounts are individual. Do not share credentials; every person who will
download or directly access the controlled files should be named as required by the
Data Access Agreement and have their own account. EGA explains its distributed
access model and Data Access Agreement process in its
[DAC guidance](https://ega-archive.org/access/data-access-committee/what-is-dac/).

## Suggested project title

**Neuronal eCLIP validation of isoform-switch networks in human brain aging and
schizophrenia**

## Suggested project description

> We will use processed eCLIP peak and signal files from control human
> iPSC-derived motor neurons to validate prespecified RNA-binding-protein
> regulatory hypotheses generated independently from human brain transcript-usage
> data. We will test whether TDP-43, NOVA1, NOVA2 and RBFOX2 binding intersects
> transcript-switch-relevant exonic and intronic regions more frequently than
> matched within-gene control regions. Results will be reported only as aggregate
> enrichment statistics and non-identifying genomic illustrations. We will not
> attempt participant identification or redistribute individual-level data. Data
> will be stored on access-controlled institutional HPC storage and handled
> according to the Data Access Agreement.

## Information to prepare

- Principal investigator and institutional affiliation
- Project title and focused scientific rationale
- Names and roles of all authorized personnel
- Institutional HPC storage location and access-control policy
- Encryption, logging, backup and incident-response practices, if requested
- Confirmation that participant re-identification will not be attempted
- Confirmation that controlled files will not be redistributed
- Intended aggregate outputs and publication plan
- Project duration and planned data-destruction or retention procedure

## After approval

1. Review the approved Data Access Agreement and record its restrictions and
   expiration date.
2. Download files through an EGA-approved mechanism using the authorized user's
   individual credentials.
3. Never place EGA credentials or tokens in the repository, YAML configuration,
   SLURM script or job log.
4. Store controlled files under an access-restricted directory within
   `inputs/raw/neuronal_clip/EGA/`; do not commit them.
5. Extend `configs/neuronal_clip.yaml` with the authorized EGAF file inventory,
   sizes and checksums before analysis.
6. Regenerate `07_rbp_regulation/_m/neuronal_clip_manifests/dataset_manifest.tsv` so the controlled-access
   state is replaced with sample- and file-level provenance.

The reproducible public-data acquisition wrapper is
`inputs/_h/download_neuronal_clip.sh`. It intentionally records EGA datasets as
controlled access and does not attempt to bypass authentication.
