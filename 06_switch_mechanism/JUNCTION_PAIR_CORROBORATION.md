# Local junction corroboration for Figure 6c

Run from the repository root:

```bash
bash 06_switch_mechanism/_h/02g.junction_pair_corroboration.sh
```

The wrapper uses `$HOME/.venvs/isograph/bin/python` if present, otherwise `python3`.
Override with `PYTHON_BIN=/path/to/python RSCRIPT_BIN=/path/to/Rscript`.
Python dependencies are numpy, pandas, pyarrow and scipy; the figure needs R packages
ggplot2, dplyr and jsonlite, and Cairo PDF support. No model fitting, GPU, network,
BAM or scheduler access is required. Inputs must be actual local files, not Git LFS pointers.

The standalone module has `--root`, `--output-dir`, `--min-total-fragments`,
`--min-donors`, `--depth-sensitivity` and `--batch-size` options. Custom output
directories can be plotted using:

```bash
Rscript manuscript/_h/junction_pair_corroboration_figure.R /path/to/repo /path/to/results /path/to/figures
```

## Scientific design fixed before execution

Use the signal-colocalization gene nominations underlying Figure 6a–b. Recover
actual transcript-pair edges from each nominated GTEx tissue's
`structure_switch_pairs.parquet`. Do not create a clique from the event table's
flattened list of transcript members. Match versioned transcript IDs exactly,
allowing pair reversal but not silently stripping versions. Retain all nomination
events and their mapping statuses in the event ledger.

Primary anatomical mappings are BA9–DLPFC, hippocampus–hippocampus,
caudate–caudate and generic cortex–DLPFC. The last is explicitly a cortex proxy.
Strict anatomy, adjacent mappings, secondary caudate mappings and pairs for which
the anchored junction is discriminating are separate sensitivities. A matched
recount pair must have two-sided discriminating junction models and carry the
anchored junction at exact coordinates after converting LeafCutter's exon-flank
coordinates to 1-based closed introns (`start + 1`, `end - 1`), with strand checked.

Adult control donors must belong to the aging bundle and have completed count QC.
Missing count runs are not zeros. Sum fragment counts over haplotypes 0, 1 and 2;
hap 0 contributes isoform evidence without allele assignment. Do not filter on
allelic significance, heterozygosity or the ASE screen's testability flag.

Primary endpoint: among donors with at least 10 pair-discriminating fragments,
the proportion with at least one fragment for each alternative. Require at least
30 covered donors for a measurable pair. Sensitivities use 5 and 20 total
fragments and at least two fragments per alternative. Report minor-form usage
continuously, without a 5% biological validity gate. Several pairs per gene are
dependent; the figure shows their values and range, not an inferential interval.
Set summaries, where supplied, weight genes equally.

Secondary descriptive concordance compares oriented junction fraction with the
original transcript estimates normalized by the available EffectiveLength field.
Calculate raw Spearman and correlation of rank residuals adjusted for age (centered
cubic), sex, manner of death, RIN, mapping rate, mitochondrial rate and five genotype
PCs. Complete-case donor counts and non-estimable statuses are explicit. No nominal
P values or replication claims are attached to these shared-read correlations.
No composition-adjusted independent-validation claim is made.

The selected junction models discriminate the two members of a pair, not necessarily
each transcript against every other transcript of the gene. Counts aggregate over
the model's junctions: they do not by themselves prove the exact colocalizing
junction was observed. BrainSEQ can overlap discovery, so this is complementary
measurement corroboration, not independent disease-effect replication.

## Inputs, outputs and reproducibility

Inputs: signal event ledger; GTEx tissue-specific structural pair ledgers;
BrainSEQ regional recount pair/junction models, count QC and fragment tables;
adult-control sample tables; original transcript estimates. The analysis reads
large count tables in batches, one region at a time, retaining only target pairs
and eligible samples. It does not modify any existing analysis outputs.

Outputs in `06_switch_mechanism/_m/junction_pair_corroboration/`:

- `params.json`: design and thresholds, written before computing outcomes.
- `target_event_pairs.*`: complete nomination, pair and region audit.
- `donor_pair_counts.*`: selected donor-level counts and concordance inputs.
- `primary_pair_results.*`, `primary_gene_region_summary.*`: panel data.
- `pair_results_all_thresholds.*`, `sensitivity_summary.*`: all sensitivity results.
- `gene_coverage.*`, `primary_coverage.*`, `summary.json`: explicit denominators.
- `FIGURE_CAPTION.md`, `RESULTS.md`: generated interpretation and result counts.
- `provenance.json`: full SHA-256 hashes of inputs/scripts, package versions,
  parameters, Git revision, dirty status and elapsed time.
- `plot_sessionInfo.txt`: plotting environment.

The plotter writes vector PDF and 400-dpi PNG to
`manuscript/_m/figures/figJunctionPairCorroboration_c.{pdf,png}`. It uses an 82-mm
column, the existing Figure 6 typography, an accessible palette, no internal title,
and no panel tag (matching the individual a/b exports for later assembly).

Validation:

```bash
PYTHONPATH=. "$HOME/.venvs/isograph/bin/python" -B -m unittest discover -s tests -p test_junction_pair_corroboration.py
```

Fixtures guard against invented pair edges, coordinate/strand mismatches, exclusion
of hap 0, accidental inclusion of cases, treating missing samples as zeros, and
reporting covariate-only concordance as residual evidence. Re-running with the same
inputs and parameters is deterministic; elapsed time and environment/Git metadata
can change. Keep the generated provenance with a frozen result release.
