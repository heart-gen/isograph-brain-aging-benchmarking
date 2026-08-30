# Session next steps + per-gene mechanistic deep-dive plan

_Working doc. Login node was slow (SLURM controller timing out on `squeue`/`sbatch`);
move to a compute/interactive node to resume. Not committed — delete or keep as you like._

Branch: `clinical-consequence-light` (2 commits, unpushed: `d56e30c`, `f84e4e7`).

---

## PART A — Operational resume (do these first, on a compute node)

### A0. Get an interactive node
```bash
interact -A bio260021p -p RM-shared --cpus-per-task=2 -t 2:00:00
cd /ocean/projects/bio260021p/kbenjamin/projects/isograph-brain-aging-benchmarking
```

### A1. Clinical-consequence array — finish the 4 failed regions
**State:** job `42404471` (pre-module-guard-fix) ran 17 tasks:
- 7 wrote outputs: brainseq/caudate, gtex/{cerebellum, cortex, frontal_cortex_ba9,
  hippocampus, hypothalamus, substantia_nigra}
- 6 legitimately **skipped** (no phenotype-significant switch genes): brainseq/hippocampus,
  gtex/{caudate_basal_ganglia, cerebellar_hemisphere, nucleus_accumbens_basal_ganglia,
  putamen_basal_ganglia, spinal_cord_cervical_c_1}
- **4 FAILED exit 127** (`module: command not found`, now fixed by the module-init guard):
  indices **0 (brainseq caudate_sczd), 3 (brainseq dlpfc), 4 (gtex amygdala),
  5 (gtex anterior_cingulate_cortex_ba24)**

A re-run of just those 4 was attempted on the login node but `sbatch` timed out with **no
jobid printed and no jobid file written** → almost certainly did NOT land. **Verify first,
then submit if clear:**
```bash
squeue -u kbenjamin -n clinical-consequence -o '%i %t %r'   # expect empty
# if empty:
sbatch --array=0,3,4,5 02_module_discovery/brainseq/_h/18.clinical_consequence.sh
```
Confirm success = each of the 4 either writes `.../clinical_consequence/clinical_consequence.parquet`
OR logs a legitimate "skipping" reason (caudate_sczd should NOT skip — it's the SCZD region).

### A2. Commit the wrapper fixes (USER ALREADY APPROVED, gated on A1 success)
Once the re-run completes cleanly, commit the module-init guard added to both wrappers:
```bash
git add 02_module_discovery/brainseq/_h/15.switch_consequence.sh \
        02_module_discovery/brainseq/_h/18.clinical_consequence.sh
git commit   # message below
```
Commit message:
```
Guard module-init in switch/clinical-consequence SLURM wrappers

Compute nodes run the array body in a non-interactive shell where the `module`
function isn't pre-loaded, so `module purge` aborted every task under `set -e`
(exit 127). Source the lmod init before use in both wrappers.
```
_Selective add only — never `git add -A`. `real_data/coloc/_m` and `real_data/ldsc/_m`
stay untracked._

### A3. Cross-region meta rollup (after all 17 accounted for)
```bash
python -m isograph_benchmark.real_data.clinical_consequence_meta
# writes real_data/_m/clinical_consequence_meta.parquet + CLINICAL_CONSEQUENCE_META.md
```
Report: median LOEUF + LOEUF Fisher p (PRIMARY anchor), ClinVar ratio>1/<1 counts,
median ratio, exon-contrast Fisher p — split by stratum (all/go_invisible/go_visible)
× scope (all_exons/cds).

---

## PART B — Per-gene mechanistic deep-dive (Task #1)

**Goal:** move beyond the general co-switch-gene review to per-gene mechanistic vignettes
for the colocalization-prioritized disease genes — integrate every layer we've built into a
single causal story per gene, and pick 2 for main-figure treatment.

### B1. The gene panel (grouped by anchoring trait)
| Gene | Trait(s) | Coloc signal (from genetic_anchoring_v2) |
|------|----------|-------------------------------------------|
| **SNCA** | LBD, PD | splicing coloc (candidate main figure) |
| TPP1 | ALS | splicing coloc |
| SCFD1 | ALS | splicing coloc |
| PGS1 | ALS | |
| PPP6R2 | ALS | |
| GGNBP2 | ALS | |
| CTSH | AD | |
| TPCN1 | AD | |
| MRPS10 | AD | |
| **MYO18A** | SCZ | (candidate main figure) |
| PBX1 | SCZ | |
| RBFA | SCZ | |
| MED15 | SCZ | |

_Verify this list against the coloc registry on the compute node before writing anything —
paths/PP4/lead-variant to pull:_ `real_data/coloc/**` lean summaries + the
`switch_bundles` / `gwas_traits` registries referenced in `project_genetic_anchoring_v2`.

### B2. The six evidence layers to join per gene (all keyed by `gene_id`)
For each gene, assemble a row pulling from artifacts already on disk:

1. **Genetic anchor** — `real_data/coloc/`: which QTL colocalizes (sQTL = splicing-led,
   the IsoGraph-unique case; vs eQTL = expression-led), PP4, lead variant + position,
   direction of effect, tissue.
2. **The switch** — IsoGraph module membership + isoform event: which module, which
   region(s), which transcripts (dominant→alternative), switch importance / loading.
   Source: per-region `isograph_vae` module tables + `isoform_events`.
3. **Coding consequence** — `switch_consequence.parquet`: structural class of the switch
   (UTR remodel / CDS remodel / NMD-routing / biotype change) and the within-gene-null
   verdict. Ties to the SWITCH_CONSEQUENCE_SUMMARY headline (productive UTR/CDS remodel,
   not decay).
4. **Regulatory logic** — `real_data/_m/rbp/`: which RBP motifs are enriched in the
   gene's module (KHDRBS1, PPIE, A1CF, U2AF2, …). **Key mechanistic test:** does the
   colocalizing sQTL lead variant fall in/near an enriched RBP motif in a switched exon?
   That is the concrete "variant → altered RBP binding → splice switch" hypothesis.
5. **Constraint / clinical** — `clinical_consequence.parquet` + `gene_constraint.parquet`:
   gene LOEUF (constrained?), any ClinVar P/LP in the switched exons (CDS scope).
6. **Literature** — known isoform biology for the gene in the disease (needs WebSearch on
   compute node w/ network, or manual): e.g., SNCA 3'UTR/aSyn isoforms in synucleinopathy.

### B3. Per-gene verdict axis
Classify each gene as:
- **splicing-led, GO-invisible** (sQTL coloc + switch drives DTU without DGE) → IsoGraph's
  unique contribution, the headline cases; or
- **expression-confounded** (eQTL coloc dominates) → honest limitation, report as such.

### B4. Proposed implementation — one joining CLI
Add `isograph_benchmark/real_data/gene_deep_dive.py`:
- input: `--gene SYMBOL` (or `--genes SNCA,MYO18A,…`), resolves symbol→`gene_id`;
- joins the six layers above across all regions where the gene switches;
- emits per-gene `real_data/_m/deep_dive/<GENE>.md` vignette + a panel-wide
  `deep_dive_panel.parquet` (one row/gene, the B2 columns) for the summary table;
- deterministic, no new heavy compute — pure joins over existing parquets.
- SLURM: trivial, single `RM-shared` task (or run interactively).

### B5. Main-figure candidates
- **SNCA** — strongest a priori (LBD+PD splicing coloc, canonical synucleinopathy gene,
  well-documented isoform biology) → mechanistic schematic: lead sQTL → RBP motif →
  exon/3'UTR switch → module → constraint.
- **MYO18A** — SCZ splicing coloc; contrast case in a different disease axis.

### B6. Deliverables
1. `deep_dive_panel.parquet` + a summary table (all 13 genes × 6 layers + verdict).
2. Two figure-quality vignettes (SNCA, MYO18A).
3. Fold the splicing-led cases into the manuscript synthesis (Manubot), alongside
   SWITCH_CONSEQUENCE_SUMMARY / RBP_REGULON_SUMMARY / (pending) CLINICAL_CONSEQUENCE_META.

---

## PART C — Remaining backlog (after A + B)
- **Manuscript synthesis** (medium): fold clinical-consequence meta + deep-dive vignettes
  into the Manubot summaries.
- **Push/PR** `clinical-consequence-light` → main (like PR #7), once A2 + deep-dive land.
- **RBP full genomic-intronic scope** (large): genome FASTA staged at
  `/ocean/projects/bio260021p/shared/resources/genomes/human/gencode-v47/fasta`.
- **External-cohort replication** (large, blocked on data): project frozen modules into
  independent disease cohorts.
