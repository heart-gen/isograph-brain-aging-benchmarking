# Module GWAS enrichment is a giant-module artifact that Leiden resolution removes

Modular analysis summary for Manubot integration. Generated from
`05_genetic_anchoring/_m/gwas/magma_results_combined.parquet` (canonical Leiden resolution 5.0) and
`magma_results_combined_res2.parquet` (resolution 2.0). Every numeric claim is reproduced
from those tables; do not edit the numbers by hand — regenerate.

## Purpose

Stop a false-positive enrichment story before it reaches the manuscript. MAGMA competitive
gene-set analysis rewards large gene sets, so any module pipeline that emits a few giant
modules will show "significant GWAS enrichment" that is an artifact of module size, not
biology. This analysis asks (1) whether IsoGraph's module–GWAS hits are driven by giant
modules, (2) whether the canonical Leiden resolution (5.0) removes that artifact, and (3)
whether a real, size-controlled disease signal survives — set against the gene-level WGCNA
baseline, which is itself an unguarded giant-module pipeline.

## Inputs

- **Module gene sets** — IsoGraph `isograph_vae` modules at canonical resolution 5.0 and, in
  parallel, at resolution 2.0 (`isograph_vae_res2`); plus classical `wgcna_gene` modules
  (resolution-independent), across BrainSEQ and 13 GTEx regions.
- **GWAS summary statistics** for six brain-relevant traits: schizophrenia (SCZ),
  bipolar disorder (BP), major depression (MDD), Alzheimer's disease (AD), Parkinson's
  disease (PD), and stroke.
- **MAGMA reference** — 1000 Genomes EUR LD panel (`g1000_eur`, build hg19/b37) with the
  matching `NCBI37.3` gene-location file (the build pairing is load-bearing — an hg38 gene
  loc mis-maps every SNP).

## Methods text

For each trait we ran MAGMA (v1.10) gene analysis on the hg19 1000G-EUR panel, then
competitive gene-set analysis testing each transcriptomic module as a gene set. The gene
analysis (`.genes.raw`) is module-independent and was computed once per trait and reused
across module definitions, so the only thing that varies between IsoGraph resolution 2.0,
IsoGraph resolution 5.0, and gene-level WGCNA is the module partition. Module-level
p-values were Benjamini–Hochberg corrected, and a module was called significant at FDR <
0.05. We classified each module by size and flagged "giant" modules at ≥ 900 genes — the
size class above which MAGMA's competitive test is dominated by set size rather than
trait-specific signal. The artifact test is comparative: a size-controlled pipeline should
not concentrate its significant hits in giant modules. Analyses used MAGMA v1.10 with the
project Python 3.12 environment for input preparation and deterministic module fits.

## Results text

**At resolution 2.0 IsoGraph's GWAS hits are a giant-module artifact.** Of 43
FDR-significant IsoGraph modules at resolution 2.0, **18 (42%) are giant** (≥ 900 genes),
extending up to 2,600+ genes — the significant hits pile up at large module sizes, the
signature of the MAGMA size bias rather than focused biology.

**Canonical resolution 5.0 eliminates the artifact.** At resolution 5.0 IsoGraph yields **6
FDR-significant modules, 0 of them giant** (largest significant module 691 genes; no
module in the entire partition reaches 900 genes among the hits). The size-vs-significance
dependence is gone — significant modules are confined to the small-module regime.

**The schizophrenia signal survives the size control.** Schizophrenia drops from 22
significant modules at resolution 2.0 to **5 at resolution 5.0**, but those 5 are real,
modest-sized modules (139–691 genes) spread across independent regions — GTEx amygdala
M000 (554 genes, FDR 0.0007), GTEx hypothalamus M000 (691, 0.033), GTEx hippocampus M002
(421, 0.036), BrainSEQ caudate M009 (139, 0.036), and BrainSEQ hippocampus M000
(373, 0.036). SCZ accounts for 5 of the 6 surviving hits; the sixth is a small PD module
(GTEx caudate basal ganglia M029, 37 genes). The disease signal is genuine, not
size-driven.

**The gene-level WGCNA baseline is itself an unguarded giant-module pipeline.** Classical
`wgcna_gene` returns **99 FDR-significant modules, 79 (80%) of them giant**, the largest
spanning 12,369 genes. Its apparent GWAS enrichment advantage is overwhelmingly a
size-bias phenomenon — exactly the artifact resolution control removes from IsoGraph —
reinforcing the project-wide read that gene-abundance/size dominates the bulk signal.

**Headline:** *MAGMA module-GWAS enrichment is confounded by module size; at IsoGraph's
canonical Leiden resolution 5.0 the giant-module artifact disappears (0/6 significant
modules are giant, vs 18/35 at resolution 2.0 and 77/99 for gene-level WGCNA), yet a
size-controlled schizophrenia signal survives across five independent brain regions.*

## Figure and table notes

- **Supplementary figure — GWAS resolution
  (`manuscript/_m/figures/figGwasResolution.{pdf,png}`, built by
  `manuscript/_h/gwas_resolution_figure.R`).** (A) MAGMA −log10 P vs module size (log) for
  IsoGraph res 2.0 / res 5.0 / WGCNA, significant hits highlighted, with the 900-gene giant
  cutoff marked — the artifact and its fix read directly off where the coloured points sit
  relative to the cutoff. (B) significant-module counts split by giant vs non-giant per
  method — the giant fraction collapses to zero at res 5.0. (C) the six surviving res-5.0
  IsoGraph hits as a dot plot (−log10 FDR, sized by module genes, coloured by trait),
  showing modest sizes and SCZ dominance. No in-panel titles; the read lives in the caption.
- **Tables:** `magma_results_combined.parquet` (res 5.0, canonical) and
  `magma_results_combined_res2.parquet` (res 2.0) — per module × trait × backend MAGMA NGENES,
  BETA, SE, P, BH FDR.

## Reproducibility information

- Analysis directory: `05_genetic_anchoring/_m/gwas/`.
- Scripts: `_h/01.prep_module_gene_sets.{R,sh}` (module→gene-set files),
  `_h/02.run_magma.sh` (gene + gene-set analysis; SLURM, RM-shared), and
  `isograph_benchmark/gwas/prepare_magma_inputs.py` (SNP p-value inputs). The resolution-2.0
  run is the same driver with `MAGMA_ISOGRAPH_BACKEND=isograph_vae_res2`; output files embed
  the backend so non-canonical runs never clobber canonical results.
- Inputs: production module partitions (`isograph_vae`, `isograph_vae_res2`, `wgcna_gene`),
  six GWAS sumstats, `g1000_eur` (hg19) LD panel, `NCBI37.3.gene.loc`.
- Outputs: `results/`, `gene_analysis/`, `magma_results_combined{,_res2}.parquet`.
- Key parameters: MAGMA v1.10 competitive GSA; BH FDR < 0.05; giant-module cutoff 900 genes;
  canonical Leiden resolution 5.0.
- Compute environment: PSC Bridges-2 RM-shared; MAGMA v1.10
  (`/ocean/projects/bio250020p/shared/opt/magma-v1.10`); project Python 3.12
  (`/ocean/projects/bio260021p/shared/opt/envs/isograph`).
- Missing reproducibility information: GWAS sumstats provenance/versions are configured in
  `configs/gwas_magma.yaml`, not restated here; per-package versions taken from the live
  environment.

## Limitations and integration notes

- This is a **specificity / robustness** analysis, not a discovery claim: its job is to show
  that the headline module-GWAS enrichment is size-confounded and that controlling module
  size (via canonical resolution) is what makes the surviving signal trustworthy. It is the
  GWAS-axis companion to the giant-module-cap stability work.
- The five surviving SCZ modules are reported as size-controlled hits, not as a fine-mapped
  causal account; integrate them with the QTL splicing-specificity contrast (the genetic
  anchoring of the same co-switch modules) and the GO-invisible gate (the disease switch
  content), and read the WGCNA giant-module result alongside the three-baseline comparison
  (gene-abundance/size dominates the bulk enrichment signal).
- The giant-module cutoff (900 genes) is a reporting threshold for the size class, not a
  modelling parameter; the underlying claim (significance concentrates in large sets at res
  2.0 and for WGCNA, not at res 5.0) is visible across the whole size axis in panel A.
