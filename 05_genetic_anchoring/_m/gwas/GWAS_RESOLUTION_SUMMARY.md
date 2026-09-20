# Module GWAS enrichment is size-confounded at every resolution, and worse in the baseline

Modular analysis summary for Manubot integration. **Hand-written** — nothing regenerates this
file, so it must be re-quoted by hand against
`05_genetic_anchoring/_m/gwas/magma_results_combined.parquet` (production, Leiden 2.0) and
`magma_results_combined_res5.parquet` (the resolution-5.0 sensitivity arm) whenever those
change. *(Re-quoted 2026-09-19 against the switching-filter re-run; it had carried the legacy
expression-filter numbers and the opposite conclusion.)*

## Purpose

Stop a false-positive enrichment story before it reaches the manuscript. MAGMA competitive
gene-set analysis rewards large gene sets, so any module pipeline that emits a few giant
modules will show "significant GWAS enrichment" that is an artifact of module size, not
biology. This analysis asks (1) whether IsoGraph's module–GWAS hits are driven by giant
modules, (2) whether Leiden resolution can remove that artifact, and (3) whether a real,
size-controlled disease signal survives — set against the gene-level WGCNA baseline, which is
itself an unguarded giant-module pipeline. The answer to (2) on the switching filter is **no**
in both directions tested, which is why the result is reported rather than engineered away.

## Inputs

- **Module gene sets** — IsoGraph `isograph_vae` modules at production resolution 2.0 and, in
  parallel, the resolution-5.0 sensitivity arm (`isograph_vae_res5`); plus classical
  `wgcna_gene` modules (resolution-independent), across BrainSEQ and 13 GTEx regions.
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
across module definitions, so the only thing that varies between IsoGraph resolution 2.0
(production), IsoGraph resolution 5.0 (sensitivity) and gene-level WGCNA is the module
partition. Module-level
p-values were Benjamini–Hochberg corrected, and a module was called significant at FDR <
0.05. We classified each module by size and flagged "giant" modules at ≥ 900 genes — the
size class above which MAGMA's competitive test is dominated by set size rather than
trait-specific signal. The artifact test is comparative: a size-controlled pipeline should
not concentrate its significant hits in giant modules. Analyses used MAGMA v1.10 with the
project Python 3.12 environment for input preparation and deterministic module fits.

## Results text

**Giant modules carry most of the significant hits at both resolutions.** On the switching
filter, production resolution 2.0 gives **37** FDR-significant IsoGraph module × trait hits,
**24 of them (65%)** in modules of ≥ 900 genes (median significant set 1,182 genes). The
resolution-5.0 sensitivity arm is **worse**, not better: **26** significant hits of which
**20 (77%)** are giant (median 1,390). Raising the resolution does not remove the size
dependence — it removes small modules' hits faster than giant ones'.

**This reverses the legacy finding, and the reversal is what moved production to 2.0.** On
the legacy expression filter the same comparison read 17/28 giant at resolution 2.0 against
0/10 at 5.0, and that contrast is what selected 5.0 as canonical. Under the switching
transcript filter the criterion points the other way, which — together with ~38% more gene
coverage — is why the PI moved production to resolution 2.0 on 2026-09-16.

**The gene-level WGCNA baseline is a far worse giant-module pipeline.** Classical
`wgcna_gene` returns **110** FDR-significant hits, **83 (75%)** of them giant, with a median
significant set of **9,906** genes and the largest spanning **12,195**. IsoGraph's median
significant set is eight times smaller. The honest comparison is therefore not "IsoGraph
avoids the artifact" but "every module pipeline on these data shows it, and the matched
abundance baseline shows it on sets an order of magnitude larger".

**Which traits survive.** At resolution 2.0 the IsoGraph hits are SCZ 20 (15 giant), BP 11
(7), MDD 3 (1), ALS 2 (1) and PD 1 (0); AD, LBD and stroke give none. The 13 non-giant hits
are spread across independent cohorts and regions — BrainSEQ caudate M003 (836 genes, SCZ,
FDR 0.004), BrainSEQ DLPFC M000 (806, SCZ, 0.005), BrainSEQ hippocampus M001 (574, MDD,
0.004), GTEx spinal cord M003 (319, BP, 0.004), GTEx cerebellum M002 (513, ALS, 0.006) and
GTEx spinal cord M006 (206, PD, 0.019) among them — so a size-controlled signal exists, it is
simply not the majority of the hits. *(Module ids are re-assigned at every fit and are valid
only against this partition.)*

**Headline:** *MAGMA module-GWAS enrichment is confounded by module size at every resolution
tested: 65% of significant IsoGraph hits at production resolution 2.0 and 77% at resolution
5.0 fall in ≥ 900-gene modules, against 75% for gene-level WGCNA on sets with a median of
9,906 genes. The giant-module effect is reported and compared rather than engineered away
(PI, 2026-09-14); the defensible claim is the minority of size-controlled hits, led by SCZ
and BP across independent regions.*

## Figure and table notes

- **Supplementary figure — GWAS resolution
  (`manuscript/_m/figures/figGwasResolution.{pdf,png}`, built by
  `manuscript/_h/gwas_resolution_figure.R`).** (A) MAGMA −log10 P vs module size (log) for
  IsoGraph production res 2.0, the res-5.0 sensitivity arm and WGCNA, significant hits
  highlighted, with the 900-gene giant cutoff marked. (B) significant-module counts split by
  giant vs non-giant per method — the giant fraction is high everywhere and highest in the
  res-5.0 arm. (C) the size-controlled (< 900-gene) hits as a dot plot (−log10 FDR, sized by
  module genes, coloured by trait). No in-panel titles; the read lives in the caption.
- **Tables:** `magma_results_combined.parquet` (production, res 2.0; `magma_results_combined_res2.parquet`
  is an identical copy under the explicit name) and `magma_results_combined_res5.parquet`
  (sensitivity arm) — per module × trait × backend MAGMA NGENES, BETA, SE, P, BH FDR.
