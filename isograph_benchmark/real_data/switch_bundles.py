"""Named bundles of IsoGraph analyses whose switch layers are pooled together.

The genetic-anchoring capstone treats two cases symmetrically:

  * ``disease``  — the single SCZD analysis (brainseq-sczd), anchored in SCZ GWAS.
  * ``aging``    — the pooled aging switch layer across all 16 aging analyses
                   (3 BrainSEQ aging regions + 13 GTEx brain tissues), anchored in
                   neurodegeneration GWAS (AD, PD, LBD, ALS[, FTD]).

Pooling the aging analyses gives one aging switch-gene set with enough QTL SNPs for a
well-powered S-LDSC heritability partition and enough loci for coloc, mirroring how the
disease case draws on the one SCZD analysis. Each entry is an ``(analysis, region)``
pair understood by :func:`isograph_benchmark.real_data.sweep_leiden._artifact_dir`.
Mirrors the analysis grid in ``02_module_discovery/brainseq/_h/13.qtl_anchoring.sh``.
"""
from __future__ import annotations

AGING_BUNDLE: list[tuple[str, str | None]] = [
    ("brainseq-aging", "caudate"),
    ("brainseq-aging", "hippocampus"),
    ("brainseq-aging", "dlpfc"),
    ("gtex-aging", "amygdala"),
    ("gtex-aging", "anterior_cingulate_cortex_ba24"),
    ("gtex-aging", "caudate_basal_ganglia"),
    ("gtex-aging", "cerebellar_hemisphere"),
    ("gtex-aging", "cerebellum"),
    ("gtex-aging", "cortex"),
    ("gtex-aging", "frontal_cortex_ba9"),
    ("gtex-aging", "hippocampus"),
    ("gtex-aging", "hypothalamus"),
    ("gtex-aging", "nucleus_accumbens_basal_ganglia"),
    ("gtex-aging", "putamen_basal_ganglia"),
    ("gtex-aging", "spinal_cord_cervical_c_1"),
    ("gtex-aging", "substantia_nigra"),
]

DISEASE_BUNDLE: list[tuple[str, str | None]] = [
    ("brainseq-sczd", None),
]

BUNDLES: dict[str, list[tuple[str, str | None]]] = {
    "aging": AGING_BUNDLE,
    "disease": DISEASE_BUNDLE,
}


def get_bundle(name: str) -> list[tuple[str, str | None]]:
    if name not in BUNDLES:
        raise SystemExit(f"Unknown bundle {name!r}; known: {sorted(BUNDLES)}")
    return BUNDLES[name]
