# IsoGraph candidate-switch follow-up plan

## Goal

Test whether the prioritized **SNCA, CTSH, and PPP6R2 transcript switches** are detectable in the disease-relevant human brain context, determine which cell type carries each switch, and establish whether the associated genetic signal changes the switch fraction independently of total gene abundance.

This plan deliberately separates three levels of evidence:

1. **Switch validation:** confirm the exact junction or transcript pair.
2. **Cell-type localization:** identify the cell population in which the switch occurs.
3. **Mechanism:** test allelic regulation and, only after validation, candidate RBP regulators.

The immediate paper-facing goal is levels 1 and 2. A positive allele-specific or RBP perturbation result would strengthen causality but is not required for the initial targeted validation.

## Recommended disease, region, and cell type

| Gene | Disease to test | Brain region | Primary cell type | Prespecified comparator | Rationale |
|---|---|---|---|---|---|
| **SNCA** | **Lewy body dementia (LBD) and Parkinson's disease (PD)** | **Frontal association cortex, prioritizing BA9** | **Cortical excitatory neurons** | Inhibitory neurons; dopaminergic neurons only as a later PD extension | The same alternative-first-exon junction colocalizes with LBD in cortex and PD in frontal cortex BA9. Cortical neurons therefore test the tissue in which the human genetic association was observed. |
| **CTSH** | **Alzheimer's disease (AD)** | **Hippocampus** | **Microglia** | Astrocytes | CTSH has the strongest candidate-level colocalization signal, with hippocampus as the top tissue. Independent functional-genomic evidence also supports a role for CTSH in human microglial biology. |
| **PPP6R2** | **Schizophrenia (SCZ)** as the primary test; ALS as a secondary extension | **Dorsolateral prefrontal cortex/frontal cortex BA9** | **Cortical excitatory neurons** | Inhibitory neurons and oligodendrocytes | BA9 is the top sQTL tissue and is directly relevant to SCZ. PPP6R2 does not currently have strong evidence for specificity to one brain cell class, so the excitatory-neuron assignment is a testable starting hypothesis rather than an established fact. |

## 1. SNCA: cortical neuronal switch in LBD and PD

### Manuscript basis

- LBD locus: **rs7680557-A**, cortex, CLPP approximately 0.038.
- PD locus: **rs1471483-C**, frontal cortex BA9, CLPP approximately 0.025.
- Both loci increase the same junction: **chr4:89,835,692–89,836,127**.
- Resolved transcript pair: **ENST00000508895 / ENST00000618500**.
- The event changes the first exon and therefore represents a **5-prime regulatory or transcript-initiation switch**, not merely a change in total SNCA expression.

### Primary hypothesis

The SNCA alternative-first-exon switch is present in cortical excitatory neurons in frontal association cortex, and its transcript ratio differs by risk-allele dosage and/or synucleinopathy status.

### Proposed test

Use human frontal cortex BA9 from genotyped control, PD, and LBD donors when available. Enrich neuronal nuclei and distinguish excitatory from inhibitory neurons by fluorescence-activated nuclei sorting or single-nucleus profiling.

The minimum molecular assay should include:

- junction-specific RT-ddPCR or quantitative RT-PCR for the shared junction;
- a second assay distinguishing ENST00000508895 from ENST00000618500;
- targeted long-read amplicon sequencing to verify full-length transcript structure;
- total SNCA abundance as a separate measurement; and
- genotype at rs7680557 and rs1471483, or credible proxy variants when necessary.

The primary endpoint is the **ENST00000508895:ENST00000618500 transcript ratio** or equivalent junction-based switch fraction. Total SNCA expression should not substitute for this endpoint.

### Cell-type interpretation

SNCA is expressed in both excitatory and inhibitory cortical neurons. Excitatory neurons are the primary group because they provide a tractable cortical disease model and have documented vulnerability in Lewy body disease, but inhibitory neurons should be measured in the same experiment as a specificity control. Midbrain dopaminergic neurons are biologically important to PD, yet they should be treated as a later extension because the reported colocalization is anchored in cortex rather than substantia nigra.

### Public-data contribution

The public frontal-cortex **SnISOr-Seq dataset (GSE178175)** can be used first to determine whether the junction or either full-length transcript is detectable in neuronal nuclei. Because it contains only two adult frontal-cortex samples, it is a feasibility and cell-class-localization dataset, not an independent disease or genetic replication cohort.

## 2. CTSH: hippocampal microglial switch in AD

### Manuscript basis

- Disease: **AD**.
- Top molecular-QTL tissue: **hippocampus**.
- Lead signal: **rs12148472-T**.
- Maximum CLPP: approximately **0.386**, the strongest of these three candidates.
- Both eQTL and sQTL evidence are present, making it important to distinguish total CTSH abundance from transcript switching.

### Primary hypothesis

The AD-associated CTSH switch occurs primarily in hippocampal microglia and is separable from the locus's effect on total CTSH expression.

### Proposed test

The preferred experiment is targeted profiling of microglia isolated from genotyped human hippocampus. If fresh or frozen sorted microglia are impractical, use iPSC-derived human microglia for assay development and variant perturbation, followed by confirmation in hippocampal nuclei or RNA.

The minimum assay should include:

- targeted RT-ddPCR or quantitative RT-PCR for the resolved CTSH event;
- targeted long-read amplicon sequencing to confirm the transcript pair;
- total CTSH abundance measured independently;
- comparison of microglia with astrocytes; and
- genotype at rs12148472 or a credible proxy.

The primary endpoint is the **CTSH switch fraction within microglia**. A secondary analysis should ask whether genotype predicts switch fraction after adjustment for total CTSH abundance.

### Optional disease-relevant perturbation

Once the switch is detectable in microglia, an amyloid-beta uptake or exposure condition can test whether disease-relevant activation changes the transcript ratio. This should remain secondary to the genotype and cell-type analyses so that a stress response is not mistaken for the genetically anchored mechanism.

### Public-data contribution

The cell-type-resolved long-read atlas generated from **six adult human hippocampi** can establish whether the CTSH isoforms are present in broad glial versus neuronal populations and help refine assay design. It cannot replace a dedicated microglial validation if the relevant transcript has low coverage.

Independent functional work has linked the AD-protective CTSH locus to human microglial expression and showed that CTSH loss alters microglial transcription and amyloid-beta phagocytosis. This makes microglia the strongest prespecified cell type, with astrocytes retained as a biologically plausible comparator.

## 3. PPP6R2: frontal-cortical neuronal switch in SCZ

### Manuscript basis

- Colocalized diseases: **SCZ and ALS**.
- Top molecular-QTL tissue: **frontal cortex BA9**.
- Lead signal: **rs76300267-G**.
- Maximum CLPP: approximately **0.065**.
- Three events resolve to switch pairs, indicating a potentially complex splicing locus.
- The locus is classified as **splicing-led** rather than primarily abundance-led.

### Primary hypothesis

PPP6R2 transcript usage in frontal cortex BA9 differs in cortical excitatory neurons as a function of the SCZ-associated haplotype, while total PPP6R2 abundance changes little or does not explain the switch.

### Why SCZ is the first disease test

SCZ is the cleanest initial disease context because the strongest tissue assignment is frontal cortex BA9. ALS remains important evidence that the locus may act across disorders, but an ALS-focused first experiment would require a different anatomical and cellular hypothesis, such as motor cortex corticospinal neurons or spinal motor neurons. That should be a second study after the BA9 event is validated.

### Proposed test

Use genotyped human BA9 tissue with excitatory-neuron, inhibitory-neuron, and oligodendrocyte nuclei resolved by sorting or single-nucleus profiling. In parallel, iPSC-derived cortical excitatory neurons can support assay development and eventual isogenic editing of rs76300267 or the credible variant set.

The minimum assay should include:

- assays for **each of the three resolved PPP6R2 switch events**, rather than selecting only the strongest event after observing the data;
- targeted long-read sequencing to determine whether the events occur in the same or different full-length transcripts;
- separate measurement of total PPP6R2 abundance;
- comparison among excitatory neurons, inhibitory neurons, and oligodendrocytes; and
- genotype-stratified switch fractions.

The primary endpoint is the switch fraction for the prespecified event with the strongest original colocalization evidence. The other two resolved events are secondary but should be retained to reveal coordinated transcript remodeling.

### Cell-type caveat

PPP6R2 is not currently established as specific to a single brain cell class. The proposal should therefore say that it will **test whether the BA9 sQTL is carried by cortical excitatory neurons**, not that excitatory neurons are already known to be the causal cell type. The inhibitory-neuron and oligodendrocyte comparisons are essential for interpreting a negative result.

### Public-data contribution

GSE178175 can be screened for PPP6R2 full-length transcripts in frontal-cortex neuronal and non-neuronal nuclei. A positive observation would guide primer design and cell-type prioritization; absence would be inconclusive because of the dataset's small sample size and incomplete transcript capture.

## Shared experimental design

### Stage 0: public-data feasibility screen

1. Query GSE178175 for the exact SNCA and PPP6R2 junctions and transcript structures in frontal-cortex nuclei.
2. Query the adult human hippocampal long-read atlas for the CTSH transcript pair and broad cell-class distribution.
3. Lock the transcript definitions and primer locations before examining disease or genotype effects.

Public data can resolve assay feasibility and sometimes cell-class localization. It should not be presented as genotype or disease validation unless the dataset has the required donors, phenotypes, and coverage.

### Stage 1: molecular-collaborator validation

For every gene:

- confirm the exact event with two orthogonal measurements, preferably junction-specific RT-ddPCR and targeted long-read amplicon sequencing;
- measure switch fraction and total gene abundance separately;
- use multiple biological donors or independent iPSC lines, not only technical replicates;
- blind genotype or disease labels during initial assay quantification where feasible; and
- predefine the primary transcript ratio and direction of effect from the manuscript.

### Stage 2: allelic and mechanistic follow-up

Only events that pass Stage 1 should advance to:

- allele-specific expression or transcript-usage analysis in heterozygous donors;
- isogenic editing of the lead or credible causal variant;
- minigene testing if the regulatory interval is compact and interpretable; and
- perturbation of candidate RBPs followed by rescue.

RBP candidates should be treated as hypotheses until the transcript switch has been reproduced in the selected cell type. The manuscript's motif enrichments alone do not establish RBP causality.

## Go/no-go criteria

An event advances when all of the following are met:

1. The exact predicted junction or full-length transcript pair is detected.
2. The switch is reproducible across biological samples or independent differentiations.
3. The measured direction agrees with the original QTL prediction when genotype information is available.
4. The switch fraction adds information beyond total gene abundance.
5. At least one disease-relevant cell type has sufficient expression and assay dynamic range for perturbation.

Failure to detect an event in a small public long-read dataset is not a no-go result. Failure in adequately powered targeted assays across the primary and comparator cell types is substantially stronger evidence against pursuing that candidate.

## Practical priority

1. **CTSH — AD, hippocampal microglia:** highest probability of a clear functional result because the colocalization is strongest and independent microglial evidence already exists.
2. **SNCA — LBD/PD, BA9 cortical excitatory neurons:** strongest disease narrative and a shared switch across two synucleinopathies; particularly valuable for the paper if the alternative-first-exon event is confirmed.
3. **PPP6R2 — SCZ, BA9 cortical excitatory neurons:** high novelty but higher cell-type and transcript-structure uncertainty; best treated as a discovery-focused third candidate.

If resources permit only one molecular experiment before submission, prioritize a targeted **SNCA or CTSH** validation depending on whether the desired emphasis is disease narrative or probability of technical success. PPP6R2 is well suited to public-data localization and assay development while the first two are being validated.

## Sources

- IsoGraph manuscript [per-gene colocalization table](../../manuscripts/isograph-brain-manuscript/content/supplementary_tables/tableS9_per_gene_colocalization.tsv) and [Results](../../manuscripts/isograph-brain-manuscript/content/03.results.md).
- [Single-nucleus long-read sequencing of adult human frontal cortex (SnISOr-Seq; GSE178175)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175).
- [Single-cell long-read isoform atlas of the adult human hippocampus](https://www.nature.com/articles/s41593-024-01616-4).
- [Cell-type-specific SNCA transcription during Lewy body formation](https://pmc.ncbi.nlm.nih.gov/articles/PMC10666428/).
- [Cell-type-specific SNCA expression and vulnerability in Lewy body diseases](https://pmc.ncbi.nlm.nih.gov/articles/PMC8357687/).
- [Functional genomic evidence for the CTSH locus in AD and human microglia](https://pmc.ncbi.nlm.nih.gov/articles/PMC10516988/).
- [Human Protein Atlas: PPP6R2 single-cell expression summary](https://www.proteinatlas.org/ENSG00000100239-PPP6R2/single%2Bcell).
