# phASER feasibility screen — caudate

Can the existing BrainSEQ phASER release support an allele-specific test of an isoform switch? A read is only isoform-informative if it carries a heterozygous site **and** sequence unique to one isoform of the pair. `gene_ae` cannot answer this: its ENST labels name genomic intervals, and reads in shared exons are counted into several transcript features at once.

## Pre-registered thresholds

- at least **30** donors
- each with at least **5** informative reads
- informative on **both** isoforms of the pair (the statistic is a contrast)

## Result

- switch pairs with isoform-unique exonic sequence: **151**
- pairs with any informative site: **37**
- **pairs passing the screen: 3** (2 genes)
- median donors per pair: 93; best pair: 367

## Why pairs fail — and it is not what the QC tables suggested

- pairs failing **only** the both-isoforms rule, with donors and depth to spare: **25**
- pairs failing **only** on donor count: 2
- median depth among those donor-rich failures: **28 reads**

The phASER QC tables report a median `gene_ae` depth of 1-2 reads, which reads like a depth problem. It is not. That figure is the whole-transcript-interval aggregation; once attention is restricted to isoform-DISCRIMINATING exonic sequence, the donors that have any signal have plenty of it. **The binding constraint is structural: for most switch pairs only ONE of the two isoforms carries unique exonic sequence containing a heterozygous site.** A pair where one isoform's sequence is a subset of the other's — a truncation, an alternative last exon — has no unique sequence on that side at all, so there is nothing to contrast against however deeply it is sequenced.

**3 pairs clear the screen** and could proceed to allele orientation. That is too few to carry a claim on its own.

### The rescue this points at

Exonic-unique segments are a ONE-SIDED discriminator. A discriminating **junction** is two-sided by construction: isoform A splices one junction where isoform B splices another, so a junction-spanning read is informative for whichever isoform it came from. Counting reads that simultaneously cross a discriminating junction and carry a phased heterozygous site would test many of the 25 pairs that fail here for want of a second side. That needs allele-aware counting from the BAMs; `junction_tables/*_SJ.out.tab` are STAR junction counts and are **not** allele-aware, so they cannot substitute.

If stage 2 is attempted, note that orientation — not the model — is the hard part: `gw_phased = 1` makes haplotype A consistent WITHIN a donor but not the same biological allele ACROSS donors, so each donor must be oriented to the risk allele before pooling or a real cis effect averages toward zero.

## Top pairs by informative donors

| gene | T1 | T2 | donors (>= min reads) | sites | isoforms informative | testable |
|---|---|---|---|---|---|---|
| ENSG00000126214 | ENST00000553325.5 | ENST00000348520.10 | 367 | 41 | 1 | no |
| ENSG00000198053 | ENST00000356025.7 | ENST00000358771.5 | 279 | 4 | 1 | no |
| ENSG00000158604 | ENST00000477639.5 | ENST00000444131.2 | 251 | 3 | 1 | no |
| ENSG00000090975 | ENST00000280562.9 | ENST00000436074.6 | 250 | 5 | 1 | no |
| ENSG00000122218 | ENST00000241704.8 | ENST00000650154.1 | 249 | 4 | 1 | no |
| ENSG00000145335 | ENST00000508895.5 | ENST00000674129.1 | 243 | 2 | 1 | no |
| ENSG00000126214 | ENST00000553286.5 | ENST00000553325.5 | 231 | 4 | 1 | no |
| ENSG00000115239 | ENST00000406687.5 | ENST00000482829.5 | 228 | 3 | 2 | **yes** |
| ENSG00000092140 | ENST00000553504.5 | ENST00000549159.1 | 227 | 10 | 2 | **yes** |
| ENSG00000163512 | ENST00000479665.6 | ENST00000420543.6 | 226 | 1 | 1 | no |
| ENSG00000163512 | ENST00000479665.6 | ENST00000334100.10 | 222 | 2 | 1 | no |
| ENSG00000126214 | ENST00000445352.8 | ENST00000553325.5 | 209 | 3 | 1 | no |
| ENSG00000145335 | ENST00000336904.7 | ENST00000506244.5 | 173 | 2 | 1 | no |
| ENSG00000092140 | ENST00000206595.11 | ENST00000549159.1 | 170 | 12 | 2 | **yes** |
| ENSG00000158604 | ENST00000457408.7 | ENST00000289577.10 | 154 | 1 | 1 | no |
| ENSG00000158604 | ENST00000289577.10 | ENST00000477639.5 | 136 | 3 | 1 | no |
| ENSG00000090975 | ENST00000320201.10 | ENST00000546049.5 | 118 | 1 | 1 | no |
| ENSG00000158604 | ENST00000457408.7 | ENST00000477639.5 | 118 | 2 | 1 | no |
| ENSG00000092140 | ENST00000648008.1 | ENST00000552488.1 | 93 | 5 | 1 | no |
| ENSG00000145335 | ENST00000394991.8 | ENST00000508895.5 | 80 | 3 | 1 | no |
| ENSG00000090975 | ENST00000320201.10 | ENST00000451868.2 | 75 | 2 | 1 | no |
| ENSG00000090975 | ENST00000546049.5 | ENST00000436074.6 | 60 | 1 | 1 | no |
| ENSG00000115239 | ENST00000263634.8 | ENST00000406687.5 | 55 | 1 | 1 | no |
| ENSG00000163512 | ENST00000420543.6 | ENST00000463512.5 | 48 | 1 | 1 | no |
| ENSG00000145335 | ENST00000508895.5 | ENST00000394986.5 | 48 | 1 | 1 | no |

