"""Per-gene mechanistic deep-dive for the colocalization-prioritized disease genes.

For each target gene this joins the genetic-anchoring and downstream-biology layers we have
already built into a single mechanistic vignette:

  1. genetic anchor      -- coloc_isoform_events_combined / coloc_direction_combined:
                            which trait, sQTL vs eQTL, tissue, CLPP, risk allele + rsID,
                            signed direction of the risk allele on the QTL.
  2. the switch          -- whether the colocalized junction maps onto an IsoGraph switch
                            pair (junction_in_switch_pair), the switch pair, GTEx concordance,
                            GO-invisibility, and independent BrainSeq replication.
  3. coding consequence  -- the resolved event's structural_consequence.
  4. regulatory logic    -- rbp_switch_calls / rbp_regulon: which RBP motifs are called
                            switched in the gene, and which of those are also enriched in the
                            gene's IsoGraph module (candidate splice regulators of the switch).
  5. constraint / clinic -- gnomAD LOEUF + missense o/e from the clinical-consequence layer.

Each gene is classified as splicing-led (an sQTL that colocalizes onto a concordant IsoGraph
switch pair -- the IsoGraph-unique, GO-invisible case), expression-led (eQTL gene-level only),
or splicing-unresolved (an sQTL that does not map onto the switch pair). Writes one markdown
vignette per gene plus a panel-wide parquet of one row per gene.

Two parts, because the layers above come from different stages. ``--part events`` needs only the
coloc layer and writes ``deep_dive_events`` into stage 05, where stage 06 reads it;
``--part panel`` needs the clinical (06) and RBP (07) layers and writes the vignettes, the panel
and the remaining tables into 08_integration.

Deterministic; pure joins over existing parquets (no heavy compute). Usage:
  python -m isograph_benchmark.real_data.gene_deep_dive --part events
  python -m isograph_benchmark.real_data.gene_deep_dive --part panel
  python -m isograph_benchmark.real_data.gene_deep_dive --part panel --genes SNCA,CTSH,PPP6R2
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from isograph_benchmark.paths import cohort_dir, ensure_dir, stage_out

# Originally hand-prioritized panel (genetic_anchoring_v2), kept for reference. The CLI
# defaults to EVERY colocalized gene so multi-locus / multi-trait cases are not dropped;
# main-figure framing is reserved for the resolved cross-disease headliners (SNCA, CTSH).
PANEL = [
    "SNCA", "TPP1", "SCFD1", "PGS1", "PPP6R2", "GGNBP2",
    "CTSH", "TPCN1", "MRPS10", "MYO18A", "PBX1", "RBFA", "MED15",
]
_Q_ENRICH = 0.05

# Layer 6 -- known isoform biology per resolved splicing-led gene (curated, with Manubot
# citekeys where a specific source is confirmed and clearly-marked placeholders otherwise).
# Kept as data so the vignettes and the literature supp table regenerate deterministically.
# `refs` are Manubot citekeys; `[citation needed: ...]` marks a real finding whose exact
# citekey still needs to be pinned before submission (never fabricate a DOI/PMID).
#
# REFRESHED 2026-09-20 onto the current 30-gene anchored set. Genes curated for the legacy
# set are kept -- curation is knowledge, not a result, and several are still in the 160-gene
# panel -- and each row now carries whether the gene is in the current anchored set.
#
# `kind` says what is known, and the three values are NOT a quality ranking:
#
#   isoform_documented  a disease-relevant isoform program is established for this gene
#   gene_documented     the gene and its disease association are established; which isoform
#                       the risk variant selects is not characterised
#   novel_candidate     no established disease-specific isoform biology was found
#
# A `novel_candidate` is a NOMINATION, not a null result and not a failed control. The layer
# exists because the method is built to surface GO-invisible, under-characterised switching;
# a gene with no isoform literature is the intended output, and its absence from the
# literature is a statement about the literature, not about the evidence here. Nothing
# downstream may treat these rows as negatives.
_LITERATURE: dict[str, dict] = {
    "SNCA": {
        "kind": "isoform_documented",
        "text": (
            "SNCA carries an extensively documented alternative-splicing program that is "
            "disease-relevant in synucleinopathy: at least four alternative 5'UTR first exons plus "
            "internal exon-3/exon-5 skipping generate transcripts differentially expressed across PD "
            "and dementia-with-Lewy-bodies brain regions, and the coding splice variants "
            "(SNCA-126/112/98) modulate alpha-synuclein aggregation kinetics. The event IsoGraph "
            "resolved here is a 5'-end (alternative first exon) choice, matching the regulatory arm "
            "of that program rather than a coding change. Not in the current anchored set: SNCA is no "
            "longer concordant at resolution 2.0."),
        "refs": ["@doi:10.3389/fgene.2019.00584", "@doi:10.3390/genes9020063"],
    },
    "DLG1": {
        "kind": "isoform_documented",
        "text": (
            "DLG1/SAP97 is a canonical alternatively-spliced synaptic scaffold: N-terminal alpha vs "
            "beta isoforms, an internal I3 insert and additional cassette exons tune its PDZ/GK "
            "synaptic function. A DLG1 splice variant is reported at reduced cortical levels in "
            "early-onset schizophrenia, and DLG1 sits in the 3q29 schizophrenia locus, so an sQTL "
            "that shifts DLG1 isoform choice is a mechanistically plausible splicing-led route to SCZ "
            "risk."),
        "refs": ["@doi:10.1038/tp.2015.154"],
    },
    "CTSH": {
        "kind": "isoform_documented",
        "text": (
            "CTSH has reported brain isoform usage differences relevant to psychiatric phenotypes. "
            "Not in the current anchored set."),
        "refs": ["@doi:10.1038/s41386-023-01542-2"],
    },
    "ARVCF": {
        "kind": "gene_documented",
        "text": (
            "ARVCF sits in the 22q11.2 locus (and carries a reported association with COMT) and is "
            "itself a modulator of pre-mRNA splicing -- it interacts with SRSF1, DDX5 and hnRNP H2 "
            "and alters alternative-splicing activity -- so a splicing-led ARVCF event is consistent "
            "with both its locus and its molecular role. Not in the current anchored set."),
        "refs": ["@doi:10.1038/sj.mp.4001586"],
    },
    "PRDM2": {
        "kind": "isoform_documented",
        "text": (
            "PRDM2/RIZ is the clearest documented isoform program in the anchored set: it is "
            "transcribed from alternative promoters into RIZ1, which carries the PR (SET-like) "
            "methyltransferase domain, and RIZ2, which lacks it, and the RIZ1:RIZ2 balance -- not "
            "total PRDM2 -- is what changes in disease states where the gene has been studied. The "
            "colocalizing ALS event here is a 5'-end choice, which is the same class of event as that "
            "switch, so the variant plausibly acts by selecting between PR-positive and PR-negative "
            "products rather than by changing gene dosage. Brain-specific and ALS-specific isoform "
            "work is not established."),
        "refs": ["[citation needed: RIZ1/RIZ2 alternative-promoter isoforms and the PR-domain balance]"],
    },
    "IFNAR2": {
        "kind": "isoform_documented",
        "text": (
            "IFNAR2 produces functionally distinct products by alternative splicing -- a full-length "
            "signalling receptor, a truncated form lacking most of the cytoplasmic domain, and a "
            "soluble form -- and the ratio between them sets the cell's type-I interferon response "
            "rather than its receptor abundance. IFNAR2 is an established Alzheimer's GWAS locus, and "
            "the anchored event is splicing-specific on the genetics (no eQTL instrument at all), so "
            "isoform choice is the mechanism this locus most plausibly acts through."),
        "refs": ["[citation needed: IFNAR2 transmembrane vs truncated vs soluble isoforms]"],
    },
    "DNAJA3": {
        "kind": "isoform_documented",
        "text": (
            "DNAJA3/Tid1 is a mitochondrial HSP40 co-chaperone with two long-standing splice forms "
            "that differ at the C-terminus and have been reported to act in opposite directions on "
            "apoptosis, so isoform choice rather than total level is the functional variable. A "
            "schizophrenia-specific isoform role is not established."),
        "refs": ["[citation needed: Tid1-L / Tid1-S splice forms and opposing apoptotic effects]"],
    },
    "TPP1": {
        "kind": "gene_documented",
        "text": (
            "TPP1/CLN2 encodes lysosomal tripeptidyl peptidase 1; its loss causes late-infantile "
            "neuronal ceroid lipofuscinosis, and it is one of the best-characterised lysosomal genes "
            "in neurodegeneration. What is not established is which TPP1 isoform an ALS risk variant "
            "selects -- the disease literature is about enzyme deficiency, not about isoform choice."),
        "refs": [],
    },
    "TMEM175": {
        "kind": "gene_documented",
        "text": (
            "TMEM175 is an established Parkinson's risk gene and a lysosomal potassium/proton channel "
            "whose coding variant M393T reduces channel function; the locus is one of the "
            "best-supported in PD. Isoform-level biology is not established, and in this panel "
            "TMEM175's own genetics are weak (sQTL PP4 ~ 0 against a large eCAVIAR CLPP), so it is "
            "carried as a documented gene with undocumented and unsupported isoform evidence."),
        "refs": [],
    },
    "CTSB": {
        "kind": "gene_documented",
        "text": (
            "Cathepsin B is a lysosomal cysteine protease and an established Parkinson's GWAS gene, "
            "functionally tied to the GBA/lysosomal axis and to alpha-synuclein degradation. Which "
            "CTSB isoform the risk variant selects is not characterised."),
        "refs": [],
    },
    "SNAP91": {
        "kind": "gene_documented",
        "text": (
            "SNAP91/AP180 assembles clathrin at the synapse and is an established schizophrenia GWAS "
            "gene with reported effects on synaptic development. Its paralogue PICALM has documented "
            "disease-relevant isoform biology in Alzheimer's, which makes an isoform-level mechanism "
            "plausible here, but SNAP91's own isoform program is not characterised."),
        "refs": [],
    },
    "NT5C2": {
        "kind": "gene_documented",
        "text": (
            "NT5C2 is a cytosolic 5'-nucleotidase, an established schizophrenia GWAS gene, and the "
            "cause of hereditary spastic paraplegia SPG45 when lost. Isoform choice in brain is not "
            "characterised."),
        "refs": [],
    },
    "SPG7": {
        "kind": "gene_documented",
        "text": (
            "SPG7/paraplegin is a subunit of the mitochondrial m-AAA protease; biallelic loss causes "
            "hereditary spastic paraplegia. The gene is well documented; its isoform biology is not."),
        "refs": [],
    },
    "VAMP2": {
        "kind": "gene_documented",
        "text": (
            "VAMP2/synaptobrevin-2 is a core SNARE of synaptic vesicle fusion and is strongly "
            "LoF-constrained. The gene is textbook; a disease-relevant isoform program is not "
            "established."),
        "refs": [],
    },
    "DOC2A": {
        "kind": "gene_documented",
        "text": (
            "DOC2A is a calcium sensor for spontaneous neurotransmitter release and sits in the "
            "16p11.2 locus. It is the most splicing-specific gene in this panel on the genetics -- "
            "strong sQTL colocalization in nine tissues with no eQTL instrument at all -- while its "
            "isoform biology is uncharacterised, which is exactly the gap this layer is meant to "
            "mark."),
        "refs": [],
    },
    "INO80E": {
        "kind": "gene_documented",
        "text": (
            "INO80E is a subunit of the INO80 chromatin-remodelling complex and lies in the 16p11.2 "
            "schizophrenia locus. Its genetics here colocalize on both modalities and are not "
            "splicing-specific; isoform biology is not established."),
        "refs": [],
    },
    "PTPRN": {
        "kind": "gene_documented",
        "text": (
            "PTPRN/IA-2 is a dense-core vesicle transmembrane protein and a well-known autoantigen in "
            "type 1 diabetes, with a documented role in neuroendocrine secretion. Its isoform usage "
            "in brain, and in ALS, is not characterised."),
        "refs": [],
    },
    "PIGQ": {
        "kind": "gene_documented",
        "text": (
            "PIGQ acts in the first step of GPI-anchor biosynthesis; biallelic variants cause a "
            "developmental and epileptic encephalopathy. GPI-anchoring is dosage-sensitive, but no "
            "isoform-level disease biology is established."),
        "refs": [],
    },
    "CTC1": {
        "kind": "gene_documented",
        "text": (
            "CTC1 is part of the CST telomere-maintenance complex; its loss causes Coats plus / "
            "cerebroretinal microangiopathy. No isoform-level disease biology is established."),
        "refs": [],
    },
    "FLCN": {
        "kind": "gene_documented",
        "text": (
            "FLCN/folliculin is the Birt-Hogg-Dube gene and a GTPase-activating protein in lysosomal "
            "amino-acid sensing upstream of mTORC1. The gene is well documented; its brain isoform "
            "biology is not, and its colocalization here is weak."),
        "refs": [],
    },
    "NADSYN1": {
        "kind": "gene_documented",
        "text": (
            "NADSYN1 completes NAD+ synthesis and sits in a locus shared with DHCR7, which "
            "complicates gene attribution at the signal. Isoform biology is not established."),
        "refs": [],
    },
    "B3GAT1": {
        "kind": "gene_documented",
        "text": (
            "B3GAT1/GlcAT-P makes the HNK-1 glycan carried by neural adhesion molecules and has "
            "documented roles in synaptic plasticity. Which isoform a schizophrenia risk variant "
            "selects is not known."),
        "refs": [],
    },
    "BAIAP3": {
        "kind": "gene_documented",
        "text": (
            "BAIAP3 is a Munc13-family protein controlling dense-core vesicle secretion, so it sits "
            "in the same secretory machinery as several other genes in this panel. Isoform biology is "
            "not established, and its usage range here touches the detection floor."),
        "refs": [],
    },
    "CRELD2": {
        "kind": "gene_documented",
        "text": (
            "CRELD2 is an ER-stress-responsive secreted protein induced through the ATF6 arm of the "
            "unfolded protein response. It has the broadest SMR support in the panel (48 of 51 "
            "probes) and no documented isoform program."),
        "refs": [],
    },
    "GSTO2": {
        "kind": "gene_documented",
        "text": (
            "GSTO2 is a glutathione S-transferase omega-family enzyme with reported associations to "
            "age-related phenotypes; its brain isoform biology is not established. The junction here "
            "does validate in the PSI catalogue, so the switch is measured even though the literature "
            "is silent."),
        "refs": [],
    },
    "TBC1D15": {
        "kind": "gene_documented",
        "text": (
            "TBC1D15 is a Rab7 GTPase-activating protein at the mitochondria-lysosome interface, a "
            "pathway central to Parkinson's mitophagy. Its PD-associated splicing-led switch has no "
            "established isoform literature -- a mechanistically suggestive nomination rather than a "
            "confirmation."),
        "refs": [],
    },
    "THAP3": {
        "kind": "novel_candidate",
        "text": (
            "THAP3 is a THAP-domain transcription factor, a family in which THAP1 is the DYT6 "
            "dystonia gene; THAP3 itself is little characterised in brain. It is one of the four "
            "splicing-specific genes here, which makes it a first-order nomination rather than a "
            "footnote."),
        "refs": [],
    },
    "PCGF3": {
        "kind": "novel_candidate",
        "text": (
            "PCGF3 is a Polycomb RING-finger subunit of a non-canonical PRC1 complex; no "
            "disease-specific isoform biology is established, and its usage range here touches the "
            "detection floor."),
        "refs": [],
    },
    "RPAIN": {
        "kind": "novel_candidate",
        "text": (
            "RPAIN/RIP is an RPA-interacting nuclear import factor with no established brain disease "
            "isoform biology. Its junction validates in the PSI catalogue and in the recount, so the "
            "switch is well measured and the literature is simply absent."),
        "refs": [],
    },
    "TARBP1": {
        "kind": "novel_candidate",
        "text": (
            "TARBP1 is a TRBP-related RNA methyltransferase with no established disease isoform "
            "biology; its usage range touches the detection floor here."),
        "refs": [],
    },
    "RBFA": {
        "kind": "novel_candidate",
        "text": (
            "RBFA is a mitoribosome assembly factor with essentially no brain disease literature; the "
            "switch is nominated here on genetics alone and is not measurable in the orthogonal "
            "assays."),
        "refs": [],
    },
    "GPR135": {
        "kind": "novel_candidate",
        "text": (
            "GPR135 is an orphan G-protein-coupled receptor with minimal functional characterisation, "
            "which makes any isoform claim premature; the colocalization here is weak."),
        "refs": [],
    },
    "PPP6R2": {
        "kind": "novel_candidate",
        "text": (
            "PPP6R2 (PP6 regulatory subunit) colocalizes as a splicing-led switch across both ALS and "
            "SCZ, the most cross-trait of the panel, but disease-specific isoform biology is not "
            "established -- a cross-trait nomination for follow-up. Both long-read and the junction "
            "recount confirm the switch, so the evidence is orthogonal even where the literature is "
            "not."),
        "refs": [],
    },
    "GGNBP2": {
        "kind": "novel_candidate",
        "text": (
            "GGNBP2/ZNF403 (17q12) is LoF-constrained and was nominated as a splicing-led switch in "
            "ALS on the legacy set; its isoform biology in neurodegeneration is uncharacterised. Not "
            "in the current anchored set."),
        "refs": [],
    },
    "RTEL1": {
        "kind": "novel_candidate",
        "text": (
            "RTEL1 (telomere-maintenance helicase; AD/SCZ locus) has documented alternative "
            "C-terminal isoforms in other tissues, but a brain disease-specific splice role is not "
            "established. Not in the current anchored set."),
        "refs": [],
    },
    "PGS1": {
        "kind": "novel_candidate",
        "text": (
            "PGS1 (mitochondrial phospholipid biosynthesis) was nominated in ALS with no established "
            "disease isoform biology. Not in the current anchored set."),
        "refs": [],
    },
    "CDIP1": {
        "kind": "novel_candidate",
        "text": (
            "CDIP1 (cell-death-inducing p53 target) was nominated as a two-event splicing-led switch "
            "in schizophrenia; its isoform biology is uncharacterised. Not in the current anchored "
            "set."),
        "refs": [],
    },
    "PRRC2B": {
        "kind": "novel_candidate",
        "text": (
            "PRRC2B is LoF-constrained and was nominated in schizophrenia with no established disease "
            "isoform literature. Not in the current anchored set."),
        "refs": [],
    },
    "TPCN1": {
        "kind": "novel_candidate",
        "text": (
            "TPCN1 (endolysosomal two-pore calcium channel; AD locus) was nominated with no "
            "established disease isoform biology. Not in the current anchored set."),
        "refs": [],
    },
}


_KIND_NOTE = {
    "isoform_documented": "A disease-relevant isoform program is established for this gene.",
    "gene_documented": ("The gene and its disease association are established; which isoform "
                        "the risk variant selects is not characterised."),
    "novel_candidate": ("No established disease-specific isoform biology was found. That is a "
                        "statement about the literature, not about the evidence here: an "
                        "under-characterised switch is what this method is built to surface, "
                        "and this row is a nomination rather than a null result."),
}


def _literature_lines(gene: str) -> list[str]:
    """Section 6 (literature) for a vignette; empty if the gene is not curated."""
    lit = _LITERATURE.get(gene)
    if lit is None:
        return []
    out = ["## 6. Literature (known isoform biology)", lit["text"], "",
           f"_Curation: {lit['kind']}._ {_KIND_NOTE[lit['kind']]}"]
    if lit["refs"]:
        out.append("")
        out.append("_References:_ " + "; ".join(lit["refs"]))
    out.append("")
    return out


def _all_coloc_genes() -> list[str]:
    ev = pd.read_parquet(stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet"))
    return sorted(ev["gene_name"].dropna().unique())


def _ens(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split(".").str[0]


def _load(part: str = "panel") -> dict:
    ev = pd.read_parquet(stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet"))
    di = pd.read_parquet(stage_out("anchoring.coloc", "coloc_direction_combined.parquet"))
    if part == "events":
        for df in (ev, di):
            if "gene" in df.columns:
                df["ens"] = _ens(df["gene"])
        return {"ev": ev, "di": di}
    calls = pd.read_parquet(stage_out("regulation", "rbp", "rbp_switch_calls.parquet"))
    regulon = pd.read_parquet(stage_out("regulation", "rbp", "rbp_regulon.parquet"))
    con = []
    for tree in ("brainseq", "gtex"):
        for f in cohort_dir(tree).glob(
                "*/_m/isograph_vae/clinical_consequence/gene_constraint.parquet"):
            con.append(pd.read_parquet(f))
    constraint = pd.concat(con, ignore_index=True) if con else pd.DataFrame(
        columns=["gene", "loeuf", "mis_oe", "region"])
    for df in (ev, di, calls, constraint):
        if "gene" in df.columns:
            df["ens"] = _ens(df["gene"])
    return {"ev": ev, "di": di, "calls": calls, "regulon": regulon, "constraint": constraint}


def _rbp_layer(gene_calls: pd.DataFrame, regulon: pd.DataFrame) -> dict:
    """Switched RBP motifs for the gene, flagged by whether they are also module-enriched."""
    switched = gene_calls[gene_calls["switched"]]
    if switched.empty:
        return {"switched_rbps": [], "module_enriched_rbps": []}
    # recurrence of each switched RBP across the gene's region/module memberships
    rec = switched.groupby("rbp")["region"].nunique().sort_values(ascending=False)
    # module-level enrichment: match the gene's (region, module_id) to the regulon
    mods = switched[["region", "module_id"]].drop_duplicates()
    reg = regulon.merge(mods, on=["region", "module_id"], how="inner")
    reg = reg[reg["q"] < _Q_ENRICH]
    # keep only RBPs that are BOTH switched in the gene AND enriched in its module
    reg = reg[reg["rbp"].isin(set(switched["rbp"]))].sort_values("q")
    enr = (reg.groupby("rbp")["q"].min().sort_values().index.tolist())
    return {
        "switched_rbps": [f"{r}({n})" for r, n in rec.items()],
        "module_enriched_rbps": enr,
    }


def _classify(g_ev: pd.DataFrame) -> tuple[str, bool]:
    """Return (verdict, go_invisible_resolved)."""
    sqtl = g_ev[g_ev["kind"] == "sQTL"]
    resolved = sqtl[(sqtl["junction_in_switch_pair"]) & (sqtl["concordant"] == True)]  # noqa: E712
    if not resolved.empty:
        go_inv = bool(resolved["go_invisible"].any())
        return "splicing-led (IsoGraph-resolved switch)", go_inv
    if (g_ev["kind"] == "eQTL").all():
        return "expression-led (eQTL gene-level)", bool(g_ev["go_invisible"].any())
    if not sqtl.empty:
        return "splicing (sQTL not resolved to switch pair)", bool(g_ev["go_invisible"].any())
    return "unclassified", False


def _gene_row(gene: str, d: dict) -> dict | None:
    ev = d["ev"][d["ev"]["gene_name"] == gene]
    if ev.empty:
        return None
    ens = ev["ens"].iloc[0]
    di = d["di"][d["di"]["ens"] == ens]
    calls = d["calls"][d["calls"]["ens"] == ens]
    con = d["constraint"][d["constraint"]["ens"] == ens]
    rbp = _rbp_layer(calls, d["regulon"])
    verdict, go_inv = _classify(ev)
    resolved = ev[(ev["kind"] == "sQTL") & (ev["junction_in_switch_pair"]) &
                  (ev["concordant"] == True)]  # noqa: E712
    top = ev.sort_values("clpp", ascending=False).iloc[0]
    return {
        "gene": gene,
        "ens": ens,
        "traits": ",".join(sorted(ev["trait"].str.upper().unique())),
        "kinds": ",".join(sorted(ev["kind"].unique())),
        "n_events": int(len(ev)),
        "max_clpp": float(ev["clpp"].max()),
        "top_tissue": str(top["tissue"]),
        "top_rsid": str(top["best_rsid"]),
        "top_risk_allele": str(top["risk_allele"]),
        "resolved_to_switch_pair": bool(not resolved.empty),
        "n_resolved_events": int(len(resolved)),
        "multi_locus": bool(len(resolved) >= 2),
        "concordant_traits": ",".join(sorted(resolved["trait"].str.upper().unique()))
                             if not resolved.empty else "",
        "brainseq_replicates": bool(ev["brainseq_replicates_switch"].eq(True).any()),
        "go_invisible": bool(go_inv),
        "loeuf": float(con["loeuf"].min()) if not con.empty else np.nan,
        "mis_oe": float(con["mis_oe"].min()) if not con.empty else np.nan,
        "n_switched_rbps": len(rbp["switched_rbps"]),
        "module_enriched_rbps": ",".join(rbp["module_enriched_rbps"][:10]),
        "verdict": verdict,
        "_ev": ev, "_di": di, "_rbp": rbp, "_con": con, "_resolved": resolved,
    }


def _vignette(row: dict) -> str:
    ev, di, rbp = row["_ev"], row["_di"], row["_rbp"]
    L = [f"# {row['gene']} — mechanistic deep-dive", ""]
    L.append(f"**Verdict:** {row['verdict']}"
             f"{' · GO-invisible' if row['go_invisible'] else ''}"
             f"{' · replicates in BrainSeq' if row['brainseq_replicates'] else ''}")
    L.append("")
    L.append(f"- **Traits:** {row['traits']}  ·  **QTL kinds:** {row['kinds']}  ·  "
             f"**max CLPP:** {row['max_clpp']:.2f} ({row['top_tissue']})")
    loeuf = "n/a" if np.isnan(row["loeuf"]) else f"{row['loeuf']:.3f}"
    L.append(f"- **Constraint:** LOEUF {loeuf}  ·  missense o/e "
             f"{'n/a' if np.isnan(row['mis_oe']) else f'{row['mis_oe']:.2f}'}")
    L.append("")
    L.append("## 1–3. Genetic anchor → switch → coding consequence")
    L.append("| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |")
    L.append("|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|")
    for r in ev.sort_values(["kind", "clpp"], ascending=[True, False]).itertuples():
        sc = (str(r.structural_consequence) or "").strip() or "—"
        L.append(f"| {r.trait} | {r.kind} | {r.tissue} | {r.clpp:.2f} | "
                 f"{'yes' if r.junction_in_switch_pair else 'no'} | "
                 f"{'yes' if r.concordant == True else ('no' if r.concordant == False else '—')} | "  # noqa: E712
                 f"{'yes' if r.go_invisible else 'no'} | {sc[:48]} | {str(r.resolved_event)[:90]} |")
    if not di.empty:
        L.append("")
        L.append("### Signed risk-allele direction (colocalized loci with allele matching)")
        L.append("| trait | tissue | rsID | risk allele | risk QTL effect | direction |")
        L.append("|-------|--------|------|-------------|-----------------|-----------|")
        for r in di.itertuples():
            L.append(f"| {r.trait} | {r.tissue} | {r.best_rsid} | {r.risk_allele} | "
                     f"{r.risk_qtl_effect} | {r.direction} |")
    L.append("")
    L.append("## 4. Regulatory logic (RBP motifs in switched exons)")
    if rbp["module_enriched_rbps"]:
        L.append(f"- **Switched *and* module-enriched (q<{_Q_ENRICH}) RBPs:** "
                 f"{', '.join(rbp['module_enriched_rbps'][:15])}")
    L.append(f"- **All switched-motif RBPs (recurrence across regions):** "
             f"{', '.join(rbp['switched_rbps'][:20]) or '—'}")
    L.append("")
    L.append("## 5. Interpretation")
    L.append(_interpretation(row))
    L.append("")
    L.extend(_literature_lines(row["gene"]))
    return "\n".join(L)


def _interpretation(row: dict) -> str:
    if row["resolved_to_switch_pair"]:
        s = (f"A splicing QTL colocalizes onto an IsoGraph switch pair for {row['gene']} "
             f"in {row['concordant_traits']}")
        if "," in row["concordant_traits"]:
            s += " — the same switch is genetically anchored across more than one trait"
        s += (". This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts "
              "through isoform choice, not gene dosage")
        if row["go_invisible"]:
            s += ", in a GO-invisible module a pathway-enrichment scan would miss"
        if row["brainseq_replicates"]:
            s += "; the switch also replicates in an independent BrainSeq cohort"
        return s + "."
    if "expression-led" in row["verdict"]:
        return (f"{row['gene']} colocalizes as an eQTL (gene-level expression), with no splicing "
                "event resolving to an IsoGraph switch pair — an honest expression-confounded case "
                "that abundance networks would also capture.")
    return (f"{row['gene']} has a colocalizing sQTL, but the junction does not map onto the "
            "IsoGraph switch pair for the tissue — splicing-associated but not resolved to a "
            "switch; a candidate for deeper transcript-level follow-up.")


def run(genes: list[str], part: str = "panel") -> pd.DataFrame:
    d = _load(part)
    if part == "events":
        # Stage 05: the per-event table needs only the coloc layer, and stage 06 reads it, so it
        # is written before (and apart from) the panel, which needs stages 06 and 07.
        out_dir = ensure_dir(stage_out("anchoring", "deep_dive"))
        ens_set = set()
        for g in genes:
            g_ev = d["ev"][d["ev"]["gene_name"] == g]
            if not g_ev.empty:
                ens_set.add(g_ev["ens"].iloc[0])
        _write_events(d, ens_set, out_dir)
        print(f"deep-dive events over {len(ens_set)} genes -> {out_dir}")
        return d["ev"][d["ev"]["ens"].isin(ens_set)]
    out_dir = ensure_dir(stage_out("integration", "deep_dive"))
    rows = []
    for g in genes:
        row = _gene_row(g, d)
        if row is None:
            print(f"  {g}: no colocalized isoform events; skipping.")
            continue
        (out_dir / f"{g}.md").write_text(_vignette(row))
        rows.append({k: v for k, v in row.items() if not k.startswith("_")})
    panel = pd.DataFrame(rows)
    # order: resolved splicing-led first, multi-locus above single, then by CLPP
    panel = panel.sort_values(["resolved_to_switch_pair", "n_resolved_events", "max_clpp"],
                              ascending=[False, False, False]).reset_index(drop=True)
    panel.to_parquet(out_dir / "deep_dive_panel.parquet", index=False)
    _write_panel_md(panel, out_dir)
    _write_supp_tables(d, set(panel["ens"]), out_dir, panel)
    print(f"deep-dive over {len(panel)} genes -> {out_dir}")
    return panel


def _emit(df: pd.DataFrame, out_dir, name: str) -> None:
    """Write a supplementary table as both parquet (analysis) and tsv (reader-facing)."""
    df.to_parquet(out_dir / f"{name}.parquet", index=False)
    df.to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)
    print(f"  supp table {name}: {len(df)} rows")


def _write_events(d: dict, ens_set: set, out_dir) -> None:
    """(1) events - one row per colocalized isoform event (anchor -> switch -> consequence),
    merged with the signed risk-allele direction where allele matching succeeded. Coloc layer
    only; written by ``--part events`` into stage 05, where stage 06 reads it.
    """
    ev = d["ev"][d["ev"]["ens"].isin(ens_set)].copy()
    # one direction row per locus (the table is keyed per-junction; dedup avoids fan-out)
    dkeep = (d["di"][["ens", "trait", "tissue", "best_rsid", "variant_id", "ref", "alt",
                      "risk_beta", "slope", "direction"]]
             .drop_duplicates(["ens", "trait", "tissue", "best_rsid"]))
    events = ev.merge(dkeep, on=["ens", "trait", "tissue", "best_rsid"], how="left")
    cols = ["gene_name", "ens", "trait", "case", "kind", "tissue", "best_rsid", "risk_allele",
            "risk_qtl_effect", "direction", "variant_id", "ref", "alt", "junction", "clpp",
            "go_invisible", "switch_pair", "junction_in_switch_pair", "structural_consequence",
            "concordant", "brainseq_region", "brainseq_replicates_switch", "resolved_event"]
    events = events[[c for c in cols if c in events.columns]].sort_values(
        ["gene_name", "trait", "clpp"], ascending=[True, True, False])
    _emit(events, out_dir, "deep_dive_events")


def _write_supp_tables(d: dict, ens_set: set, out_dir,
                       panel: pd.DataFrame | None = None) -> None:
    """Machine-readable per-gene tables so readers can reconstruct any gene's deep-dive.

    (2) rbp        - per gene, RBP motifs both switched in the gene and enriched in its module;
    (3) exons      - per gene/region/exon: switched vs constitutive, CDS overlap, ClinVar P/LP;
    (4) literature - curated isoform biology for the resolved splicing-led genes.
    The per-event table (1) is ``_write_events``, written by ``--part events`` into stage 05.
    """
    # (2) per-gene switched + module-enriched RBP regulators
    ens2sym = dict(zip(d["ev"]["ens"], d["ev"]["gene_name"]))
    calls = d["calls"][(d["calls"]["ens"].isin(ens_set)) & (d["calls"]["switched"])]
    reg = d["regulon"][d["regulon"]["q"] < _Q_ENRICH][
        ["region", "module_id", "rbp", "enrichment", "q"]]
    rbp = calls.merge(reg, on=["region", "module_id", "rbp"], how="inner")
    rbp["gene_name"] = rbp["ens"].map(ens2sym)
    rbp = rbp[["gene_name", "ens", "region", "module_id", "go_invisible", "rbp",
               "enrichment", "q"]].sort_values(["gene_name", "q"])
    _emit(rbp, out_dir, "deep_dive_rbp")

    # (3) per-gene/exon clinical annotation (SNCA-style read for every gene)
    frames = []
    for tree in ("brainseq", "gtex"):
        for f in cohort_dir(tree).glob(
                "*/_m/isograph_vae/clinical_consequence/exon_clinvar.parquet"):
            x = pd.read_parquet(f)
            x["ens"] = _ens(x["gene"])
            x = x[x["ens"].isin(ens_set)]
            if not x.empty:
                x["region"] = f.parents[2].name if "region" not in x.columns else x["region"]
                frames.append(x)
    if frames:
        exons = pd.concat(frames, ignore_index=True)
        exons["gene_name"] = exons["ens"].map(ens2sym)
        keep = ["gene_name", "ens", "region", "chrom", "start", "end", "length", "switched",
                "cds_overlap", "n_plp", "n_clinvar", "go_invisible"]
        exons = exons[[c for c in keep if c in exons.columns]].sort_values(
            ["gene_name", "region", "start"])
        _emit(exons, out_dir, "deep_dive_exon_clinical")

    # (4) curated literature layer (known isoform biology).
    # `curation` is what is KNOWN about the gene, never how good the evidence here is: a
    # novel_candidate is a nomination the literature has not reached, which is the output
    # this method exists to produce. `in_anchored_set` says whether the gene is currently
    # splicing-led, so a reader can see the layer's coverage without re-deriving it.
    sym_set = set(ens2sym.values())
    anchored = (set(panel.loc[panel["verdict"].str.startswith("splicing-led"), "gene"])
                if panel is not None and "verdict" in panel.columns else set())
    lit_rows = [
        {"gene_name": g,
         "literature": v["text"],
         "references": "; ".join(v["refs"]) if v["refs"] else "",
         "curation": v["kind"],
         "in_anchored_set": g in anchored}
        for g, v in _LITERATURE.items() if g in sym_set
    ]
    if lit_rows:
        lit = pd.DataFrame(lit_rows).sort_values("gene_name").reset_index(drop=True)
        _emit(lit, out_dir, "deep_dive_literature")


def _write_panel_md(panel: pd.DataFrame, out_dir) -> None:
    L = ["# Per-gene mechanistic deep-dive — panel summary", "",
         "Every colocalized disease gene, one row each (ranked splicing-led first, multi-locus "
         "above single). `resolved events` = colocalizing sQTLs that map onto a concordant "
         "IsoGraph switch pair (the splicing-led, DTU-without-DGE class); `multi-locus` flags "
         "genes with >=2 such events (multiple significant colocalizations, incl. cross-trait). "
         "Main-figure framing is a PI decision taken from the anchored gene summary "
         "(`08_integration/_m/anchored_gene_summary/`), not from this table; as of 2026-09-20 "
         "Fig 4A provisionally draws PRDM2. The remainder are supporting vignettes. See "
         "`<GENE>.md` for each.", "",
         "| gene | traits | kinds | max CLPP | LOEUF | resolved events | multi-locus | concordant traits | BrainSeq rep | GO-inv | verdict |",
         "|------|--------|-------|----------|-------|-----------------|-------------|-------------------|--------------|--------|---------|"]
    for r in panel.itertuples():
        loeuf = "n/a" if pd.isna(r.loeuf) else f"{r.loeuf:.2f}"
        L.append(f"| {r.gene} | {r.traits} | {r.kinds} | {r.max_clpp:.2f} | {loeuf} | "
                 f"{r.n_resolved_events} | {'yes' if r.multi_locus else '—'} | "
                 f"{r.concordant_traits or '—'} | "
                 f"{'yes' if r.brainseq_replicates else 'no'} | "
                 f"{'yes' if r.go_invisible else 'no'} | {r.verdict} |")
    (out_dir / "DEEP_DIVE_PANEL.md").write_text("\n".join(L) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Per-gene mechanistic deep-dive of coloc genes.")
    ap.add_argument("--genes", default="",
                    help="comma-separated gene symbols (default: every colocalized gene).")
    ap.add_argument("--part", choices=("events", "panel"), required=True,
                    help="events: per-event table from the coloc layer (stage 05; stage 06 reads "
                         "it). panel: vignettes, panel, RBP/exon/literature tables (stage 08; "
                         "needs the 06 clinical and 07 RBP layers).")
    args = ap.parse_args()
    genes = [g.strip() for g in args.genes.split(",") if g.strip()] or _all_coloc_genes()
    run(genes, args.part)


if __name__ == "__main__":
    main()
