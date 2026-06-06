# Giant-Module Collapse: Diagnosis and Options

## Problem

All genes collapse into one giant module (M000) even after increasing `latent_dim` to 32
and adding `alpha_switch_grid` calibration. Root causes:

1. **Abundance-abundance edges are universally dense**: All 27,746 edges exceed even the
   max grid threshold (0.90). They capture global expression variation (cell-type
   composition, sequencing depth), not gene-specific regulation. Any path through these
   bridges every module together.

2. **Connected-components is fragile to bridges**: Even a few hub genes correlated with
   genes in two otherwise distinct communities merge them into one component. Leiden/Louvain
   community detection handles this correctly; plain connected-components does not.

---

## Options

### Option A — Disable abundance-abundance edges (`allow_abundance_abundance=False`)

**What it does**: Excludes abundance↔abundance edges from the graph entirely. Switch-switch
and cross-channel (switch↔abundance) edges are kept.

**Why it helps**: The abundance-abundance edges are all > 0.90 correlation — they reflect
global expression (not gene-specific) and are the primary bridge collapsing all communities.

**Risk**: Loses abundance co-regulation signal, but since the current signal is degenerate
(all above any reasonable threshold), this is likely a net gain for community specificity.

**Parameters**: `allow_abundance_abundance=False`, `alpha_switch=0.6` (fixed, no grid sweep),
`alpha_abundance_grid` removed — with abundance-abundance edges disabled, the grid calibration
has no stopping criterion and wastes 7 graph-projection calls. Cross-channel edges fall back
to the main `alpha` threshold (default 0.70).

**Status**: Implemented 2026-05-13. Giant module persists with Option A alone; Options B+C added.

---

### Option B — Residual similarity before thresholding

**What it does**: After VAE reconstruction, subtract the mean reconstructed feature vector
across genes (removes the dominant global expression axis) before computing gene-gene
Pearson correlation. This is analogous to removing batch effects or the first PC.

**Implementation**: 2–3 lines of numpy in `vae.py`, before calling `_gene_similarity`.
Alternatively, regress out the top-k PCs of the reconstruction matrix.

**Why it helps**: The VAE latent space captures global expression trends in the first few
dimensions. When reconstructing 37k features, genes that share global co-variation appear
spuriously correlated even if their gene-specific patterns differ. Centering by gene removes
this shared baseline.

**Risk**: Could remove biologically meaningful co-regulation if it aligns with global trends.
Centering by gene (subtract each gene's mean) is safer than removing PCs.

**Status**: Implemented 2026-05-13 — double-centering added at `vae.py:265` (`X = X - X.mean(axis=0, keepdims=True)`).

---

### Option C — Community detection: Leiden/Louvain

**What it does**: Replace `nx.connected_components` in IsoGraph's module assignment with
Leiden community detection (via `leidenalg`). Leiden maximizes modularity within a graph,
finding communities even in dense, well-connected graphs.

**Why it helps**: Connected components merges any two communities linked by even a single
bridge gene. Leiden correctly separates functionally distinct communities that share a few
cross-connections. The resolution parameter γ controls granularity: higher γ → more,
smaller modules.

**Implementation**: Changes needed in IsoGraph's graph-to-module logic
(`isograph/models/multiplex.py` or wherever `nx.connected_components` is called for
module assignment). Requires `leidenalg` package (available via conda).

**Risk**: Changes the fundamental module-detection algorithm; results will differ from
all prior runs. Needs benchmarking on synthetic data to validate.

**Status**: Implemented 2026-05-13. `leiden_resolution: float | None` added to `VaeModelConfig`; set to `1.0` in all real-data configs. `_module_table` updated in all 4 model files with igraph fallback to connected_components when `leiden_resolution=None`.

---

## Recommendation

1. **Test Option A first** (cheap, already parameterized): if modules become interpretable,
   abundance-abundance edges were the problem and the fix is permanent.
2. **Apply Option B** if A alone doesn't fully resolve collapse (complementary, not exclusive).
3. **Plan Option C** as the correct long-term solution regardless of A/B outcome — it makes
   the pipeline robust to any future edge-density issues.
