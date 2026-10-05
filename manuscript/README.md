# Manuscript

Analysis-side figure/table builders and the outputs they read.

The manuscript text, final display numbering, legends and submission assembly live
in the [companion manuscript repository](https://github.com/heart-gen/isograph-brain-manuscript).
This directory includes historical and intermediate displays as well as submission outputs.
Numbers are never hand-edited here: each item regenerates from committed ledgers.

- `_h/` — figure and table builders (one per display item).
- `_m/figures/` — all real-data figure panels (PDF + PNG).
- `_m/main_tables/`, `_m/supp_tables/` — CSV + Markdown tables.

## Rebuilding

Figures use the `rnaseq` R environment; tables use the `isograph` python environment.

```bash
Rscript manuscript/_h/qtl_specificity_figure.R
python manuscript/_h/assemble_main_tables.py
python manuscript/_h/assemble_supp_tables.py
```

`_m/manuscript_supplement/` holds the supplementary table and Data S1–S28 exports.
Use the companion manuscript for their final captions and ordering.
