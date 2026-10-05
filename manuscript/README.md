# Manuscript

Every display item in the paper, its builder, and the analysis outputs it reads.
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
