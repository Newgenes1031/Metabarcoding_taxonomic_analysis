# metabarcoding_taxonomy

A Python package to automate ASV taxonomy filtering, statistics calculation and QIIME2-style cumulative barplots for metabarcoding data.

## Features

- **Filter & truncate** ASV taxonomy columns by a variety of criteria (missing ranks, `incertae`, `_sp`, numeric IDs, `uncultured`/`unidentified`/`candidum`/`Candidatus`/`metagenome`, special characters…)
- **Compute “unclassified” ratios** (per-sample percentage of reads in filtered ASVs) across taxonomic levels
- **Compute retained-taxa count ratios** (per-level percentage of ASV columns retained after filtering)
- **Plot**  
  - Well-classified proportions (1 – unclassified) per sample  
  - Retained-taxa (%) per level  
  - **Cumulative (stacked) barplots** very similar to QIIME2’s barplot output (Top N taxa + “Other” = 100%)

All outputs (CSV tables and PDF figures) are written into the same directory as your input `level-N.csv` files.

---

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/Newgenes1031/metabarcoding_taxonomy.git
