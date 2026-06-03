# Reviewer Guide

This repository contains the code used to develop the Cell-GPS workflow described in the manuscript.

## Naming

- Python distribution: `Cell-GPS`
- Conda-forge distribution: `cell-gps`
- Python import package: `cellgps`
- Legacy Python compatibility namespace: `sfplot` (not a separate distribution; retained for existing scripts)
- Python GitHub repository: `https://github.com/hutaobo/cellgps`
- R package/repository: `cellgpsr`, hosted at `https://github.com/hutaobo/cellgpsr`
- Windows executable: `cellgps.exe`, latest open v2 DOI `https://doi.org/10.5281/zenodo.19482685`; version-series route `https://zenodo.org/records/17859173`

## Where to start

If you only have a few minutes, read the files in this order:

1. `README.md`
2. `src/sfplot/analysis/searcher_findee_score.py`
3. `src/sfplot/preprocessing/data_processing.py`
4. `src/sfplot/analysis/tbc_analysis.py`

## What each area contains

- `src/cellgps/`
  Recommended Python import namespace.
- `src/sfplot/`
  Legacy compatibility namespace that currently hosts implementation modules; new code should import through `cellgps`.
- `docs/`
  Lightweight package documentation.
- `Cell-GPS manuscript code/`
  Curated manuscript figure and supplementary table notebooks.
- `examples/`
  Compact package usage examples and small example data files.
- `packaging/`
  Historical packaging recipes and Windows executable build assets.

## Minimal install

```bash
git clone https://github.com/hutaobo/cellgps.git
cd cellgps
pip install -e .
```

## Manuscript validation scope

The manuscript-validated scope is the COSTE/SSS workflow and the figure/table
analyses mapped in `Cell-GPS manuscript code/` and
`docs/cellgps_science_manuscript_code_inventory.md`. Ligand-receptor topology,
pathway topology, Visium helpers, GUI entry points and other convenience APIs
are included for reuse and development, but should be treated as optional or
exploratory unless a manuscript notebook or documentation page explicitly maps
them to a reported analysis.

## Minimal input contract

For the plain DataFrame workflow, the minimal required columns are:

- `x`
- `y`
- `celltype`

With those columns you can run:

```python
from cellgps import compute_cophenetic_distances_from_df

row_coph, col_coph = compute_cophenetic_distances_from_df(
    df,
    x_col="x",
    y_col="y",
    celltype_col="celltype",
)
```

## Xenium workflow

The Xenium workflow expects a standard Xenium output directory and uses:

- `load_xenium_data`
- `load_xenium_table_bundle`
- `compute_cophenetic_distances_from_adata`
- `transcript_by_cell_analysis`

## Loader note

- `load_xenium_data` uses `pyXenium.io` for standard Xenium output directories.
- `load_xenium_table_bundle` uses the same `pyXenium.io` reader with explicit `cells.parquet`, official `*_cell_groups.csv`, and `cell_feature_matrix.h5` artifacts.

## Precomputed anchors

- The LR and pathway topology extensions can reuse an existing `cellgps_tbc_formal_wta/results` directory.
- When such a directory is available, `t_and_c_result_*.csv` and `StructureMap_table_*.csv` are treated as the preferred gene-level topology anchors before any recomputation.

## Reproducibility note

Raw datasets are not included in this repository because of size and data-distribution constraints. The repository is intended to expose the analysis code, plotting functions, and manuscript-facing workflow structure used in the study.
