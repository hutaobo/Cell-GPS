# Cell-GPS

[![PyPI version](https://img.shields.io/pypi/v/Cell-GPS.svg)](https://pypi.org/project/Cell-GPS/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/cell-gps.svg)](https://anaconda.org/conda-forge/cell-gps)
[![Python versions](https://img.shields.io/pypi/pyversions/Cell-GPS.svg)](https://pypi.org/project/Cell-GPS/)
[![License](https://img.shields.io/pypi/l/Cell-GPS.svg)](LICENSE)
[![Upload Python Package](https://github.com/hutaobo/cellgps/actions/workflows/python-publish.yml/badge.svg)](https://github.com/hutaobo/cellgps/actions/workflows/python-publish.yml)

`Cell-GPS` is the Python package and reference implementation for Cophenetic Spatial Topology Embedding (COSTE), a spatial topology analysis framework for spatial omics data. New Python code should import `cellgps`; the historical `sfplot` import namespace is retained for backward compatibility.

This repository is maintained as both the installable Python package and the code companion for the Cell-GPS/COSTE bioRxiv preprint.

## Preprint and manuscript code

Cell-GPS/COSTE is described in the associated bioRxiv preprint:

> Long M, Hu T, Sountoulidis A, Samakovlis C, Nilsson M. Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics. bioRxiv. 2026. doi: [10.64898/2026.05.26.727847](https://doi.org/10.64898/2026.05.26.727847)

Versioned bioRxiv page: <https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1>

If you use the Python package, the Windows executable, the R companion package, or the manuscript figure/table code, please cite this preprint.

The code used to generate the preprint figures and supplementary tables is organized in [`Cell-GPS manuscript code/`](https://github.com/hutaobo/cellgps/tree/main/Cell-GPS%20manuscript%20code):

- [`main_figures/`](https://github.com/hutaobo/cellgps/tree/main/Cell-GPS%20manuscript%20code/main_figures): notebooks for main figure analyses.
- [`supplementary_figures/`](https://github.com/hutaobo/cellgps/tree/main/Cell-GPS%20manuscript%20code/supplementary_figures): notebooks for supplementary figure analyses.
- [`supplementary_tables/`](https://github.com/hutaobo/cellgps/tree/main/Cell-GPS%20manuscript%20code/supplementary_tables): notebooks for supplementary table analyses.

The notebooks are intentionally output-free and preserve the original manuscript data paths where those paths were required for reproduction. A detailed mapping from manuscript results to source code is available in [`docs/cellgps_science_manuscript_code_inventory.md`](https://github.com/hutaobo/cellgps/blob/main/docs/cellgps_science_manuscript_code_inventory.md).

## Package Names

- Python distribution: `Cell-GPS`
- Python import package: `cellgps`
- R package/repository: `cellgpsr`
- Windows executable: `cellgps.exe`

The Python package is hosted at `https://github.com/hutaobo/cellgps`. The R package is hosted separately at `https://github.com/hutaobo/cellgpsr`. The Windows single-file executable is distributed through Zenodo: <https://zenodo.org/records/17859173>.

## What Cell-GPS does

- Computes searcher-findee distance matrices from spatial coordinates and cell labels.
- Builds cophenetic distance matrices and StructureMap heatmaps from `AnnData` objects or plain `pandas` tables.
- Loads 10x Xenium outputs and prepares them for downstream spatial analysis.
- Provides Xenium loaders backed by `pyXenium.io`, including standard Xenium folders and table bundles with `cells.parquet`, official `*_cell_groups.csv`, and `cell_feature_matrix.h5`.
- Supports transcript-by-cell analysis for locating transcripts relative to cell types.
- Includes memory-optimized workflows for large datasets.
- Provides plotting utilities such as clustered heatmaps, circular dendrograms, and related summary figures.

## Repository layout

- `src/cellgps/`: recommended Python import namespace.
- `src/sfplot/`: legacy compatibility namespace and current implementation modules.
- `tests/`: package tests and smoke checks.
- `docs/`: Sphinx documentation.
- `Cell-GPS manuscript code/`: curated preprint figure and table notebooks.
- `sfplot-manuscript/`: legacy manuscript-specific notebooks, figures, and derived outputs.
- `benchmarking/`: benchmarking-related material.
- `segmentation_methods/`: supporting segmentation workflows.

## Installation

Install from PyPI:

```bash
pip install Cell-GPS
```

Install from conda-forge:

```bash
conda install -c conda-forge cell-gps
```

Install directly from GitHub:

```bash
pip install git+https://github.com/hutaobo/Cell-GPS.git
```

For local development or reviewer inspection:

```bash
git clone https://github.com/hutaobo/Cell-GPS.git
cd cellgps
pip install -e .
```

The package requires Python 3.9 or later.

## Quick start from a coordinate table

The minimal input is a table with spatial coordinates and a cell-type column.

```python
import pandas as pd
from cellgps import compute_cophenetic_distances_from_df, plot_cophenetic_heatmap

df = pd.DataFrame(
    {
        "x": [0, 1, 5, 6],
        "y": [0, 1, 5, 6],
        "celltype": ["A", "A", "B", "B"],
    }
)

row_coph, col_coph = compute_cophenetic_distances_from_df(
    df=df,
    x_col="x",
    y_col="y",
    celltype_col="celltype",
)

plot_cophenetic_heatmap(
    row_coph,
    matrix_name="row_coph",
    output_dir="output",
    output_filename="StructureMap_example.pdf",
    sample="Example",
)
```

## Quick start from Xenium output

```python
from cellgps import load_xenium_data, load_xenium_table_bundle, compute_cophenetic_distances_from_adata

# Standard Xenium folder through pyXenium.io
adata = load_xenium_data("/path/to/xenium/run", normalize=False)

# Explicit table-bundle route used for the Atera Xenium benchmark
adata = load_xenium_table_bundle("/path/to/xenium/run", normalize=False)

row_coph, col_coph = compute_cophenetic_distances_from_adata(
    adata,
    cluster_col="Cluster",
    output_dir="output",
)
```

## Useful public entry points

- `load_xenium_data`: load and preprocess Xenium data.
- `load_xenium_table_bundle`: load Xenium data from `cells.parquet` + `*_cell_groups.csv` + `cell_feature_matrix.h5` through `pyXenium.io`.
- `compute_cophenetic_distances_from_df`: compute structure matrices from a coordinate table.
- `compute_weighted_searcher_findee_distance_matrix_from_df`: weighted searcher-findee kernel for entity, pathway, or LR analysis.
- `compute_weighted_cophenetic_distances_from_df`: weighted StructureMap wrapper over the weighted kernel.
- `compute_cophenetic_distances_from_adata`: compute structure matrices from `AnnData`.
- `compute_entity_to_cell_topology`: generalize `t_and_c` from transcripts to arbitrary weighted entities.
- `compute_entity_structuremap`: build StructureMap-style topology among arbitrary weighted entities.
- `plot_cophenetic_heatmap`: generate StructureMap and related clustered heatmaps.
- `transcript_by_cell_analysis`: analyze transcript-to-cell spatial structure at scale.
- `ligand_receptor_topology_analysis`: score sender→receiver ligand-receptor candidates using topology, structure compatibility, and local contact.
- `ligand_receptor_target_consistency`: add a NicheNet-style downstream target-consistency layer.
- `compute_pathway_activity_matrix`: compute rank-based or weighted pathway activities per cell.
- `pathway_topology_analysis`: analyze pathway-to-cell and pathway-to-pathway spatial topology.
- `compute_cophenetic_distances_from_df_memory_opt`: memory-aware alternative for large tables.
- `plot_circular_dendrogram_pycirclize`: circular dendrogram visualization.

## Notes for reviewers

- The recommended Python import namespace is `cellgps`; the legacy `sfplot` namespace remains available for older scripts.
- The curated figure and table notebooks for the bioRxiv preprint are kept in `Cell-GPS manuscript code/`.
- `sfplot-manuscript/` is retained as a legacy manuscript workspace with additional notebooks, generated figure assets, and intermediate outputs.
- Raw experimental datasets are not bundled in this repository because of size and distribution constraints. The code expects standard spatial omics outputs such as Xenium folders or tabular coordinate inputs.
- When a `cellgps_tbc_formal_wta/results`-style directory is already available, the LR and pathway topology extensions are designed to reuse its `t_and_c_result_*.csv` and `StructureMap_table_*.csv` outputs as the preferred gene-level topology anchors before falling back to recomputation.
- Xenium loading depends on `pyXenium>=0.4.3`. Visium helpers remain optional through the separate `Cell-GPS[visium]` extra.
- A short repository walkthrough is available in [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md).

## Documentation

Sphinx documentation sources are available in `docs/`.

## Citation

If you use Cell-GPS or reuse the manuscript analysis code, please cite:

> Long M, Hu T, Sountoulidis A, Samakovlis C, Nilsson M. Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics. bioRxiv. 2026. doi: [10.64898/2026.05.26.727847](https://doi.org/10.64898/2026.05.26.727847)

```bibtex
@article{long2026cellgps,
  title = {Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics},
  author = {Long, Mengping and Hu, Taobo and Sountoulidis, Alexandros and Samakovlis, Christos and Nilsson, Mats},
  journal = {bioRxiv},
  year = {2026},
  doi = {10.64898/2026.05.26.727847},
  url = {https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1}
}
```

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
