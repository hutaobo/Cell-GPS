# Cell-GPS

[![PyPI version](https://img.shields.io/pypi/v/Cell-GPS.svg)](https://pypi.org/project/Cell-GPS/)
[![Documentation Status](https://readthedocs.org/projects/cell-gps/badge/?version=latest)](https://cell-gps.readthedocs.io/en/latest/?badge=latest)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/cell-gps.svg)](https://anaconda.org/conda-forge/cell-gps)
[![Python versions](https://img.shields.io/pypi/pyversions/Cell-GPS.svg)](https://pypi.org/project/Cell-GPS/)
[![License](https://img.shields.io/pypi/l/Cell-GPS.svg)](LICENSE)
[![Upload Python Package](https://github.com/hutaobo/Cell-GPS/actions/workflows/python-publish.yml/badge.svg)](https://github.com/hutaobo/Cell-GPS/actions/workflows/python-publish.yml)

`Cell-GPS` is the Python package and reference implementation for Cophenetic Spatial Topology Embedding (COSTE), a spatial topology analysis framework for spatial omics data.

This repository is maintained as both the installable Python package and the code companion for the Cell-GPS/COSTE bioRxiv preprint.

## Preprint and manuscript code

Cell-GPS/COSTE is described in the associated bioRxiv preprint:

> Long M, Hu T, Sountoulidis A, Samakovlis C, Nilsson M. Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics. bioRxiv. 2026. doi: [10.64898/2026.05.26.727847](https://doi.org/10.64898/2026.05.26.727847)

Versioned bioRxiv page: <https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1>

If you use the Python package, the Windows executable, the R companion package, or the manuscript figure/table code, please cite this preprint.

The code used to generate the preprint figures and supplementary tables is organized in [`Cell-GPS manuscript/`](https://github.com/hutaobo/Cell-GPS/tree/main/Cell-GPS%20manuscript):

- [`main_figures/`](https://github.com/hutaobo/Cell-GPS/tree/main/Cell-GPS%20manuscript/main_figures): notebooks for main figure analyses.
- [`supplementary_figures/`](https://github.com/hutaobo/Cell-GPS/tree/main/Cell-GPS%20manuscript/supplementary_figures): notebooks for supplementary figure analyses.
- [`supplementary_tables/`](https://github.com/hutaobo/Cell-GPS/tree/main/Cell-GPS%20manuscript/supplementary_tables): notebooks for supplementary table analyses.

The notebooks are intentionally output-free and preserve the original manuscript data paths where those paths were required for reproduction. A detailed mapping from manuscript results to source code is available in [`docs/cellgps_science_manuscript_code_inventory.md`](https://github.com/hutaobo/Cell-GPS/blob/main/docs/cellgps_science_manuscript_code_inventory.md).

## COSTE extension manuscripts

Beyond the main Cell-GPS/COSTE preprint above, the [`manuscripts/`](https://github.com/hutaobo/Cell-GPS/tree/main/manuscripts) directory collects the analysis and figure-generation code for follow-up manuscripts that build on the COSTE method. Each subfolder is a self-contained study with its own README describing the dataset, scripts, and reproduction steps:

- [`manuscripts/atera_urna/`](https://github.com/hutaobo/Cell-GPS/tree/main/manuscripts/atera_urna): Atera uRNA-COSTE study — segmentation-bias check using unassigned RNA (uRNA, transcripts not assigned to a segmented cell) and uRNA-COSTE concordance analysis.
- [`manuscripts/segmentation_free_transcript_domains/`](https://github.com/hutaobo/Cell-GPS/tree/main/manuscripts/segmentation_free_transcript_domains): segmentation-free spatial domain discovery from transcript points alone.

For details on any extension manuscript, read the README inside its subfolder.

## Package Names

- Python distribution: `Cell-GPS`
- Conda-forge distribution: `cell-gps`
- Python import package: `cellgps`
- R package/repository: `cellgpsr`
- Windows executable: `cellgps.exe`

The Python package is hosted at `https://github.com/hutaobo/Cell-GPS`. The R package is hosted separately at `https://github.com/hutaobo/cellgpsr`. The Windows single-file executable is distributed through Zenodo; use the latest open v2 DOI <https://doi.org/10.5281/zenodo.19482685> or the version-series route <https://zenodo.org/records/17859173>.

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
- `src/cellgps/pp`, `src/cellgps/tl`, `src/cellgps/pl`: scverse-style aliases for preprocessing, analysis, and plotting APIs.
- `tests/`: package tests and smoke checks.
- `docs/`: Sphinx documentation.
- `docs/project/`: project notes, changelog, authors, and reviewer guide.
- `Cell-GPS manuscript/`: curated preprint figure and table notebooks for the main Cell-GPS/COSTE bioRxiv preprint.
- `manuscripts/`: analysis code for follow-up COSTE extension manuscripts, one self-contained study per subfolder (see each subfolder's README).
- `examples/`: compact usage examples and small example data files.
- `packaging/conda-recipe/`: archived local conda recipe retained for reference.
- `packaging/pyinstaller/`: Windows executable build scripts and PyInstaller assets.

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
- `ligand_receptor_topology_analysis`: score sender->receiver ligand-receptor candidates using topology, structure compatibility, and local contact.
- `ligand_receptor_target_consistency`: add a NicheNet-style downstream target-consistency layer.
- `compute_pathway_activity_matrix`: compute rank-based or weighted pathway activities per cell.
- `pathway_topology_analysis`: analyze pathway-to-cell and pathway-to-pathway spatial topology.
- `compute_cophenetic_distances_from_df_memory_opt`: memory-aware alternative for large tables.
- `plot_circular_dendrogram_pycirclize`: circular dendrogram visualization.

## Validation scope

The manuscript-validated scope is the COSTE/SSS workflow and the figure/table analyses mapped in `Cell-GPS manuscript/` and `docs/cellgps_science_manuscript_code_inventory.md`. Ligand-receptor topology, pathway topology, Visium helpers, GUI entry points and other convenience APIs are included for reuse and development, but should be treated as optional or exploratory unless a manuscript notebook or documentation page explicitly maps them to a reported analysis.

## Notes for reviewers

- The curated figure and table notebooks for the bioRxiv preprint are kept in `Cell-GPS manuscript/`.
- Raw experimental datasets are not bundled in this repository because of size and distribution constraints. The code expects standard spatial omics outputs such as Xenium folders or tabular coordinate inputs.
- Conda-forge packages the upstream Python distribution as `cell-gps`.
- When a `cellgps_tbc_formal_wta/results`-style directory is already available, the LR and pathway topology extensions are designed to reuse its `t_and_c_result_*.csv` and `StructureMap_table_*.csv` outputs as the preferred gene-level topology anchors before falling back to recomputation.
- Xenium loading depends on `pyXenium>=0.4.3`. Visium helpers remain optional through the separate `Cell-GPS[visium]` extra.
- A short repository walkthrough is available in [docs/project/REVIEWER_GUIDE.md](docs/project/REVIEWER_GUIDE.md).

## Documentation

Read the Docs documentation is available at <https://cell-gps.readthedocs.io/en/latest/>.
The manuscript-focused pages introduce the bioRxiv preprint, explain how to use the curated figure/table notebooks, and map each figure and supplementary table to its GitHub code location.

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
