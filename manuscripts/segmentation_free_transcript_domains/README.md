# Segmentation-Free Spatial Domain Discovery from Transcript Points Alone

Analysis and figure-generation code for the study *"Segmentation-Free Spatial
Domain Discovery from Transcript Points Alone"* (manuscript under preparation).

The workflow discovers molecular tissue domains directly from transcript
coordinates, without cell segmentation: a full-gene di-sim screen ranks
informative genes, COSTE-derived cophenetic structure organizes them into gene
modules, and a streaming HistoSeg-style step assigns spatial domains to
individual transcript points.

## Dataset

Public 10x Genomics Atera In Situ Gene Expression preview, FFPE human breast
cancer: https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer

## Code (`code/`)

Analysis pipeline:
- `atera_breast_disim_coste_benchmark.py` — full-gene di-sim screen + COSTE module selection
- `atera_breast_cellfree_transcript_pipeline.py` — cell-free transcript-grid prototype
- `atera_breast_transcript_point_histoseg.py` — point-level HistoSeg domain assignment
- `atera_breast_point_vs_cell_histoseg_comparison.py` — point-vs-cell domain comparison
- `atera_breast_nature_pathology_figures.py` — H&E-aligned pathology figures

Figure / table generation:
- `_svg_relayout.py`, `_panel_d_evidence.py`, `relayout_supplementary_figures.py`,
  `relayout_main_figure_sources.py` — figure assembly / relayout
- `make_supp_fig6_sensitivity.py` — domain-number sensitivity figure
- `make_supp_table_s5.py` — sensitivity summary table
- `regen_overlap_heatmap.py` — point-vs-cell overlap heatmap

## Dependencies

Python 3.13 with numpy, scipy, scikit-learn, pandas, matplotlib, PyMuPDF (`fitz`),
Pillow, scikit-image, zarr, tifffile, psutil.

## Related software

- COSTE / Cell-GPS: https://github.com/hutaobo/Cell-GPS
- HistoSeg: https://github.com/hutaobo/HistoSeg

The manuscript text is not included in this directory; it will be available with
the published article.
