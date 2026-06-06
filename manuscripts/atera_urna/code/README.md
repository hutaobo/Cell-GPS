# Atera uRNA-COSTE code

This folder contains the manuscript-specific analysis and asset-generation scripts for the Atera uRNA-COSTE study.

## Scripts

- `atera_wta_breast_urna_segmentation_bias.py`: large transcript-level uRNA extraction, gene-level count summaries and uRNA-COSTE concordance analysis.
- `extract_spatial_example_data.py`: extraction of real uRNA and cell-centroid coordinate windows used for spatial validation.
- `make_manuscript_assets.py`: generation of manuscript-level summary tables, graphical abstract, Figure 1 and Supplementary Figure 1 from benchmark outputs.
- `make_spatial_evidence_figure.py`: generation of the multi-panel spatial evidence figure from extracted coordinate windows.
- `make_supplementary_table_workbooks.py`: formatting of Supplementary Tables S1-S8 as Excel workbooks.
- `make_docx.py`: conversion of the manuscript Markdown draft to a formatted Word document.

## Expected layout

The scripts are written to run from the repository root with this folder located at:

`manuscripts/atera_urna/code/`

The large benchmark output folder is expected at:

`benchmarking/atera_wta_breast_urna_segmentation_bias_results/`

The raw Atera transcript parquet file is not redistributed in this repository.

## Typical regeneration order

```powershell
python manuscripts\atera_urna\code\make_manuscript_assets.py
python manuscripts\atera_urna\code\make_supplementary_table_workbooks.py
python manuscripts\atera_urna\code\make_spatial_evidence_figure.py
python manuscripts\atera_urna\code\make_docx.py
```
