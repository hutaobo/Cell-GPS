# Atera uRNA/COSTE Nature Methods draft package

This folder contains a Nature Methods style draft package for the Atera WTA breast uRNA segmentation-bias analysis.

## Files

- `manuscript_draft.md`: main Article-style manuscript draft.
- `supplementary_information.md`: supplementary results, table descriptions and limitations.
- `cover_letter_draft.md`: short cover letter draft for Nature Methods.
- `nature_methods_compliance_checklist.md`: checked format requirements and current evidence gaps.
- `make_manuscript_assets.py`: script that regenerates all manuscript tables and figures from the committed benchmark output.
- `extract_spatial_example_data.py`: A100-side extractor for real uRNA/cell-coordinate windows used in Figure 5.
- `make_spatial_evidence_figure.py`: builds the Nature-style multi-panel spatial evidence figure.
- `paper_level_summary.md`: machine-generated numerical summary.
- `spatial_examples/`: real extracted uRNA and cell-centroid points for the Figure 5 examples.
- `figures/`: graphical abstract and draft main/extended figures.
- `tables/`: derived tables used by the draft.

## Source analysis

The manuscript assets are derived from:

`benchmarking/atera_wta_breast_urna_segmentation_bias_results/`

The source analysis was run on the A100 server using high-quality Atera WTA breast transcripts and the previously generated all-transcript COSTE/SSS reference table.

## Regeneration

Run from the repository root:

```powershell
python manuscripts\atera_urna_nature_methods\make_manuscript_assets.py
python manuscripts\atera_urna_nature_methods\make_spatial_evidence_figure.py
```

This does not rerun the large A100 transcript-coordinate analysis. It only converts the committed CSV, JSON and extracted spatial-window outputs into publication-level tables and figures.

To refresh the real spatial example windows on A100, copy `extract_spatial_example_data.py` to the A100 output folder and run:

```bash
/home/taobo.hu/miniconda3/envs/pyXenium/bin/python extract_spatial_example_data.py \
  --output-dir /data/taobo.hu/atera_urna_segmentation_bias_20260604/spatial_examples_nature \
  --window-um 850 --max-urna-points 7000 --max-other-cells 3500
```
