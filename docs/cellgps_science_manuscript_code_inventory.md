# CellGPS Science manuscript code inventory

Date: 2026-05-28

Manuscript target:

`Y:\long\projects\Packagedevelopment\CellGPS\CellGPS manuscript for Science`

Primary manuscript inspected:

`Y:\long\projects\Packagedevelopment\CellGPS\CellGPS manuscript for Science\manuscript new revision.docx`

## Search scope

I inspected the manuscript text, extracted figures, captions, local repository files, the mapped Y drive, relevant local user folders, PDC, and A100.

Visible Windows drives during the search were `C:`, `D:`, `E:`, `F:`, `G:`, `H:`, `X:`, `Y:`, and `Z:`. The most relevant local/Y locations were searched by content with source-only notebook parsing to avoid false hits from embedded image output. A broader filename search over the visible drives was also started; it timed out after producing partial hits, but the additional useful hits were only local cached/output copies under `C:\Users\taobo.hu\Downloads` and WPS cache.

Remote hosts checked:

- PDC: `pdc-cmd` / `dardel.pdc.kth.se`, home `/cfs/klemming/home/h/hutaobo`. No related CellGPS manuscript code was found in the searched home tree.
- A100: `sscb-a100.scilifelab.se`, home `/home/taobo.hu`, data roots `/data/taobo.hu` and `/mnt/taobo.hu`. This host contains the strongest missing evidence for the synthetic benchmark and DST-GNN.

Confidence labels:

- Direct: code explicitly loads the dataset, computes COSTE/SSS/StructureMap, or writes named manuscript outputs.
- Strong inference: code/output names and figure content match the manuscript panel, but the final assembled figure was likely made manually from generated PDFs/PNGs.
- Supporting: package/helper code or copied outputs used by direct notebooks/scripts.

## High-level conclusion

The Science manuscript directory itself is mainly a submission/assembly folder. It contains Word/PDF manuscripts, final figure PDFs/PNGs, supplementary files, and `sfplot_COSTE_README.txt`, but not the main computational notebooks/scripts.

The code used for the manuscript is distributed across:

- Current local working tree: `D:\GitHub\sfplot`
- Ignored local manuscript workspace: `D:\GitHub\sfplot\sfplot-manuscript`
- Original Y-drive CellGPS project: `Y:\long\projects\Packagedevelopment\CellGPS`
- Lung fibrosis/Vannan project: `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis`
- TNBC/Ali project: `Y:\long\publication_datasets\Ali_TNBC_2023\NTPublic\sfplot`
- 10x Xenium t-by-c outputs: `Y:\long\10X_datasets\Xenium\Xenium_5K\t_by_c_result`
- A100 synthetic/DST-GNN workspaces under `/data/taobo.hu`

## Core method/package code

Direct/supporting code used across figures:

- `D:\GitHub\sfplot\src\sfplot\analysis\searcher_findee_score.py`
  - Contains `compute_cophenetic_distances_from_adata`, `compute_searcher_findee_distance_matrix_from_df`, `compute_cophenetic_from_distance_matrix`, `compute_cophenetic_distances_from_df`, and `plot_cophenetic_heatmap`.
- `D:\GitHub\sfplot\src\sfplot\analysis\compute_cophenetic_distances_from_df_memory_opt.py`
  - Memory-aware implementation used for large transcript/cell merged point clouds.
- `D:\GitHub\sfplot\src\sfplot\analysis\tbc_analysis.py`
- `D:\GitHub\sfplot\src\sfplot\analysis\tbc_analysis_serial.py`
  - Transcript-by-cell analysis: loads Xenium data/transcripts, computes cell-type StructureMap, then recomputes gene/transcript rows.
- `D:\GitHub\sfplot\src\sfplot\analysis\compute_col_dendrogram_scores.py`
- `D:\GitHub\sfplot\src\sfplot\plotting\circle_heatmap.py`
- `D:\GitHub\sfplot\src\sfplot\plotting\circular_dendrogram.py`
- Original package copy: `Y:\long\projects\Packagedevelopment\CellGPS\sfplot\src\sfplot`
- R implementation/helper copy:
  - `Y:\long\projects\Packagedevelopment\CellGPS\sfplotR_manuscript\R`
  - `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\R`

## Figure/result mapping

### Main Figure 1 and Supplementary Figures 1-6

Topic from manuscript: synthetic modular, nested, and random spatial benchmarks; comparison of COSTE, Squidpy, Giotto, and ANE.

Direct code on A100:

- `/data/taobo.hu/projects/mouse_pup/benchmarking_complete_english.ipynb`
  - Generates modular/nested/random synthetic datasets.
  - Computes COSTE/SSS with `compute_cophenetic_distances_from_df`.
  - Runs Squidpy NE, Giotto proximity enrichment, and ANE.
  - Writes method heatmaps such as `Zscore_modular.pdf`, `Zscore_nested.pdf`, `Zscore_random.pdf`, `Giotto_*`, `ANE_*`.
- `/data/taobo.hu/projects/mouse_pup/benchmarking_complete_modular.ipynb`
  - Generates modular variants and writes `benchmarking_StructureMap_of_synthetic_modular_1.pdf`, `_2.pdf`, `_3.pdf`, plus method-specific modular heatmaps.
- `/data/taobo.hu/projects/mouse_pup/benchmarking_complete_nested.ipynb`
  - Generates the 12 nested variants named in the supplement: `nested_radius_1-4`, `nested_thickness_1-4`, `nested_centers_1-4`.
  - Writes `benchmarking_StructureMap_of_nested_*.pdf`, `COSTE_all_nested.pdf`, `Squidpy_Z_all_nested.pdf`, `Giotto_nested_12panel.pdf`, `ANE_nested_12panel.pdf`, and `twelve_nested_patterns.pdf`.
- `/data/taobo.hu/projects/mouse_pup/benchmarking_synthetic.ipynb`
  - Earlier/general synthetic benchmark notebook for six cell types A-F.
- `/data/taobo.hu/projects/mouse_pup/coste_sss_benchmark_synthetic.ipynb`
  - Related simplified benchmark notebook with COSTE-like SSS, permutation NE, analytical NE, and random baseline.

Outputs on A100:

- `/data/taobo.hu/projects/mouse_pup/COSTE_all_nested.pdf`
- `/data/taobo.hu/projects/mouse_pup/Squidpy_Z_all_nested.pdf`
- `/data/taobo.hu/projects/mouse_pup/Giotto_nested_12panel.pdf`
- `/data/taobo.hu/projects/mouse_pup/ANE_nested_12panel.pdf`
- `/data/taobo.hu/projects/mouse_pup/benchmarking_StructureMap_of_synthetic_*.pdf`
- `/data/taobo.hu/projects/mouse_pup/benchmarking_StructureMap_of_nested_*.pdf`

Local supporting/cached copies:

- `C:\Users\taobo.hu\Downloads\synthetic_nested_benchmark.ipynb`
  - Only a minimal/partial notebook was found locally, not the canonical source.
- `C:\Users\taobo.hu\Downloads\Zscore_nested.pdf`
- `C:\Users\taobo.hu\Downloads\nested_pattern.pdf`
- WPS cache contains output copies including `COSTE_all_nested.pdf` and `twelve_nested_patterns.pdf`.

Confidence: Direct for A100 notebooks and outputs. Strong inference that final Figure 1 and Supplementary Figures 1-6 were assembled from these generated PDFs.

### Main Figure 2 and Supplementary Figures 7-9

Topic from manuscript: neonatal mouse pup Xenium, hierarchical spatial structures S1-S6, retina/keratinocyte/fibroblast panels, and comparison with Squidpy/Giotto/ANE.

Direct/local code:

- `D:\GitHub\sfplot\sfplot-manuscript\figures\fig. 1.ipynb`
  - Despite the filename, this notebook loads `/mnt/taobo.hu/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Mouse_Pup_FFPE_outs`.
  - Computes `StructureMap_of_Mouse_Pup_FFPE.pdf`.
  - Generates mouse pup spatial scatter panels and Squidpy neighborhood enrichment panels with different `n_neighs`/radius settings.
- `D:\GitHub\sfplot\sfplot-manuscript\t_and_c\plot.ipynb`
- `D:\GitHub\sfplot\sfplot-manuscript\t_and_c\分细胞的CellGPS.ipynb`
  - Additional mouse pup StructureMap/spatial panel notebooks.
- `D:\GitHub\sfplot\sfplot-manuscript\figures\suppl. fig. 4 mouse brain.ipynb`
  - Mouse brain comparison/t-by-c notebook, relevant to the supplementary 5K Xenium comparison panels rather than the main mouse pup panel.
- `D:\GitHub\sfplot\sfplot-manuscript\t_and_c\plot_brain.ipynb`

Direct Y-drive code:

- `Y:\long\projects\Packagedevelopment\CellGPS\CellGPS manuscript\figures\meta-structures.R`
  - Reads `Y:/long/10X_datasets/Xenium/Xenium_5K/t_by_c_result/t_by_c_Xenium_Prime_Mouse_Pup_FFPE_outs/`
  - Reads `t_and_c_result_Xenium_Prime_Mouse_Pup_FFPE_outs.csv` and `transcript_table_percentage_Xenium_Prime_Mouse_Pup_FFPE_outs.csv`.
  - Writes `gene_module Mouse Pup.csv`.
- `Y:\long\projects\Packagedevelopment\CellGPS\CellGPS manuscript\figures\suppl. fig. 4.R`
  - Reads mouse brain t-by-c result tables and `C:/Users/taobo.hu/Downloads/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv`.

Data/output directories:

- `Y:\long\10X_datasets\Xenium\Xenium_5K\Xenium_Prime_Mouse_Pup_FFPE_outs`
- `Y:\long\10X_datasets\Xenium\Xenium_5K\t_by_c_result\t_by_c_Xenium_Prime_Mouse_Pup_FFPE_outs`
- `D:\GitHub\sfplot\sfplot-manuscript\t_and_c\t_by_c_Xenium_Prime_Mouse_Pup_FFPE_outs`
- `/data/taobo.hu/projects/mouse_pup`

Performance/comparison code for Supplementary Table 1 and Supplementary Figures 8-9:

- `/data/taobo.hu/projects/mouse_pup/benchmarking_mouse_pup.ipynb`
  - Benchmarks COSTE, Squidpy NE, Giotto, and ANE on the mouse pup dataset.
  - Tracks CPU time, wall time, peak Python memory, peak RSS, page faults, context switches, and throughput.
- `/data/taobo.hu/projects/mouse_pup/tmp_coste_perf.py`
- `/data/taobo.hu/projects/mouse_pup/tmp_squidpy_ne.py`

Confidence: Direct for StructureMap, t-by-c, and benchmark computation. Strong inference for final assembled Figure 2/Supplementary Figures 7-9 panels.

### Main Figure 3 and Supplementary Figure 11

Topic from manuscript: systemic sclerosis/pulmonary fibrosis, healthy vs fibrotic topology, AT1/AT2/capillary SSS, TRU remodeling score, regional TRS.

Direct Y-drive code:

- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\notebook\Expression Distance Similarity.ipynb`
  - Loads Xenium samples from `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Xenium`.
  - Calls `load_xenium_data`, `compute_cophenetic_distances_from_adata`, and `plot_cophenetic_heatmap`.
  - Writes `cophenetic_distances_searcher_D_score_in_all_samples.csv`.
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\example\plot.R`
  - Reads `GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds`.
  - Loops through HD/LA/MA samples.
  - Calls `compute_cophenetic_distances_from_df`.
  - Writes `cophenetic_distances_searcher_D_score_in_all_samples.csv`.
  - Plots AT1/AT2/Capillary panels and airway/RASC/myofibroblast panels.
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\example\分区域的Structure Score.R`
  - Reads `Y:/long/publication_datasets/Vannan_2023_Lung_Fibrosis/Xenium/HE_annotations/cells_partitioned_by_annotation.csv`.
  - Computes regional AT2-to-Capillary SSS and writes `分区域的Structure Score.csv`.
  - Saves `fig. 1h.pdf`.
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\notebook\regional_plot.ipynb`
  - Computes hexagon-level regional AT2-Capillary values for `VUILD110LA`.
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\R\TRU_remodelling_activity.R`
  - Uses sample `VUILD110LA`.
  - Computes a TRU remodeling activity overlay using AT2 and Capillary distances.
  - Saves `C:/Users/taobo.hu/Downloads/TRU_remodelling_activity.pdf`.

Supporting local/output copies:

- `C:\Users\taobo.hu\Downloads\Vannan_2023_Lung_Fibrosis.pdf`
- `C:\Users\taobo.hu\Downloads\StructureMap_of_SSc_2_1_1_raw.pdf`
- `C:\Users\taobo.hu\Downloads\t_by_c_SSc_2_1_1_raw`

Confidence: Direct for matrix generation, AT1/AT2/Capillary SSS, regional TRS, and TRU overlay.

### Main Figure 4

Topic from manuscript: fibrosis progression, pathological structures, airway remodeling, and DST-GNN.

Direct code for COSTE/SSS inputs:

- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\notebook\Expression Distance Similarity.ipynb`
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\example\plot.R`
- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR\R\Structure-Based Niche Analysis.R`
  - Reads the same `GSE250346` Seurat RDS and analyzes sample `VUILD110LA`.
  - Writes `clara_cluster_assignment.csv` and `Clustering_Comparison.pdf`.

Direct/recovered DST-GNN code:

- `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\notebook\GNN modelling.ipynb`
  - Reads `../Rcode/sfplotR/example/cophenetic_distances_searcher_D_score_in_all_samples.csv`.
  - Uses PyTorch/PyTorch Geometric, builds temporal graph data, defines `DiffusionGNN`, trains the model, and ranks node/edge changes.
- `/data/taobo.hu/tmp/dst_gnn_public`
  - Cleaned public DST-GNN implementation reconstructed from manuscript notebooks.
- `/data/taobo.hu/tmp/dst_gnn_release_repo`
  - Release repository with reproducibility metadata/checksum support.
- `/data/taobo.hu/tmp/dst_gnn_public/scripts/run_dst_gnn.py`
- `/data/taobo.hu/tmp/dst_gnn_release_repo/scripts/run_dst_gnn.py`
  - CLI input: a flattened CSV with `row`, `column`, `value`, `sample`, `group`.
  - Expected input matches `cophenetic_distances_searcher_D_score_in_all_samples.csv`.

DST-GNN output locations:

- `/data/taobo.hu/tmp/dst_gnn_public_outputs`
- `/data/taobo.hu/tmp/dst_gnn_public_outputs_explainer`
  - Contains `observed_T1_sss.csv`, `observed_T2_sss.csv`, `observed_T3_sss.csv`, predicted next-stage SSS matrices, `top_dynamic_nodes.csv`, `top_edge_changes.csv`, `training_history.csv`, and optional explainer output.
  - Summary indicates 47 nodes, stages `T1/T2/T3`, and sample counts `10/15/20`.

Confidence: Direct for DST-GNN code and outputs. Strong inference that the final Figure 4 DST-GNN panels were assembled from the ranked node/edge outputs plus manuscript-specific plotting.

### Main Figure 5 and Supplementary Figure 12

Topic from manuscript: segmentation-free pleural transcript/cell StructureMap, circular gene hierarchy, module genes, IL10/CCL4, CCL19/CCR7, and landmark clusters.

Direct local code:

- `D:\GitHub\sfplot\sfplot-manuscript\suppl. fig. 8 new.ipynb`
  - Uses `/data/taobo.hu/ILD/SSc_2_1_1_raw`.
  - Builds `t_by_c_SSc_2_1_1_Landmark`.
  - Computes segment-free merged transcript/cell StructureMap using `compute_cophenetic_distances_from_df`.
  - Writes `StructureMap_of_SSc_2_1_1_Landmark_Segment_free.pdf`.
  - Writes `StructureMap_table_SSc_2_1_1_Landmark_Segment_free.csv`.
  - Generates circular dendrograms and SSc gene spatial panels.
- `D:\GitHub\sfplot\sfplot-manuscript\fig. 3 new.ipynb`
  - Loads `/data/taobo.hu/ILD/SSc_2_1_1_raw`.
  - Extracts cells/transcripts from MEX/zarr inputs.
  - Produces IL10/CCL4 and selected cluster spatial panels, including `SSc_2_1_1_clusters3_6_26__IL10_CCL4.pdf`.
- `D:\GitHub\sfplot\sfplot-manuscript\t_by_c_SSc_2_1_1_Landmark`
  - Contains `t_and_c_result_SSc_2_1_1.csv`, `StructureMap_table_SSc_2_1_1.csv`, `StructureMap_table_SSc_2_1_1_Landmark_Segment_free.csv`, `GeneCounts_byCluster_SSc_2_1_1.csv`, and PDFs.

Direct/supporting A100 locations:

- `/data/taobo.hu/projects/ILD/SSc_2_1_1_raw`
- `/data/taobo.hu/projects/ILD/t_by_c_SSc_2_1_1_raw`
- `/data/taobo.hu/projects/sfplot-manuscript/t_by_c_SSc_2_1_1_Landmark`
- `/data/taobo.hu/projects/sfplot-manuscript/StructureMap_of_SSc_2_1_1_Landmark_Segment_free.pdf`
- `/data/taobo.hu/projects/sfplot-manuscript/SSc_2_1_1_clusters3_6_26__IL10_CCL4.pdf`

Supporting package code:

- `D:\GitHub\sfplot\src\sfplot\analysis\tbc_analysis.py`
- `D:\GitHub\sfplot\src\sfplot\analysis\tbc_analysis_serial.py`
- `D:\GitHub\sfplot\src\sfplot\analysis\compute_cophenetic_distances_from_df_memory_opt.py`
- `D:\GitHub\sfplot\src\sfplot\plotting\circular_dendrogram.py`
- `D:\GitHub\sfplot\src\sfplot\plotting\circle_heatmap.py`

Confidence: Direct for transcript/cell StructureMap, circular hierarchy, and gene spatial panels. Strong inference that final Figure 5 was assembled from these outputs.

### Supplementary Figure 10

Topic from manuscript: human lymph node StructureMap and spatial map.

Direct Y-drive code:

- `Y:\long\projects\Packagedevelopment\CellGPS\CellGPS manuscript\figures\LN点图.R`
  - Reads `Y:/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/cells.csv.gz`.
  - Reads `Y:/long/10X_datasets/Xenium/Xenium_5K/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_cell_types.csv`.
  - Sources `Y:/long/projects/NPC/rcode/compute_cophenetic_distances_from_df.R`.
  - Sources `Y:/long/projects/Pakagedevelopment/CellGPS/sfplotR/R/plot_cophenetic_heatmap.R` (note the path typo `Pakagedevelopment` appears in the script).
  - Saves `~/Downloads/cells_raster.pdf`.

Data/output directories:

- `Y:\long\10X_datasets\Xenium\Xenium_5K\Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs`
- `Y:\long\10X_datasets\Xenium\Xenium_5K\t_by_c_result\t_by_c_Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs`

Local notebook note:

- `D:\GitHub\sfplot\sfplot-manuscript\淋巴结.ipynb` exists but contained no relevant source cells in the parsed copy.

Confidence: Direct for the lymph node spatial/StructureMap script; output assembly likely manual.

### Supplementary Figures 13-14

Topic from manuscript: TNBC spatial biomarkers, pCR/RD comparisons, treatment arms/phases, and 12 subgroup heatmaps.

Direct Y-drive code:

- `Y:\long\publication_datasets\Ali_TNBC_2023\NTPublic\sfplot\sfplot_TNBC.ipynb`
  - Starts in `~/moldia/long/publication_datasets/Ali_TNBC_2023/NTPublic/`.
  - Sources `code/header.R`.
  - Loads curated cells with `getCells(curated = TRUE, allCells = TRUE)`.
  - Loads `data/derived/clinical.csv`.
  - Sources Vannan/sfplotR helper functions for `compute_cophenetic_distances_from_df` and `plot_cophenetic_heatmap`.
  - Writes per-patient and subgroup Searcher/Findee heatmaps.
  - Computes marker-pair Wilcoxon tests and writes `wilcox_test_results.csv` / `wilcox_results_TNBC.csv`.

Direct local copies:

- `D:\GitHub\sfplot\sfplot-manuscript\fig. 4.ipynb`
  - Local/A100 copy of TNBC subgroup and biomarker boxplot analysis.
  - Writes `PD-L1_APCs_cophenetic_heatmap.pdf`, `Baseline_C&I_pCR_vs_RD.pdf`, `Baseline_C_pCR_vs_RD.pdf`, and `wilcox_results_TNBC.csv`.
- `D:\GitHub\sfplot\sfplot-manuscript\suppl. fig. 9.ipynb`
  - Per-patient and subgroup TNBC heatmap generation.
- `D:\GitHub\sfplot\sfplot-manuscript\Findee`
- `D:\GitHub\sfplot\sfplot-manuscript\Searcher`
- `D:\GitHub\sfplot\sfplot-manuscript\output`

A100 mirror:

- `/data/taobo.hu/Ali_TNBC_2023/NTPublic/sfplot/sfplot_TNBC.ipynb`
- `/data/taobo.hu/projects/sfplot-manuscript/fig. 4.ipynb`
- `/data/taobo.hu/projects/sfplot-manuscript/suppl. fig. 9.ipynb`

Confidence: Direct for the TNBC heatmaps, pCR/RD marker-pair statistics, and Supplementary Figures 13-14.

## Supplementary tables

### Supplementary Table 1

Topic: runtime/memory performance for mouse pup.

Direct code:

- `/data/taobo.hu/projects/mouse_pup/benchmarking_mouse_pup.ipynb`
- `/data/taobo.hu/projects/mouse_pup/tmp_coste_perf.py`
- `/data/taobo.hu/projects/mouse_pup/tmp_squidpy_ne.py`

The notebook explicitly benchmarks CPU time, wall time, peak Python memory, peak RSS, page faults, context switches, and throughput for COSTE, Squidpy, Giotto, and ANE.

Confidence: Direct for benchmarking computations. The exact table file in the manuscript directory appears assembled separately.

### Supplementary Table 2

Topic: transcript-by-cell SSS in systemic sclerosis lung fibrosis.

Direct/supporting code and outputs:

- `D:\GitHub\sfplot\sfplot-manuscript\fig. 3 new.ipynb`
- `D:\GitHub\sfplot\sfplot-manuscript\suppl. fig. 8 new.ipynb`
- `D:\GitHub\sfplot\sfplot-manuscript\t_by_c_SSc_2_1_1_Landmark`
- `/data/taobo.hu/projects/ILD/t_by_c_SSc_2_1_1_raw`
- `C:\Users\taobo.hu\Downloads\t_by_c_SSc_2_1_1_raw`

Related table-generation helper:

- `D:\GitHub\sfplot\sfplot-manuscript\t_and_c\generate_low_sss_extremes_tables.py`

Confidence: Direct for the underlying t-by-c SSS tables; strong inference for final Excel table assembly.

### Supplementary Table 3

Topic: transcript SSS using landmark clusters.

Direct/supporting code and outputs:

- `D:\GitHub\sfplot\sfplot-manuscript\suppl. fig. 8 new.ipynb`
- `D:\GitHub\sfplot\sfplot-manuscript\t_by_c_SSc_2_1_1_Landmark\StructureMap_table_SSc_2_1_1_Landmark_Segment_free.csv`
- `D:\GitHub\sfplot\sfplot-manuscript\t_by_c_SSc_2_1_1_Landmark\t_and_c_result_SSc_2_1_1.csv`
- `/data/taobo.hu/projects/sfplot-manuscript/t_by_c_SSc_2_1_1_Landmark`

Confidence: Direct for landmark/segment-free output tables; strong inference for final Excel table assembly.

## Additional related local locations

These paths contain output copies or partial notebooks that are relevant but should not be treated as canonical source:

- `C:\Users\taobo.hu\Downloads\synthetic_nested_benchmark.ipynb`
- `C:\Users\taobo.hu\Downloads\Zscore_nested.pdf`
- `C:\Users\taobo.hu\Downloads\t_and_c_result_Xenium_Prime_Mouse_Pup_FFPE_outs.csv`
- `C:\Users\taobo.hu\Downloads\t_by_c_SSc_2_1_1_raw`
- `C:\Users\taobo.hu\Downloads\t_by_c_SSc_2_1_1_Landmark`
- `C:\Users\taobo.hu\Downloads\wilcox_results_TNBC.csv`
- WPS cache output PDFs including `COSTE_all_nested.pdf`, `StructureMap_of_Mouse_Pup.pdf`, and `twelve_nested_patterns.pdf`

## Not found or unresolved

- No relevant PDC code was found in the searched PDC home tree.
- The Science manuscript directory does not contain the original computational notebooks for most panels.
- Final multi-panel figure composition appears to have been assembled from generated PDFs/PNGs rather than from one reproducible script per final figure.
- Some old/local notebook filenames do not match the current manuscript figure numbering. Example: `D:\GitHub\sfplot\sfplot-manuscript\figures\fig. 1.ipynb` is mouse pup Figure 2 content, not synthetic Figure 1 content.
- `D:\GitHub\sfplot\sfplot-manuscript\淋巴结.ipynb` is present but did not contain relevant source cells in the parsed local copy; the direct lymph node code is the Y-drive R script.

## Recommended source-of-truth set

For future reproduction, treat these as the canonical code set:

- Core method: `D:\GitHub\sfplot\src\sfplot`
- Synthetic benchmark: A100 `/data/taobo.hu/projects/mouse_pup/benchmarking_complete_english.ipynb`, `benchmarking_complete_modular.ipynb`, `benchmarking_complete_nested.ipynb`
- Mouse pup and 5K Xenium panels: `D:\GitHub\sfplot\sfplot-manuscript\figures`, `D:\GitHub\sfplot\sfplot-manuscript\t_and_c`, and `Y:\long\10X_datasets\Xenium\Xenium_5K\t_by_c_result`
- Lung fibrosis COSTE: `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\notebook` and `Y:\long\publication_datasets\Vannan_2023_Lung_Fibrosis\Rcode\sfplotR`
- DST-GNN: A100 `/data/taobo.hu/tmp/dst_gnn_release_repo` plus output folders under `/data/taobo.hu/tmp/dst_gnn_public_outputs*`
- SSc/pleura transcript-by-cell: `D:\GitHub\sfplot\sfplot-manuscript\suppl. fig. 8 new.ipynb`, `D:\GitHub\sfplot\sfplot-manuscript\fig. 3 new.ipynb`, and `/data/taobo.hu/projects/ILD`
- TNBC: `Y:\long\publication_datasets\Ali_TNBC_2023\NTPublic\sfplot\sfplot_TNBC.ipynb` and its local/A100 mirrors
