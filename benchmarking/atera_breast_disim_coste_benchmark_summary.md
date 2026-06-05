# Atera WTA Breast di-sim/COSTE Benchmark

Run dates: 2026-06-04 to 2026-06-05

Remote output root:
`/data/taobo.hu/atera_breast_disim_coste_benchmark`

Input:
`/data/taobo.hu/SpatialPerturb/inputs/xenium_wta_breast`

Dataset facts recorded from Xenium metadata:

| field | value |
| --- | ---: |
| sample | WTA Preview, WTA Breast Cancer (FFPE) |
| chemistry | Atera v1 |
| cells | 170057 |
| predesigned gene targets | 18028 |
| median genes per cell | 1543 |
| median transcripts per cell | 2116 |
| high-quality decoded transcripts | 624095990 |

## Methods

The benchmark compares two readouts from the same directed matrix.

- COSTE: cluster the directed Searcher->Findee matrix in row and column directions, then convert each hierarchy to a normalized cophenetic distance matrix.
- PNAS di-sim: transform the same directed distance matrix to an affinity matrix, build the regularized directed graph Laplacian, compute top left/right singular vectors, row-normalize each embedding, and k-means cluster row and column roles.

For gene-level runs, the script uses the full 18028-gene panel as the feature universe, selects the top expressed genes for the benchmark size, bins 170057 cells into a 96 x 96 spatial grid, and constructs a directed gene->gene COSTE-style matrix as the expression-weighted mean distance from each source-gene spatial distribution to each target-gene spatial centroid.

## Results

| run | matrix | matrix build sec | COSTE sec | di-sim sec | peak build GB | peak COSTE GB | peak di-sim GB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| celltype_full | 20 x 20 | 3.92 | 0.15 | 0.33 | 0.23 | 0.19 | 0.19 |
| gene_top1024_grid96 | 1024 x 1024 | 31.72 | 0.86 | 3.40 | 5.30 | 0.90 | 0.92 |
| gene_top2048_grid96 | 2048 x 2048 | 41.82 | 22.25 | 4.08 | 5.82 | 1.04 | 1.01 |
| gene_top4096_grid96 | 4096 x 4096 | 44.39 | 189.03 | 11.48 | 6.19 | 1.63 | 1.70 |

| run | COSTE row corr | COSTE col corr | row COSTE-vs-di-sim Pearson | col COSTE-vs-di-sim Pearson | di-sim row/col ARI |
| --- | ---: | ---: | ---: | ---: | ---: |
| celltype_full | 0.906 | 0.976 | 0.048 | 0.258 | -0.040 |
| gene_top1024_grid96 | 0.808 | 0.963 | 0.443 | 0.355 | 0.319 |
| gene_top2048_grid96 | 0.826 | 0.940 | 0.440 | 0.391 | 0.271 |
| gene_top4096_grid96 | 0.825 | 0.936 | 0.371 | 0.317 | 0.234 |

## Interpretation

The two methods share the same mathematical object: an asymmetric matrix with separate row/source and column/target roles. di-sim factorizes a degree-regularized version of that matrix into low-dimensional row and column embeddings; COSTE turns row and column profiles into hierarchical cophenetic metric spaces. The benchmark reflects that relationship: the gene-level di-sim embedding distances correlate with COSTE cophenetic distances, but they are not interchangeable. COSTE preserves a hierarchy over the full directed profiles, while di-sim compresses those profiles into a small number of singular-vector dimensions and flat k-means labels.

Runtime separates cleanly by step. Matrix construction is dominated by reading and aggregating the full Xenium HDF5 matrix. Once the directed matrix exists, di-sim is much faster than exact COSTE cophenetic clustering at 2048 and 4096 genes. COSTE's cost grows sharply because exact row/column pairwise distances and hierarchical linkage scale with the square of item count and profile length; di-sim's truncated SVD stays comparatively cheap for the tested ranks.

## Full di-sim Then Selected COSTE

Output:
`/data/taobo.hu/atera_breast_disim_coste_benchmark/full_disim_18028_select1024_min5000_grid96_k16`

This run used all 18028 genes for di-sim, then selected 1024 genes for COSTE from the di-sim importance ranking. The selected-COSTE eligibility filter required total Xenium counts >= 5000 to avoid choosing unstable low-count genes.

Selection score:
`0.7 * di-sim row/column asymmetry percentile + 0.3 * spectral leverage percentile`

| step | items | seconds | peak GB |
| --- | ---: | ---: | ---: |
| matrix build | 18028 | 52.60 | 6.98 |
| full di-sim | 18028 | 60.20 | 14.02 |
| selected COSTE cophenetic | 1024 | 0.86 | 2.13 |

| metric | value |
| --- | ---: |
| eligible genes after count filter | 12339 |
| selected genes for COSTE | 1024 |
| full di-sim row/column ARI | 0.188 |
| full di-sim embedding asymmetry RMSE | 0.321 |
| selected COSTE row cophenetic corr | 0.774 |
| selected COSTE col cophenetic corr | 0.911 |
| selected row COSTE-vs-di-sim Pearson | 0.265 |
| selected col COSTE-vs-di-sim Pearson | 0.179 |

Top selected genes by the filtered di-sim importance score included `SLC30A8`, `CLIC6`, `HSPB8`, `GRIA2`, `NIBAN1`, `PIP`, `SERPINA6`, `MYBPC1`, `TAT`, `KCNQ3`, `MSMB`, `BMPER`, and `SERPINA1`.

Interpretation: full di-sim is practical for all 18028 genes on A100 when the full directed matrix is stored as float32 `.npy` rather than CSV. The di-sim-selected COSTE subset gives a focused hierarchical readout over genes with strong directional source/target behavior. The low row/column ARI and moderate asymmetry RMSE indicate that many genes have different source-like and target-like roles, which is exactly the signal di-sim is designed to expose before COSTE is used for detailed cophenetic structure.

## COSTE-Seeded HistoSeg Domains

Output:
`/data/taobo.hu/atera_breast_coste1024_histoseg/coste1024_modules16_domains12`

This run used the 1024 genes selected from full di-sim, clustered them into 16 COSTE gene modules using the combined row/column cophenetic matrix, projected module expression back to all 170057 cells, smoothed module scores over spatial nearest neighbors, and clustered cells into 12 HistoSeg-style spatial domains.

| field | value |
| --- | ---: |
| selected COSTE genes | 1024 |
| COSTE gene modules | 16 |
| HistoSeg-style domains | 12 |
| cells assigned | 170057 |
| spatial smoothing neighbors | 16 |
| runtime seconds | 23.35 |

Main outputs:

- `histoseg_domains.parquet`: cell-level domain calls with `histoseg_structure_id`, `histoseg_structure_name`, and approximate boundary distance.
- `histoseg_domain_summary.csv`: domain sizes, dominant cell types, top COSTE module, and representative genes.
- `coste_gene_modules.csv`: 1024 selected genes assigned to COSTE modules.
- `coste_histoseg_domains.png/svg`: spatial domain preview.

Selected domain examples:

| domain | cells | top module genes | dominant cell type |
| --- | ---: | --- | --- |
| COSTE-HistoSeg-02 | 29984 | SERPINA6/KCNQ3/MSMB/BMPER/SERPINA1 | 11q13 Invasive Tumor Cells |
| COSTE-HistoSeg-04 | 9080 | HSPB8/PIP | Luminal-like Amorphous DCIS Cells |
| COSTE-HistoSeg-08 | 35333 | CCL22/SAA2/MS4A1/LTB/IGLC7 | CAFs, DCIS Associated |
| COSTE-HistoSeg-10 | 2576 | SLC30A8 | 11q13 Invasive Tumor Cells |

Interpretation: this is a COSTE-driven molecular HistoSeg, not an H&E-image HistoSeg run. It uses the COSTE-selected genes as spatial molecular landmarks, so the domains should be interpreted as expression/topology-informed tissue regions. The result separates broad stromal/immune-rich areas, invasive tumor-rich regions, and smaller luminal/apocrine-like islands while preserving a cell-level output format compatible with downstream histoseg-style domain summaries.

## Cell-Free Transcript-Grid Prototype and di-sim/COSTE Input

Output:
`/data/taobo.hu/atera_breast_cellfree_transcript_pipeline/full_transcript_grid96_qv20_select1024_geneonly`

Input transcript zarr:
`/data/taobo.hu/pyxenium_lr_benchmark_2026-04/data/source_cache/breast/WTA_Preview_FFPE_Breast_Cancer_outs/spatialdata.zarr/points/transcripts`

This run does not use cells, cell segmentation, or the Xenium `cell_feature_matrix`. It scans transcript points directly, keeps valid transcripts with QV >= 20, removes control/codeword/intergenic identities, aggregates the remaining transcript points into a 96 x 96 spatial grid, and then runs the same analysis logic:

1. full 18028-gene di-sim on the transcript-grid directed matrix;
2. select 1024 genes by filtered di-sim importance score;
3. run COSTE/Co-phonetic on the selected 1024 genes;
4. cluster COSTE genes into 16 modules;
5. project module scores back to occupied transcript-grid bins and cluster 12 HistoSeg-style spatial domains.

| field | value |
| --- | ---: |
| transcript rows scanned | 740442119 |
| valid QV >= 20 transcripts | 624480342 |
| gene-expression transcripts after control/intergenic filtering | 624095990 |
| non-gene identities seen | 9050 |
| genes retained | 18028 |
| occupied 96 x 96 bins | 9072 |
| selected COSTE genes | 1024 |
| COSTE gene modules | 16 |
| HistoSeg-style domains | 12 |

| step | items | seconds | peak GB |
| --- | ---: | ---: | ---: |
| transcript metadata scan | 18028 | 97.87 | 0.98 |
| transcript grid aggregation | 18028 | 62.27 | 1.34 |
| cell-free directed matrix | 18028 | 4.03 | 3.63 |
| full di-sim | 18028 | 28.28 | 7.27 |
| selected COSTE cophenetic | 1024 | 0.82 | 3.69 |
| cell-free HistoSeg | 9072 bins | 2.74 | 3.74 |

| metric | value |
| --- | ---: |
| selected COSTE row cophenetic corr | 0.739 |
| selected COSTE col cophenetic corr | 0.901 |
| full di-sim row/column ARI | 0.176 |
| full di-sim embedding asymmetry RMSE | 0.319 |
| eligible genes after count filter | 12925 |

Selected cell-free domain examples:

| domain | bins | transcripts | top module genes |
| --- | ---: | ---: | --- |
| CellFree-HistoSeg-01 | 828 | 114511008 | GRIA2/KCNJ3/ELOVL2/SEZ6L2/IGHG3 |
| CellFree-HistoSeg-03 | 1134 | 95508360 | CLIC6/NIBAN1/CACNG4/PRKACB/S100G |
| CellFree-HistoSeg-04 | 819 | 134560208 | CXCL13 |
| CellFree-HistoSeg-10 | 143 | 20046884 | SLC30A8 |
| CellFree-HistoSeg-12 | 1137 | 115088864 | HBB/SERPINA6/MSMB/KCNQ3/SERPINA1 |

## Cell-Free vs Cell-Based Performance

| workflow | comparable step | items | seconds | peak GB |
| --- | --- | ---: | ---: | ---: |
| cell-free transcript | scan transcript metadata | 18028 genes | 97.87 | 0.98 |
| cell-free transcript | aggregate transcript grid | 18028 genes | 62.27 | 1.34 |
| cell-free transcript | build directed matrix | 18028 genes | 4.03 | 3.63 |
| cell-free transcript | full di-sim | 18028 genes | 28.28 | 7.27 |
| cell-free transcript | selected COSTE | 1024 genes | 0.82 | 3.69 |
| cell-free grid prototype | HistoSeg domains, superseded | 9072 bins | 2.74 | 3.74 |
| cell-free transcript-point | collect HistoSeg training points | 750000 points | 36.16 | 0.39 |
| cell-free transcript-point | fit point HistoSeg | 750000 points | 24.69 | 0.76 |
| cell-free transcript-point | assign selected transcript points | 60906532 points | 356.29 | 0.94 |
| cell-based | build directed matrix from cells | 18028 genes | 52.60 | 6.98 |
| cell-based | full di-sim | 18028 genes | 60.20 | 14.02 |
| cell-based | selected COSTE | 1024 genes | 0.86 | 2.13 |
| cell-based | HistoSeg domains | 170057 cells | 23.35 | n/a |

Interpretation: the cell-free path spends most of its time reading transcript points and then assigning point domains. The old 9072-bin HistoSeg timing is fast but is only a coarse prototype. The corrected transcript-point HistoSeg assigns 60906532 selected-gene transcript points in 356.29 seconds after a 60.85-second sample-and-fit stage, while keeping peak RSS below 1 GB. That is slower than the 170057-cell HistoSeg run because it assigns about 358x more spatial objects, but it matches the intended cell-free unit: transcript points from the genes already selected by di-sim/COSTE.

The cell-free PNG preview is a valid spatial map, but it is a grid-domain map rather than a cell-domain map. It preserves the tissue outline and broad domain structure, with more block-like regions because each call belongs to a 96 x 96 grid bin instead of an individual cell. For that reason, the HistoSeg part of this run should be treated as a coarse prototype only. The corrected cell-free HistoSeg result below uses transcript points directly.

## Corrected Transcript-Point HistoSeg

Output:
`/data/taobo.hu/atera_breast_transcript_point_histoseg/coste1024_point_knn_train750k_assign_all`

Script:
`benchmarking/atera_breast_transcript_point_histoseg.py`

This corrected run does not build a spatial grid. It takes the 1024 genes selected by full di-sim and organized by COSTE modules, then works directly on the transcript points from those selected genes. The HistoSeg-style feature for a point is the local composition of COSTE modules in its transcript-point neighborhood, plus a weak spatial coordinate term. The model trains on a 750000-point reservoir sample and then streams through the transcript table again to assign every selected-gene transcript point.

| field | value |
| --- | ---: |
| uses spatial grid | false |
| HistoSeg unit | transcript point |
| selected COSTE genes | 1024 |
| selected-gene transcript points | 60906532 |
| training transcript points | 750000 |
| preview transcript points | 500000 |
| COSTE gene modules | 16 |
| HistoSeg-style domains | 12 |
| local neighborhood k | 128 |
| assignment neighbors | 15 |
| written parquet parts | 371 |

| step | items | seconds | peak GB |
| --- | ---: | ---: | ---: |
| collect transcript-point training sample | 750000 | 36.16 | 0.39 |
| fit transcript-point HistoSeg | 750000 | 24.69 | 0.76 |
| assign selected transcript points | 60906532 | 356.29 | 0.94 |

Selected transcript-point domain examples:

| domain | transcript points | fraction | top module genes |
| --- | ---: | ---: | --- |
| TranscriptPoint-HistoSeg-01 | 22106213 | 0.363 | HBB/SERPINA6/MSMB/KCNQ3/SERPINA1 |
| TranscriptPoint-HistoSeg-02 | 5400966 | 0.089 | RERGL/BMPER/IGF2/DPT/CCL22 |
| TranscriptPoint-HistoSeg-06 | 5965642 | 0.098 | RERGL/BMPER/IGF2/DPT/CCL22 |
| TranscriptPoint-HistoSeg-09 | 5203948 | 0.085 | HBB/SERPINA6/MSMB/KCNQ3/SERPINA1 |
| TranscriptPoint-HistoSeg-12 | 10367916 | 0.170 | HBB/SERPINA6/MSMB/KCNQ3/SERPINA1 |

Main outputs:

- `transcript_point_domain_chunks/part_*.parquet`: point-level domain calls for all 60906532 selected-gene transcript points.
- `transcript_point_histoseg_domain_summary.csv`: domain sizes and dominant COSTE module summaries.
- `transcript_point_histoseg_training_domains.parquet`: 750000 training transcript points with fitted domain labels.
- `transcript_point_histoseg_preview_points.parquet/csv`: 500000 point preview used for plotting.
- `transcript_point_histoseg_domains.png/svg`: no-grid transcript-point spatial preview.

Interpretation: this is the corrected cell-free HistoSeg formulation. It does not try to segment all 624M transcripts, because di-sim has already selected the informative genes. Instead, it segments the 60.9M transcript points belonging to the 1024 selected genes. Runtime is dominated by streaming assignment and parquet writing, not memory; peak RSS stayed below 1 GB. Compared with the grid prototype, the point map is no longer blocky and preserves transcript-level molecular structure, but it is naturally sparse where the selected genes have little expression.

## Transcript-Point vs Cell-Based COSTE-HistoSeg Comparison

Output:
`/data/taobo.hu/atera_breast_histoseg_comparison/point_vs_cell_coste1024`

Script:
`benchmarking/atera_breast_point_vs_cell_histoseg_comparison.py`

Comparison unit: each selected-gene transcript point from the no-grid point-level HistoSeg run was matched to the nearest cell-based COSTE-HistoSeg cell centroid. The overlap is transcript-weighted and uses no spatial grid.

| field | value |
| --- | ---: |
| transcript points compared | 60906532 |
| cell-based cells used as reference | 170057 |
| point domains | 12 |
| cell-based domains | 12 |
| nearest-cell mean distance | 5.12 um |
| nearest-cell median distance | 4.37 um |
| nearest-cell 95th percentile distance | 10.99 um |
| overlap matching time | 224.96 sec |
| peak RSS | 0.36 GB |

| metric | value |
| --- | ---: |
| adjusted Rand index | 0.261 |
| normalized mutual information | 0.392 |
| point-to-cell purity | 0.450 |
| cell-to-point purity | 0.623 |
| weighted mean best Jaccard | 0.314 |

High-overlap domain examples:

| point domain | best cell-based domain | point-to-cell fraction | Jaccard | pathology interpretation |
| --- | --- | ---: | ---: | --- |
| TranscriptPoint-HistoSeg-10 | COSTE-HistoSeg-10 | 0.903 | 0.176 | invasive/tumor-rich -> invasive/tumor-rich |
| TranscriptPoint-HistoSeg-03 | COSTE-HistoSeg-03 | 0.790 | 0.776 | immune/stroma signature over apocrine region |
| TranscriptPoint-HistoSeg-07 | COSTE-HistoSeg-12 | 0.685 | 0.353 | tumor-rich -> tumor/stroma interface |
| TranscriptPoint-HistoSeg-06 | COSTE-HistoSeg-08 | 0.627 | 0.511 | immune/stroma interface -> immune/stroma interface |
| TranscriptPoint-HistoSeg-01 | COSTE-HistoSeg-02 | 0.402 | 0.338 | invasive/tumor-rich -> invasive/tumor-rich; shared SERPINA6/KCNQ3/MSMB/SERPINA1 signature |

Interpretation: the maps are related but not identical. ARI/NMI are moderate, while point-to-cell purity and the heatmap show clear dominant matches for several domains. This means the point-level HistoSeg is not simply reproducing the cell-level COSTE-HistoSeg; it often splits broad cell-based regions into transcript-level molecular micro-regions. Pathology interpretation is strongest when both overlap and marker logic agree, such as the SERPINA6/MSMB/KCNQ3/SERPINA1 tumor-rich domains and the CCL22/SAA2/MS4A1/LTB/IGLC7 immune/stroma-associated domain. It is weaker or mixed when the selected transcript signature falls inside a different cell-level context, such as RERGL/BMPER/IGF2/DPT/CCL22 points overlapping luminal/DCIS or apocrine cell-domain regions. These labels are molecular/pathology heuristics, not formal H&E diagnoses.

Main outputs:

- `point_vs_cell_overlap_counts.csv`: 12 x 12 transcript-weighted overlap counts.
- `point_vs_cell_overlap_point_fraction.csv`: row-normalized overlap, used to read point-domain purity.
- `point_vs_cell_best_domain_matches.csv`: best cell-domain match, Jaccard, nearest dominant cell type, and marker interpretation for each point domain.
- `point_vs_cell_pathology_interpretation.md`: concise pathology interpretation report.
- `point_vs_cell_overlap_heatmap.png/svg`: overlap heatmap.
- `point_vs_cell_spatial_comparison.png/svg`: side-by-side cell map, point map, and point preview colored by nearest cell-domain label.

## Alignment-Based Nature-Style Pathology Explanation Figures

Output:
`/data/taobo.hu/atera_breast_histoseg_nature_figures/pathology_explanations_aligned`

Script:
`benchmarking/atera_breast_nature_pathology_figures.py`

The H&E overlays now use the official 10x Xenium Explorer alignment file:
`WTA_Preview_FFPE_Breast_Cancer_he_alignment.csv`. The remote copy used by the script is
`/data/taobo.hu/pyxenium_lazyslide_breast_wta_20260507/data/WTA_Preview_FFPE_Breast_Cancer_he_alignment.csv`,
paired with the official H&E OME image and keypoint file. Alignment QC against 12 official keypoints gives mean error 12.24 px, median error 14.42 px, and max error 21.77 px. Using 0.2125 um per fixed image pixel, that is mean error 2.60 um and max error 4.63 um, which is appropriate for transcript/cell overlay interpretation at this scale.

Figure logic:

| panel type | purpose |
| --- | --- |
| story flow | turns each pathology question into a connected five-step point-vs-cell visual argument |
| global H&E overlap | shows tissue-wide location of each point-domain pathology signal |
| local H&E overlap | zooms into the densest representative point-domain/cell-domain overlap region |
| full-data statistics | uses all 60906532 selected-gene transcript points, not the visual preview |
| marker/interpretation panel | links overlap, COSTE marker genes, and pathology interpretation |

Core pathology explanations:

| explanation | point domains and cell-domain context | statistical support | interpretation |
| --- | --- | --- | --- |
| tumor-rich molecular signal | SERPINA6/MSMB/KCNQ3/SERPINA1-rich point domains overlap invasive/tumor-rich COSTE-HistoSeg contexts | P10->C10 has 0.903 point-to-cell fraction; P07->C12 has 0.685; P01->C02 has 0.402 | the transcript-point map clearly recovers tumor-rich molecular regions and further splits broad tumor cell domains into transcript-level subdomains |
| immune/stroma interface | CCL22/DPT/IGF2/BMPER/RERGL-rich point domains overlap immune/stromal interface contexts | P06->C08 has 0.627 point-to-cell fraction and 0.511 Jaccard; P03->C03 has 0.790 and 0.776 | the point-level map captures an immune/stromal interface signal, rather than a pure tumor-only signal |
| mixed luminal/apocrine microenvironment | RERGL/BMPER/IGF2/DPT/CCL22-like point domains fall inside luminal/DCIS C04 and apocrine C03 cell-domain contexts | P02->C04 has 0.494; P11->C04 has 0.541; P03->C03 has 0.790 | these regions are best interpreted as local microenvironment mixing inside epithelial pathological contexts, not as a single pure cell type zone |

Main exported figures:

- `nature_pathology_story_tumor_rich_steps.png/svg/pdf/tiff`: five-step story figure for the tumor-rich question, moving from P10 point-domain global, to C10 cell-domain global, to P10/C10 local overlap, top full-data support, and pathology answer.
- `nature_pathology_story_immune_stroma_steps.png/svg/pdf/tiff`: five-step story figure for the immune/stroma interface question.
- `nature_pathology_story_mixed_luminal_apocrine_steps.png/svg/pdf/tiff`: five-step story figure for the mixed luminal/DCIS and apocrine microenvironment question, with separate local zooms for both local contexts.
- `nature_pathology_triptych_summary.png/svg/pdf/tiff`: one compact three-row summary of the three pathology explanations.
- `nature_pathology_tumor_rich_plate.png/svg/pdf/tiff`: global H&E, local H&E, full-data statistics, and marker interpretation for the tumor-rich signal.
- `nature_pathology_immune_stroma_plate.png/svg/pdf/tiff`: the same four-panel structure for the immune/stroma interface signal.
- `nature_pathology_mixed_luminal_apocrine_plate.png/svg/pdf/tiff`: the same structure for the mixed luminal/DCIS and apocrine microenvironment signal.
- `nature_pathology_source_domain_statistics.csv`: source rows from the full point-vs-cell domain statistics used in the figures.
- `nature_pathology_figure_manifest.json`: input paths, alignment QC, output paths, and ROI coordinates.
- `nature_pathology_figure_contract.json`: figure intent and image-integrity notes.

Important caveat: the H&E overlay panels use a 500000-point preview plus cell centroids so the figure remains readable. The bar plots and interpretation statistics use the full 60906532 selected-gene transcript-point comparison. No spatial grid is used in this pathology figure workflow.

Story-flow design note: the three `nature_pathology_story_*_steps` figures are intended for presentation or manuscript narrative. Each one asks a single pathology question and then follows the same evidence chain: example point-domain global on official-aligned H&E, matched cell-domain global on the same H&E context, densest local point/cell overlap, top full-data point-vs-cell support, and a final pathology answer. This makes the point-domain and cell-domain evidence separable before the zoom panel, while keeping the three findings visually consistent.

Nature-format QA for story-flow figures: the three `nature_pathology_story_*_steps` outputs are exported at exact double-column size, 183.0 x 91.9 mm for PNG/PDF/SVG and 183.0 x 92.0 mm for TIFF. PNG files are 3242 x 1629 px at 450 dpi; TIFF files are 4322 x 2173 px at 600 dpi. Story-flow text is kept in the 5-7 pt range, SVG text remains editable, PDF font embedding uses TrueType text, and the H&E panels are top-anchored to reduce unnecessary whitespace. The right-side pathology-answer panels and bottom explanation boxes were shortened to avoid clipping at final size.
