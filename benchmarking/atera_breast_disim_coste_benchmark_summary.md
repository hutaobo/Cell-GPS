# Atera WTA Breast di-sim/COSTE Benchmark

Run date: 2026-06-04

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
