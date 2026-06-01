=====================
Figure and Table Code
=====================

This page maps the bioRxiv manuscript figures and supplementary tables to the
curated notebooks in the GitHub repository. Each notebook is intentionally
output-free and contains comments describing the original data paths,
provenance, and expected generated files.

Main Figures
============

.. list-table::
   :header-rows: 1
   :widths: 12 28 35 25

   * - Result
     - Topic
     - Curated notebook
     - Reproduction notes
   * - Figure 1
     - Synthetic modular, nested, and random spatial benchmarks; COSTE,
       Squidpy, Giotto, and ANE comparison.
     - `Figure_1_synthetic_benchmark.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/main_figures/Figure_1_synthetic_benchmark.ipynb>`__
     - Canonical source came from A100 benchmark notebooks. Synthetic data
       generation is included, but some final panels were assembled from PDFs.
   * - Figure 2
     - Neonatal mouse pup Xenium StructureMap and hierarchical structures.
     - `Figure_2_mouse_pup_xenium.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/main_figures/Figure_2_mouse_pup_xenium.ipynb>`__
     - Requires mouse pup Xenium outputs and t-by-c result tables.
   * - Figure 3
     - Lung fibrosis COSTE analysis, SSS, regional TRS, and TRU remodeling.
     - `Figure_3_lung_fibrosis_TRU.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/main_figures/Figure_3_lung_fibrosis_TRU.ipynb>`__
     - Requires Vannan lung fibrosis processed data and regional annotations.
   * - Figure 4
     - Fibrosis progression and DST-GNN modeling from SSS matrices.
     - `Figure_4_lung_fibrosis_DST_GNN.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/main_figures/Figure_4_lung_fibrosis_DST_GNN.ipynb>`__
     - Requires flattened SSS tables and recovered DST-GNN workspace.
   * - Figure 5
     - Segment-free SSc pleura transcript/cell StructureMap, circular
       hierarchy, and selected gene panels.
     - `Figure_5_SSc_pleura_segment_free.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/main_figures/Figure_5_SSc_pleura_segment_free.ipynb>`__
     - Requires SSc pleura Xenium/transcript data and t-by-c outputs.

Supplementary Figures
=====================

.. list-table::
   :header-rows: 1
   :widths: 14 31 35 20

   * - Result
     - Topic
     - Curated notebook
     - Notes
   * - Supplementary Figure 1
     - Synthetic modular benchmark.
     - `Supp_Fig_01_synthetic_modular.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_01_synthetic_modular.ipynb>`__
     - A100 modular benchmark source.
   * - Supplementary Figure 2
     - Synthetic nested spatial patterns.
     - `Supp_Fig_02_synthetic_nested_patterns.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_02_synthetic_nested_patterns.ipynb>`__
     - A100 nested benchmark source.
   * - Supplementary Figure 3
     - COSTE nested-pattern heatmaps.
     - `Supp_Fig_03_COSTE_nested_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_03_COSTE_nested_heatmaps.ipynb>`__
     - Uses generated COSTE benchmark outputs.
   * - Supplementary Figure 4
     - Squidpy nested-pattern heatmaps.
     - `Supp_Fig_04_Squidpy_nested_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_04_Squidpy_nested_heatmaps.ipynb>`__
     - Uses generated Squidpy benchmark outputs.
   * - Supplementary Figure 5
     - Giotto nested-pattern heatmaps.
     - `Supp_Fig_05_Giotto_nested_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_05_Giotto_nested_heatmaps.ipynb>`__
     - Uses Giotto benchmark source.
   * - Supplementary Figure 6
     - ANE nested-pattern heatmaps.
     - `Supp_Fig_06_ANE_nested_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_06_ANE_nested_heatmaps.ipynb>`__
     - Uses analytical neighborhood enrichment source.
   * - Supplementary Figure 7
     - Mouse pup hierarchical structures.
     - `Supp_Fig_07_mouse_pup_structures.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_07_mouse_pup_structures.ipynb>`__
     - Requires mouse pup t-by-c and spatial outputs.
   * - Supplementary Figure 8
     - Mouse pup method comparison.
     - `Supp_Fig_08_mouse_pup_method_comparison.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_08_mouse_pup_method_comparison.ipynb>`__
     - Uses runtime/method benchmark code.
   * - Supplementary Figure 9
     - Squidpy parameter sensitivity.
     - `Supp_Fig_09_squidpy_parameter_sensitivity.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_09_squidpy_parameter_sensitivity.ipynb>`__
     - Compares neighbor/radius settings.
   * - Supplementary Figure 10
     - Human lymph node StructureMap and spatial map.
     - `Supp_Fig_10_lymph_node_structuremap.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_10_lymph_node_structuremap.ipynb>`__
     - Direct source is a Y-drive R script plus Xenium lymph node outputs.
   * - Supplementary Figure 11
     - Lung fibrosis unclustered SSS heatmaps.
     - `Supp_Fig_11_lung_fibrosis_unclustered_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_11_lung_fibrosis_unclustered_heatmaps.ipynb>`__
     - Requires lung fibrosis SSS matrices.
   * - Supplementary Figure 12
     - SSc pleura transcript/cell panels.
     - `Supp_Fig_12_SSc_pleura_transcript_cell.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_12_SSc_pleura_transcript_cell.ipynb>`__
     - Uses SSc transcript and cell-level outputs.
   * - Supplementary Figure 13
     - TNBC spatial biomarker statistics.
     - `Supp_Fig_13_TNBC_spatial_biomarkers.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_13_TNBC_spatial_biomarkers.ipynb>`__
     - Requires Ali TNBC clinical and spatial data.
   * - Supplementary Figure 14
     - TNBC subgroup heatmaps.
     - `Supp_Fig_14_TNBC_subgroup_heatmaps.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_figures/Supp_Fig_14_TNBC_subgroup_heatmaps.ipynb>`__
     - Per-patient and subgroup Searcher/Findee heatmaps.

Supplementary Tables
====================

.. list-table::
   :header-rows: 1
   :widths: 16 34 34 16

   * - Result
     - Topic
     - Curated notebook
     - Notes
   * - Supplementary Table 1
     - Runtime and memory benchmarking for COSTE, Squidpy, Giotto, and ANE.
     - `Supp_Table_1_runtime_memory.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_tables/Supp_Table_1_runtime_memory.ipynb>`__
     - Uses A100 mouse pup benchmark scripts and performance counters.
   * - Supplementary Table 2
     - SSc transcript-by-cell SSS.
     - `Supp_Table_2_SSc_transcript_by_cell.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_tables/Supp_Table_2_SSc_transcript_by_cell.ipynb>`__
     - Requires SSc transcript-by-cell result tables.
   * - Supplementary Table 3
     - SSc landmark transcript SSS.
     - `Supp_Table_3_SSc_landmark_transcript_SSS.ipynb <https://github.com/hutaobo/cellgps/blob/main/Cell-GPS%20manuscript%20code/supplementary_tables/Supp_Table_3_SSc_landmark_transcript_SSS.ipynb>`__
     - Requires landmark/segment-free output tables.

Detailed Provenance
===================

For a longer inventory of original local, Y-drive, and A100 locations, see
``docs/cellgps_science_manuscript_code_inventory.md`` in the GitHub
repository:

   https://github.com/hutaobo/cellgps/blob/main/docs/cellgps_science_manuscript_code_inventory.md
