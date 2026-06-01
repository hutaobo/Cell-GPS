===================
Manuscript Overview
===================

This page introduces the Cell-GPS/COSTE preprint and explains how the
manuscript code in this repository is organized.

Preprint
========

Cell-GPS/COSTE is described in the bioRxiv preprint:

   Long M, Hu T, Sountoulidis A, Samakovlis C, Nilsson M.
   Cophenetic Spatial Topology Embedding reveals multiscale tissue
   architecture in spatial omics. bioRxiv. 2026.
   doi: https://doi.org/10.64898/2026.05.26.727847

The versioned bioRxiv page is:

   https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1

The preprint presents Cophenetic Spatial Topology Embedding (COSTE), a
spatial topology framework that converts directed nearest-neighbor distance
profiles between cell populations, transcripts, or other spatial entities
into a hierarchical metric space. The resulting StructureMap and Spatial
Separation Score (SSS) summarize how tissue components are arranged across
the full analyzed field, without requiring a user-selected neighborhood
radius.

What the Manuscript Shows
=========================

The manuscript uses synthetic spatial benchmarks, mouse pup Xenium data,
systemic sclerosis and pulmonary fibrosis spatial transcriptomics, and
triple-negative breast cancer data to show how COSTE captures tissue-scale
architecture that can be difficult to summarize with local neighborhood
enrichment alone.

The main figure flow is:

* Figure 1: synthetic modular and nested spatial benchmarks.
* Figure 2: neonatal mouse pup Xenium tissue structures.
* Figure 3: pulmonary fibrosis and TRU remodeling analysis.
* Figure 4: fibrosis progression and DST-GNN analysis.
* Figure 5: segment-free pleural transcript/cell StructureMap analysis.

The supplementary figures document method comparisons, parameter sensitivity,
lymph node analysis, unclustered fibrosis heatmaps, transcript/cell panels,
and TNBC biomarker heatmaps. Supplementary tables cover runtime/memory
benchmarks and SSc transcript-by-cell SSS tables.

Repository Companion Code
=========================

The curated manuscript notebooks are stored in:

   https://github.com/hutaobo/cellgps/tree/main/Cell-GPS%20manuscript%20code

The folder is organized by manuscript result:

* ``main_figures/`` contains one notebook per main figure.
* ``supplementary_figures/`` contains one notebook per supplementary figure.
* ``supplementary_tables/`` contains one notebook per supplementary table.

The notebooks are intentionally output-free. They preserve the original data
paths used during manuscript preparation and add explanatory comments so that
readers can understand the analysis and rerun the parts for which they have
access to the corresponding datasets.

Citation
========

If you use the Python package, Windows executable, R companion package, or
the manuscript figure/table code, please cite the preprint:

.. code-block:: bibtex

   @article{long2026cellgps,
     title = {Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics},
     author = {Long, Mengping and Hu, Taobo and Sountoulidis, Alexandros and Samakovlis, Christos and Nilsson, Mats},
     journal = {bioRxiv},
     year = {2026},
     doi = {10.64898/2026.05.26.727847},
     url = {https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1}
   }
