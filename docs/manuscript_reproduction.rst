===========================
Reproducing Manuscript Code
===========================

This page describes a practical workflow for using the Cell-GPS manuscript
code. The goal is to make the computational source for each figure and table
easy to inspect, rerun where data access permits, and adapt to related spatial
omics datasets.

Quick Start
===========

Install the package and clone the repository:

.. code-block:: console

   pip install Cell-GPS
   git clone https://github.com/hutaobo/cellgps.git
   cd cellgps

For editable local development:

.. code-block:: console

   pip install -e .

The manuscript notebooks are under:

.. code-block:: text

   Cell-GPS manuscript code/
     main_figures/
     supplementary_figures/
     supplementary_tables/

Open the notebook corresponding to the figure or table of interest, inspect
the first markdown cell for provenance, then adapt the data paths if needed.

Data Availability and Paths
===========================

Large raw spatial omics datasets are not bundled with the repository. The
notebooks therefore preserve the original paths used during analysis, such as
``Y:\long\...`` on Windows workstations and ``/data/taobo.hu/...`` or
``/mnt/taobo.hu/...`` on A100/Linux systems.

When rerunning a notebook, use one of these strategies:

* Run on the same workstation or server where the original paths are mounted.
* Replace the original path variables with local copies of the same datasets.
* Use the notebook as a code template and substitute your own Xenium,
  transcript-coordinate, or cell-coordinate data.

The final multi-panel manuscript figures were usually assembled from generated
PDF/PNG panels. The notebooks focus on the computational source for those
panels rather than reproducing Illustrator or PowerPoint assembly exactly.

Recommended Reading Order
=========================

For readers who want to understand the method before rerunning full
manuscript analyses:

1. Read :doc:`manuscript_overview`.
2. Run the simple coordinate-table example in :doc:`usage`.
3. Inspect ``Figure_1_synthetic_benchmark.ipynb`` for the benchmark logic.
4. Inspect the figure-specific notebook listed in :doc:`manuscript_code_index`.
5. Use ``docs/cellgps_science_manuscript_code_inventory.md`` for deeper
   provenance notes when a notebook refers to legacy local or server paths.

Core Package Entry Points
=========================

Most manuscript notebooks use one or more of these public entry points:

* ``compute_cophenetic_distances_from_df`` for coordinate-table inputs.
* ``compute_cophenetic_distances_from_adata`` for AnnData/Xenium workflows.
* ``plot_cophenetic_heatmap`` for StructureMap heatmaps.
* ``transcript_by_cell_analysis`` for transcript-to-cell spatial topology.
* ``compute_cophenetic_distances_from_df_memory_opt`` for large tables.
* ``plot_circular_dendrogram_pycirclize`` for circular hierarchy plots.

Reproducibility Status
======================

The curated notebooks are designed to expose the analysis logic and code
provenance. They are not all guaranteed to execute on a fresh public machine
because some require large datasets, local manuscript outputs, or server-side
intermediate files. This is the expected status by category:

* Synthetic benchmark notebooks are the easiest to rerun when dependencies are
  installed, because they generate benchmark patterns programmatically.
* Mouse pup, lymph node, pulmonary fibrosis, SSc, and TNBC notebooks require
  access to the corresponding public or local processed spatial omics data.
* DST-GNN notebooks require the recovered DST-GNN release workspace and its
  flattened SSS input tables.
* Supplementary table notebooks require the intermediate benchmark or
  transcript-by-cell result tables used for the manuscript.

Validation
==========

For a local documentation check:

.. code-block:: console

   python -m sphinx -b html docs docs/_build/html

For notebook validation without rerunning heavy computations, use
``nbformat`` to confirm that the notebook files are structurally valid.
GitHub rendering also expects notebooks to include stable cell IDs, which
these curated notebooks now include.
