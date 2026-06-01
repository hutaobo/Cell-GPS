========
Cell-GPS
========

``Cell-GPS`` is a Python package for spatial topology analysis in spatial omics
data. Cell-GPS is used to compute cophenetic distance-based structure maps,
analyze cell-cell and transcript-cell relationships, and visualize multiscale
tissue organization.

Package names
-------------

* Python distribution: ``Cell-GPS``
* Python import package: ``cellgps``
* R package/repository: ``cellgpsr``
* Windows executable: ``cellgps.exe``
* Windows release: https://zenodo.org/records/19482685

Key features
------------

* Compute cophenetic distance matrices from ``AnnData`` objects or coordinate tables.
* Generate StructureMap heatmaps and circular dendrograms.
* Load and preprocess 10x Xenium outputs.
* Run transcript-by-cell analysis for large spatial datasets.
* Support memory-optimized workflows for larger coordinate tables.

Installation
------------

Install from GitHub:

.. code-block:: bash

   pip install git+https://github.com/hutaobo/Cell-GPS.git

For local inspection:

.. code-block:: bash

   git clone https://github.com/hutaobo/Cell-GPS.git
   cd cellgps
   pip install -e .

Reviewer note
-------------

The recommended import namespace is ``src/cellgps/``. Legacy compatibility modules
remain in ``src/sfplot/``. Curated manuscript figure and supplementary table
notebooks are stored in ``Cell-GPS manuscript code/``. See
``docs/project/REVIEWER_GUIDE.md`` for a short walkthrough of the repository.
