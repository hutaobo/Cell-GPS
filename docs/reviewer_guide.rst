==============
Reviewer Guide
==============

This repository contains the code used to develop the Cell-GPS
workflow described in the manuscript.

Package naming
--------------

* Python distribution: ``Cell-GPS``
* Python import package: ``cellgps``
* Python GitHub repository: https://github.com/hutaobo/cellgps
* R package/repository: ``cellgpsr`` at https://github.com/hutaobo/cellgpsr
* Windows executable: ``cellgps.exe`` at https://zenodo.org/records/17859173

Suggested reading order
-----------------------

* ``README.md``
* ``src/sfplot/analysis/searcher_findee_score.py``
* ``src/sfplot/preprocessing/data_processing.py``
* ``src/sfplot/tbc_analysis.py``

Repository map
--------------

* ``src/cellgps/``: recommended Python import namespace
* ``src/sfplot/``: legacy compatibility namespace and current implementation modules
* ``docs/``: package documentation
* ``Cell-GPS manuscript code/``: curated manuscript figure and table notebooks
* ``examples/``: compact package usage examples and small example data files
* ``packaging/``: historical packaging recipes and Windows executable build assets

Minimal install
---------------

.. code-block:: console

   $ git clone https://github.com/hutaobo/cellgps.git
   $ cd cellgps
   $ pip install -e .

Manuscript validation scope
---------------------------

The manuscript-validated scope is the COSTE/SSS workflow and the figure/table
analyses mapped in ``Cell-GPS manuscript code/`` and
``docs/cellgps_science_manuscript_code_inventory.md``. Ligand-receptor
topology, pathway topology, Visium helpers, GUI entry points and other
convenience APIs are included for reuse and development, but should be treated
as optional or exploratory unless a manuscript notebook or documentation page
explicitly maps them to a reported analysis.

Input contract
--------------

For the DataFrame workflow, the minimal required columns are ``x``, ``y``, and
``celltype``.

Reproducibility note
--------------------

Raw datasets are not bundled in this repository because of size and
data-distribution constraints. The repository is intended to expose the code,
analysis structure, and plotting workflows used in the study.
