==============
Reviewer Guide
==============

This repository contains the code used to develop the Cell-GPS
workflow described in the manuscript.

Package naming
--------------

* Python distribution: ``Cell-GPS``
* Python import package: ``cellgps``
* Python GitHub repository: https://github.com/hutaobo/Cell-GPS
* R package/repository: ``cellgpsr`` at https://github.com/hutaobo/cellgpsr
* Windows executable: ``cellgps.exe`` at latest open v2 DOI https://doi.org/10.5281/zenodo.19482685; version-series route https://zenodo.org/records/17859173

Suggested reading order
-----------------------

* ``README.md``
* ``src/cellgps/analysis/searcher_findee_score.py``
* ``src/cellgps/preprocessing/data_processing.py``
* ``src/cellgps/analysis/tbc_analysis.py``

Repository map
--------------

* ``src/cellgps/``: Python package and implementation modules
* ``docs/``: package documentation
* ``Cell-GPS manuscript/``: curated manuscript figure and table notebooks
* ``examples/``: compact package usage examples and small example data files
* ``packaging/``: historical packaging recipes and Windows executable build assets

Minimal install
---------------

.. code-block:: console

   $ git clone https://github.com/hutaobo/Cell-GPS.git
   $ cd cellgps
   $ pip install -e .

Manuscript validation scope
---------------------------

The manuscript-validated scope is the COSTE/SSS workflow and the figure/table
analyses mapped in ``Cell-GPS manuscript/`` and
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
