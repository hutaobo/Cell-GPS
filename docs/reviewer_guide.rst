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
* ``sfplot-manuscript/``: manuscript-facing notebooks and outputs
* ``benchmarking/``: benchmarking material
* ``docs/``: package documentation

Minimal install
---------------

.. code-block:: console

   $ git clone https://github.com/hutaobo/cellgps.git
   $ cd cellgps
   $ pip install -e .

Input contract
--------------

For the DataFrame workflow, the minimal required columns are ``x``, ``y``, and
``celltype``.

Reproducibility note
--------------------

Raw datasets are not bundled in this repository because of size and
data-distribution constraints. The repository is intended to expose the code,
analysis structure, and plotting workflows used in the study.
