.. highlight:: shell

============
Installation
============

From PyPI
---------

Install the Python package:

.. code-block:: console

   $ pip install Cell-GPS

From GitHub
-----------

Install directly from the public repository:

.. code-block:: console

   $ pip install git+https://github.com/hutaobo/cellgps.git

For local inspection
--------------------

Clone the repository and install it in editable mode:

.. code-block:: console

   $ git clone https://github.com/hutaobo/cellgps.git
   $ cd cellgps
   $ pip install -e .

Requirements
------------

* Python 3.9 or newer
* Standard scientific Python stack, including ``numpy``, ``pandas``,
  ``matplotlib``, ``seaborn``, ``scanpy``, ``scikit-learn``, and
  ``pyXenium`` for Xenium I/O

Optional extras
---------------

Install Xenium support:

.. code-block:: console

   $ pip install -e .[xenium]

Install Visium support:

.. code-block:: console

   $ pip install -e .[visium]

Related Packages
----------------

* Python repository: https://github.com/hutaobo/cellgps
* R package repository: https://github.com/hutaobo/cellgpsr
* Windows executable release: https://zenodo.org/records/17859173
