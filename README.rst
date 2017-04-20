
|Build Status|

HyDe: hybridization detection using phylogenetic invariants
-----------------------------------------------------------

Installation
------------

Requirements:
~~~~~~~~~~~~~

-  Numpy
-  Pandas

HyDe is a Python package that wraps a primary C++ program for running a
hybridization detection analysis. We have provided a ``Makefile`` that
will compile the C++ executable (``hyde_cpp``) and will then install the
HyDe Python package using ``pip``. To ensure that the necessary
dependencies are available, we suggest using a Python distribution such
as `Miniconda <https://conda.io/miniconda.html>`__.

.. code:: bash

    # If you don't have numpy or pandas, but do have pip (available with Miniconda)
    pip install numpy
    pip install pandas

    # Clone HyDe repository from GitHub
    git clone https://github.com/pblischak/HyDe.git
    cd HyDe

    # If you don't have OpenMP for multithreading,
    # simply compile using `make`
    make

    # If you do have OpenMP for multithreading,
    # Set the OPENMP flag to 'yes'
    make OPENMP=yes

    # Now install HyDe Python package
    pip install .

    # Test the installation
    make test

Running HyDe from the Command Line
----------------------------------

Type ``run_hyde.py -h`` for options.

.. code:: bash

    run_hyde.py -i <infile> -m <map-file> -o <outgroup> \
                -n <num-ind> -t <num-taxa> -s <num-sites> \
                -j <num-threads> --prefix <prefix>

Running HyDe in Python
----------------------

.. code:: py

    import hyde.main as hd

    # Run a hyde analysis
    hd.run_hyde(...)

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

