
|Build Status| |Documentation|

HyDe: hybridization detection using phylogenetic invariants
-----------------------------------------------------------

Installation
------------

Requirements:
~~~~~~~~~~~~~

-  Python 2.7
-  Numpy
-  Pandas
-  C++ compiler (g++ >= v4.8 preferred)

HyDe is a Python package that wraps a primary C++ program for running a
hybridization detection analysis. We have provided a ``Makefile`` that
will compile the C++ executable (``hyde_cpp``) and will then install the
HyDe Python package using ``pip``. To ensure that the necessary
dependencies are available, we suggest using a Python distribution such
as `Miniconda <https://conda.io/miniconda.html>`__.

.. code:: bash

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

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest
   
