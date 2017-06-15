
|Build Status| |Documentation|

HyDe: hybridization detection using phylogenetic invariants
-----------------------------------------------------------

`Read the Docs <http://hybridization-detection.rtfd.io/>`__
-----------------------------------------------------------

Installation
------------

Requirements:
~~~~~~~~~~~~~

-  Python 2.7
-  Cython
-  Numpy
-  Scipy
-  Pandas
-  Matplotlib
-  Seaborn
-  C++ compiler (g++ >= v4.8 preferred)

HyDe is a Python package that wraps a primary C++ program for running a
hybridization detection analysis. We have provided a ``Makefile`` that
will compile the C++ executable (``hyde_cpp``) and will then install the
HyDe Python package using the ``setup.py`` file. To ensure that the necessary
dependencies are available, we suggest using a Python distribution such
as `Miniconda <https://conda.io/miniconda.html>`__.

.. code:: bash

    # To install dependencies
    pip install cython numpy scipy pandas matplotlib seaborn

    # Clone HyDe repository from GitHub
    git clone https://github.com/pblischak/HyDe.git
    cd HyDe

    # Compile hyde_cpp using `make`
    make

    # Now install HyDe Python package
    python setup.py install

    # Test the installation
    make test

Running HyDe from the Command Line
----------------------------------

Type ``run_hyde.py -h`` for options.

.. code:: bash

    run_hyde.py -i <infile> -m <map-file> -o <outgroup> \
                -n <num-ind> -t <num-taxa> -s <num-sites> \
                --prefix <prefix>

Running HyDe in Python
----------------------

.. code:: py

    import hyde as hd

    # Run a hyde analysis
    hd.run_hyde(...)

    # import a data set for testing particular triples
    data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 250000)

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest
