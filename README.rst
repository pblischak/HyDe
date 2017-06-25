
|Build Status| |Documentation|  |PyPI Badge|  |Gitter|

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
-  C++ compiler (g++ >= v4.8)

HyDe is a software package that detects hybridization in phylogenomic
data sets using phylogenetic invariants. The primary interface for HyDe is a Python
module called ``phyde`` (**P**\ ythonic **Hy**\ bridization **De**\ tection).
``phyde`` provides a suite of tools for performing hypothesis tests on triples of taxa
to detect hybridization. It also has built in functions to wrap calls to the pure C++ version
of HyDe, ``hyde_cpp``. We have provided a ``Makefile`` that
will compile the ``hyde_cpp`` C++ executable and will then install the
``phyde`` Python package using the ``setup.py`` file. To ensure that the necessary
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

    # Now install phyde module
    python setup.py install

    # Test the installation
    make test

The ``phyde`` module is also hosted on the Python Package Index (PyPI), and can be installed directly using
``pip``.

.. code:: py

  # Install from PyPI with pip
  pip install phyde

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

    import phyde as hd

    # Run a hyde analysis
    res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=100)

    # import a data set for testing particular triples
    data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest

.. |PyPI Badge| image:: https://badge.fury.io/py/phyde.svg
   :target: https://pypi.python.org/pypi/phyde

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-HyDe/Lobby
