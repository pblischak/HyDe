
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
-  Python Modules:

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

.. code:: bash

  # Install from PyPI with pip
  pip install phyde

Example HyDe analysis in Python
-------------------------------

.. code:: py

  # Import the modules that we'll need
  from __future__ import print_function
  import phyde as hd

  # Run a full HyDe analysis to see if there is any evidence for hybridization.
  # This can be done using the run_hyde() function.
  # If you're using ipython, type `hd.run_hyde?` to see more information.
  res_all = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)

  # res_all is an object of class HydeResult. The main results are stored in this object
  # as a dictionary that is called res. We can filter the results using dictionary comprehensions.
  # We will filter for significant results with hybridization parameter estimates (gamma)
  # between 0 and 1.
  filter_res = {k:v for k,v in res_all.res.items() if v['Pvalue'] < 0.025 and v['Gamma'] > 0.0 and v['Gamma'] < 1.0}

  # We can now print these results to a file.
  outfile = open("hyde-filtered-out.txt", 'wa')
  for k,v in filtered_res.items():
    print(k[0], '\t', k[1], '\t', k[2], '\t', sep='', end='', file=outfile)
    for ki,vi in v.items():
      print(ki, ":", vi, '\t', sep='', end='', file=outfile)
    print('\n', end='', file=outfile)
  outfile.close()

  # Next we'll test each individual in each hybrid population.
  # First we have to read the data into Python using the HydeData class
  dat = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)

  # Then we'll run the tests using the text_individuals() function.
  # We'll use a dictionary comprehension again. The significat triples
  # are stored at tuples (p1, hyb, p2).
  ind_tests = {k: dat.test_individuals(k[0], k[1], k[2]) for k in filter_res.keys()}
  print(ind_tests)

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest

.. |PyPI Badge| image:: https://badge.fury.io/py/phyde.svg
   :target: https://pypi.python.org/pypi/phyde

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-HyDe/Lobby
