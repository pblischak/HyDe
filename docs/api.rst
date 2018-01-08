.. include:: links.rst

.. _API:

API Reference
=============

**This page is still in the works. Sorry!**

**Submodule**: ``core``
-----------------------

The ``core`` submodule provides the main Python interface for hybridization detection analyses
and data exploration through the use of implemented classes.

**File**: ``data.pyx``
^^^^^^^^^^^^^^^^^^^^^^

The ``data.pyx`` source file is written in Cython, a superset of the Python language
that allows the use of C/C++ code to speed up computationally intensive tasks
(and many other things). This file is automatically translated into a C++ source
file and is compiled into a shared library that can be called from Python.

.. autoclass:: phyde.HydeData
  :members:

**File**: ``result.py``
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: phyde.HydeResult
  :members:

  **Constructor**: ``HydeResult(infile)``

    Example:

    .. code:: py

      import phyde as hd
      res = hd.HydeResult("hyde-out.txt")

**File**: ``bootstrap.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: phyde.Bootstrap
  :members:
  :special-members: __call__

  **Constructor**: ``Bootstrap(bootfile)``

    Read in bootstrap replicates from a HyDe analysis.

    Example:

    .. code:: py

      import phyde as hd
      boot = hd.Bootstrap("hyde-boot.txt")

  **Method**: ``gamma(tripl)``

    Return the gamma values for all bootstrap replicates for a triplet of taxa. Triples
    are passed to the function as a tuple with the names of the three taxa for that
    particular hypothesis test.



**Submodule**: ``visualize``
----------------------------

The ``visualize`` submodule uses the ``matplotlib`` and ``seaborn`` packages to provide basic
plotting of results.

**File**: ``viz.py``
^^^^^^^^^^^^^^^^^^^^

The ``viz.py`` file provides simple functions for plotting different distributions
of bootstrap replicate variables. It uses the ``matplotlib`` and ``seaborn`` plotting
libraries.

.. autofunction:: phyde.visualize.viz.density

.. autofunction:: phyde.visualize.viz.dist

.. autofunction:: phyde.visualize.viz.violinplot
