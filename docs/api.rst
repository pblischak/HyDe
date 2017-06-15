import phyde.. include:: links.rst

.. _API:

API Reference
=============

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

**Class**: ``HydeData``

  The ``HydeData`` class is an extension type that implements methods for reading
  in DNA sequence data and taxon maps for hybridization detection analyses.

  **Constructor**: ``HydeData(data, map, outgroup, num_ind, num_taxa, num_sites)``

    The constructor for the ``HydeData`` class.

    **Example**

    .. code:: py

      import phyde as hd
      data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)

  **Method**: ``test_triple(p1, hyb, p2)``

    This is the main method for testing for hybridization at the population level.

    **Example**

    .. code:: py

      import phyde as hd
      data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
      res = data.test_triple("sp1", "sp2", "sp3")

  **Method**: ``test_individual(p1, hyb, p2)``

    This is the main method for testing each individual in the specified hybrid
    population.

    **Example**

    .. code:: py

      import phyde as hd
      data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
      res = data.test_individual("sp1", "sp2", "sp3")

**File**: ``result.py``
^^^^^^^^^^^^^^^^^^^^^^^

**Class**: ``HydeResult``

  **Constructor**: ``HydeResult(infile)``

    **Example**

    .. code:: py

      import phyde as hd
      res = hd.HydeResult("hyde-out.txt")

**File**: ``bootstrap.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Class**: ``Bootstrap``

  **Constructor**: ``Bootstrap(bootfile)``

    Read in bootstrap replicates from a HyDe analysis.

    **Example**

    .. code:: py

      import phyde as hd
      boot = hd.Bootstrap("hyde-boot.txt")

  **Method**: ``gamma(tripl)``

    Return the gamma values for all bootstrap replicates for a triplet of taxa. Triples
    are passed to the function as a tuple with the names of the three taxa for that
    particular hypothesis test.

    **Example**

    .. code:: py

      import phyde as hd
      boot = hd.Bootstrap("hyde-boot.txt")
      boot.gamma(("sp1", "sp2", "sp3"))

**Submodule**: ``analyze``
--------------------------

The ``analyze`` submodule provides wrapper functions for running a full HyDe
analysis on all possible triples using the ``hyde_cpp`` executable. It also reads
the results back into Python using the classes from the ``core`` submodule so that
they are available from within Python.

**File**: ``main.py``
^^^^^^^^^^^^^^^^^^^^^

  **Function**: ``run_hyde(data, map, outgroup, num_ind, num_taxa, num_sites, boot_reps=0)``

**Submodule**: ``visualize``
----------------------------

The ``visualize`` submodule uses the ``matplotlib`` package to provide basic
plotting of results.

**File**: ``viz.py``
^^^^^^^^^^^^^^^^^^^^

**Class**: ``HydeViz``

The HydeViz class provides a simplified interface for plotting different distributions
of bootstrap replicate variables. It uses the ``matplotlib`` and ``seaborn`` plotting
libraries.

  **Constructor**: ``HydeViz(boot_obj)``

    **Example**

    .. code:: py

      import phyde as hd
      boot = hd.Bootstrap("hyde-boot.txt")
      vz   = hd.HydeViz(boot)



  **Function**: ``density(arr, **kwargs)``
