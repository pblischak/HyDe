.. include:: links.rst

.. _API:

API Reference
=============

``core`` module
---------------

The ``core`` module provides the Python interface for hybridization detection analyses
and data exploration through the use implemented classes.

``data.pyx``
^^^^^^^^^^^^

**Class**: ``HydeData``

  The ``HydeData`` class implements methods for reading in DNA sequence data and
  taxon maps for hybridization detection analyses.

  **Constructor**: ``HydeData(data, map, outgroup, num_ind, num_taxa, num_sites)``

    The constructor for the ``HydeData`` class

      **Example**

      .. code:: py

        import hyde as hd
        data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)

  **Method**: ``test_triple(p1, hyb, p2)``

    This is the main method for testing for hybridization at the population level.

    **Example**

    .. code:: py

      import hyde as hd
      data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
      res = data.test_triple("sp1", "sp2", "sp3")

  **Method**: ``test_individual(p1, hyb, p2)``

    This is the main method for testing each individual in the specified hybrid
    population.

    **Example**

    .. code:: py

      import hyde as hd
      data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
      res = data.test_individual("sp1", "sp2", "sp3")

``result.py``
^^^^^^^^^^^^^

**Class**: ``HydeResult``

  **Constructor**: ``HydeResult(infile)``

    **Example**:

    .. code:: py

      import hyde as hd
      res = hd.HydeResult("hyde-out.txt")

``bootstrap.py``
^^^^^^^^^^^^^^^^

  ``Bootstrap`` class

``analyze`` module
------------------

``main.py``
^^^^^^^^^^^

  **Function**: ``run_hyde(data, map, outgroup, num_ind, num_taxa, num_sites, boot_reps=0)``

``visualize`` module
--------------------

``viz.py``
^^^^^^^^^^

  **Function**: ``plt()``
