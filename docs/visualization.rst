.. include:: links.rst

.. _Visualization:

Visualization
=============

Basic plotting
--------------
.. code:: py

  import phyde as hd

  # Read in bootstrap replicates for plotting
  boot = hd.Bootstrap("hyde-boot.txt")


  vz = hd.HydeViz(boot)
  hd.viz.density(boot('Gamma', 'sp1', 'sp2', 'sp3'))

Advanced plotting
-----------------

For plotting that goes beyond the basic methods implemented in the ``HydeViz`` class,
direct interaction with the `Seaborn <https://seaborn.pydata.org/>`_ package is required.
