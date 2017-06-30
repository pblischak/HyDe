.. include:: links.rst

.. _Visualization:

Visualization
=============

**This page is still in the works. Sorry!**

Basic plotting
--------------
.. code:: py

  import phyde as hd

  # Read in bootstrap replicates for plotting
  boot = hd.Bootstrap("hyde-boot.txt")

  # Make a density plot of the bootstrap reps for Gamma
  hd.viz.density(boot, 'Gamma', 'sp1', 'sp2', 'sp3', title="Bootstrap Dist. of Gamma")

Advanced plotting
-----------------

For plotting that goes beyond the basic methods implemented in ``phyde``,
direct interaction with the `Seaborn <https://seaborn.pydata.org/>`_ package is required.
