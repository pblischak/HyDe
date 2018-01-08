.. include:: links.rst

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

.. autoclass:: phyde.HydeData
  :members:

**File**: ``result.py``
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: phyde.HydeResult
  :members:
  :special-members: __call__

**File**: ``bootstrap.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: phyde.Bootstrap
  :members:
  :special-members: __call__

**Submodule**: ``visualize``
----------------------------

The ``visualize`` submodule uses the ``matplotlib`` and ``seaborn`` packages to provide basic
plotting of results stored by the ``phyde.Bootstrap`` class.

**File**: ``viz.py``
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: phyde.visualize.viz.density

.. autofunction:: phyde.visualize.viz.dist

.. autofunction:: phyde.visualize.viz.violinplot
