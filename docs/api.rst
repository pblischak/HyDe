.. include:: links.rst

.. _API:

API Reference
=============

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
