.. HyDe documentation master file, created by
   sphinx-quickstart on Sun Apr 23 10:50:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HyDe: hybridization detection using phylogenetic invariants
===========================================================

|Build Status| |Documentation|

HyDe is a Python package that wraps a primary C++ program for running a
hybridization detection analysis. We have provided a ``Makefile`` that
will compile the C++ executable (``hyde_cpp``) and will then install the
HyDe Python package using the ``setup.py`` file. To ensure that the necessary
dependencies are available, we suggest using a Python distribution such
as `Miniconda <https://conda.io/miniconda.html>`__.

Documentation
=============

.. toctree::
   :maxdepth: 1

   installation.rst
   input_files.rst
   workflow.rst
   explore.rst
   analyze.rst
   visualization.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest
