.. HyDe documentation master file, created by
   sphinx-quickstart on Sun Apr 23 10:50:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HyDe: hybridization detection using phylogenetic invariants
===========================================================

|Build Status| |Documentation| |PyPI Badge| |Gitter|

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

A tutorial on using HyDe is available on GitHub: `HyDe Tutorial <http://nbviewer.jupyter.org/github/pblischak/evol2017/blob/master/HyDe.ipynb>`__.

Features
--------

- Conduct hypothesis tests using multiple individuals per taxon.
- Test each individual within a putative hybrid lineage to assess if
  hybridization is uniform.
- Test all possible triples of taxa and process results from within Python.
- Bootstrap individuals within taxa to assess patterns of hybrid speciation vs.
  introgression.
- Visualize the distributions of various quantities (Test Statistic, Hybridization Parameter, D-Statistic)
  using bootstrap replicates.
- Calculate the D-Statistic (ABBA-BABA) using site pattern counts returned during
  a hypothesis test.

Documentation
=============

.. toctree::
   :maxdepth: 1

   installation.rst
   input_files.rst
   analyze.rst
   visualization.rst
   api.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest

.. |PyPI Badge| image:: https://badge.fury.io/py/phyde.svg
   :target: https://pypi.python.org/pypi/phyde

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-HyDe/Lobby
