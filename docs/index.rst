.. HyDe documentation master file, created by
   sphinx-quickstart on Sun Apr 23 10:50:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HyDe: Hybridization Detection Using Phylogenetic Invariants
===========================================================

|Build Status| |Documentation| |PyPI Badge|

HyDe is a software package that detects hybridization in phylogenomic
data sets using phylogenetic invariants. The primary interface for HyDe is a Python
module called ``phyde`` (**P**\ ythonic **Hy**\ bridization **De**\ tection).
``phyde`` provides a suite of tools for performing hypothesis tests on triples of taxa
to detect hybridization. To ensure that the necessary
dependencies for ``phyde`` are available, we suggest using a Python distribution such
as `Miniconda <https://conda.io/miniconda.html>`__.

To facilitate analyses using the Python module, three scripts are provided to
conduct hybridization detection analyses directly from the command line:

- ``run_hyde.py``: runs a standard hybridization detection analysis on all triples
  in all directions. Results are also filtered based on whether there is significant
  evidence for hybridization.
- ``individual_hyde.py``: tests each individual within a putative hybrid population
  using a list of specified triples specified.
- ``bootstrap_hyde.py``: conducts bootstrap resampling of the individuals within
  the putative hybrid lineages for each specified triple.

These last two scripts need to be given a three column table of triples
(P1, Hybrid, P2) that you wish to test:

.. code::

  sp1 sp2 sp3
  sp1 sp3 sp4
  sp3 sp4 sp5
  .
  .
  .

You can also use a results file from a previous analysis as a triples file.
For example, you can use the filtered results from the ``run_hyde.py`` script so that
you only run analyses on triples that have significant levels of hybridization.
If you only have a few hypotheses that you want to test, then you can also pass
a triples file to ``run_hyde.py`` and it will only test those hypotheses rather than
testing everything.

Multithreaded versions of these scripts are also available (``run_hyde_mp.py``,
``individual_hyde_mp.py``, and ``bootstrap_hyde_mp.py``).
Make sure you have the ``multiprocess`` module installed before you use them:
``pip install multiprocess``.

Getting Help
------------

If you have questions about running HyDe, please feel free to use the
**gitter chatroom** to get help:

|Gitter|

If you have a problem while running HyDe and you think it may be a bug,
please consider filing an issue on GitHub:

|HyDe Issues|

Features
--------

- Conduct hypothesis tests using multiple individuals per population.
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
   gui.rst
   api.rst

.. |Build Status| image:: https://travis-ci.org/pblischak/HyDe.svg?branch=master
   :target: https://travis-ci.org/pblischak/HyDe

.. |Documentation| image:: https://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io/en/latest/?badge=latest

.. |PyPI Badge| image:: https://img.shields.io/pypi/v/phyde.svg
   :target: https://pypi.python.org/pypi/phyde

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-HyDe/Lobby

.. |HyDe Issues| image:: https://img.shields.io/badge/HyDe-Issues-blue.svg
   :target: https://github.com/pblischak/HyDe/issues
