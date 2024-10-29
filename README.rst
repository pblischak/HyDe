
|Documentation|  |PyPI Badge|

HyDe: Hybridization Detection Using Phylogenetic Invariants
-----------------------------------------------------------

**Publication:**

Blischak, P. D., J. Chifman, A. D. Wolfe, and L. S. Kubatko. 2018.
HyDe: A Python Package For Genome-Scale Hybridization Detection.
*Systematic Biology*, doi: `https://doi.org/10.1093/sysbio/syy023 <https://doi.org/10.1093/sysbio/syy023>`__.

`Read the Docs <https://hybridization-detection.rtfd.io/>`__
------------------------------------------------------------

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
Make sure you have the ``multiprocess`` module installed before you use them.

Graphical User Interface (Beta)
-------------------------------

We have also developed a graphical user interface (GUI) for **single-threaded**
hybridization detection analyses with HyDe. Users can supply file names and other
options in the GUI window and hit ``Run`` or ``<Return>`` to run analyses.

.. image:: docs/img/gui.png
   :align: center

Getting Help
------------

If you have questions about running HyDe, please feel free to use the
**gitter chatroom** to get help:

|Gitter|

If you have a problem while running HyDe and you think it may be a bug,
please consider filing an issue:

|HyDe Issues|

Installation
------------

Requirements:
~~~~~~~~~~~~~

-  Python 3.9+
-  Python Modules:

   -  cython
   -  numpy
   -  multiprocess
   -  progress (optional: needed for GUI)

-  C++ compiler

.. code:: bash

    # To install dependencies -- add 'progress' if you're using the GUI
    python3 -m pip install -r requirements.txt

    # Clone HyDe repository from GitHub
    git clone https://github.com/pblischak/HyDe.git
    cd HyDe

    # Now install phyde module
    python3 -m pip install .

    # Test the installation
    make test

    # Test multithreaded scripts
    make test_threads

The ``phyde`` module is also hosted on the Python Package Index (PyPI), and can be installed directly using
``pip``.

.. code:: bash

  # Install from PyPI with pip
  pip install phyde

Documentation for analyzing data using HyDe can be found `here <http://hybridization-detection.readthedocs.io/analyze.html>`_.

.. |Documentation| image:: http://readthedocs.org/projects/hybridization-detection/badge/?version=latest
   :target: http://hybridization-detection.readthedocs.io

.. |PyPI Badge| image:: https://img.shields.io/pypi/v/phyde.svg
   :target: https://pypi.python.org/pypi/phyde

.. |Gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pblischak-HyDe/Lobby

.. |HyDe Issues| image:: https://img.shields.io/badge/HyDe-Issues-blue.svg
   :target: https://github.com/pblischak/HyDe/issues
