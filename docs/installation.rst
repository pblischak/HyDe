.. include:: links.rst

.. _Installation:

Installation
============

HyDe requires Python v3.6+, a C++ compiler, and several Python modules (listed below).

Miniconda
---------

We recommend using a Python distribution such as Miniconda to make it easier
to manage modules and environments.

.. code:: bash

  # Get Miniconda for your operating system (Mac or Linux)
  # Answer yes to the questions the Installer asks
  # These commands will download Python 3.6 for Mac OSX
  curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh

Required Python Packages
------------------------

With Miniconda installed, we can use ``pip`` (or ``conda``) to install all of the Python
modules that HyDe requires (Cython, Numpy, and Multiprocess).

.. code:: bash

  # Install packages with pip
  conda install cython numpy multiprocess

Installing HyDe
---------------

There are two ways that HyDe can be installed once you have all of the required Python modules:
(1) install from PyPI using pip or (2) clone from GitHub and install manually.

To install from PyPI, all we need to do is type the following command:

.. code:: bash

  #installing from PyPI
  pip install phyde

Next, we'll take a look at how to install HyDe by cloning it from GitHub.
The commands below take you through every step to accomplish this:

.. code:: bash

    # Clone the HyDe repo from GitHub and cd into the repo
    git clone https://github.com/pblischak/HyDe.git
    cd HyDe/

    # Builds and installs phyde module
    python3 -m pip install .

    # Test installation
    make test

Step-By-Step
^^^^^^^^^^^^

After cloning HyDe from GitHub and moving into the main HyDe directory, the next two
steps accomplish the following tasks:

  #. ``python3 -m pip install .``: this will build and install the HyDe Python
     module, including the compilation of any Cython files (C++ compiler required).
  #. ``make test``: this will test the installation by running a series of commands
     designed to check that the installation was completed successfully. The main
     tests are in the ``test.py`` script in the ``test/`` folder.
