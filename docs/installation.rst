.. include:: links.rst

.. _Installation:

Installation
============

HyDe requires Python v2.7, a C++ compiler (C++11 compliant), and several Python modules (listed below).

Miniconda
---------

We recommend using a Python distribution such as Miniconda to make it easier
to manage modules and environments.

.. code:: bash

  # Get Miniconda for your operating system (Mac or Linux)
  # Answer yes to the questions the Installer asks
  # These commands will download Python 2.7 for Mac OSX
  curl -O https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
  bash Miniconda2-latest-MacOSX-x86_64.sh

Required Python packages
------------------------

With Miniconda installed, we can use ``pip`` or ``conda`` to install all of the Python
modules that HyDe requires (Cython, Numpy, Pandas, Scipy, Matplotlib, Seaborn, and Multiprocess).

.. code:: bash

  # Install packages with conda or pip command
  conda install cython numpy pandas scipy matplotlib seaborn

  # use pip to install multiprocess
  pip install multiprocess

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

    # Compile hyde_cpp
    make

    # Builds and installs phyde module
    python setup.py install

    # Test installation
    make test

Step-by-step
^^^^^^^^^^^^

After cloning HyDe from GitHub and moving into the main HyDe directory, the three
steps that follow accomplish the following tasks:

  #. ``make``: this will compile the C++ code in the ``src/`` folder and will
     will make the ``hyde_cpp`` executable.
  #. ``python setup.py install``: this will build and install the HyDe Python
     module, including the compilation of any Cython files.
  #. ``make test``: this will test the installation by running a series of commands
     designed to check that the installation was completed successfully. The main
     tests are in the ``test.py`` script in the ``test/`` folder.
