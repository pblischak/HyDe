.. include:: links.rst

.. _Installation:

Installation
============

Miniconda
---------

.. code:: bash

  # Get Miniconda for your operating system (Mac or Linux)
  # Answer yes to the questions the Installer asks
  # These commands will download Python 2.7 for Mac OSX
  curl -O https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
  bash Miniconda2-latest-MacOSX-x86_64.sh

Required Python packages
------------------------

.. code:: bash

  # Install packages with conda or pip command
  conda install cython numpy pandas scipy matplotlib seaborn
  # pip install cython numpy pandas scipy matplotlib seaborn

Installing HyDe
---------------

We'll install the HyDe module by cloning it from GitHub using the commands below:

.. code:: bash

    # Clone the HyDe repo from GitHub and cd into the repo
    git clone https://github.com/pblischak/HyDe.git
    cd HyDe/

    # Compiles hyde_cpp
    make

    # Builds and installs phyde module
    python setup.py install

    # Tests installation
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
     designed to check that the installation was completed succesfully. The main
     tests are in the ``test.py`` script in the ``test/`` folder.
