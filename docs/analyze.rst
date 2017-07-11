.. include:: links.rst

.. _Analyze:

Analyzing Data with HyDe
========================

HyDe can be used to analyze data in a couple different ways. Analyses can be run from the command line
using either the ``run_hyde.py`` script or the ``hyde_cpp`` executable itself. Analyses can also be run
from within Python by importing the ``phyde`` module, reading in your sequence data and map file, and
running hypothesis tests at the population or individual level
(``test_triple()`` and ``test_individual()``, respectively). We will start by going through
how to analyze data from within Python, then we will cover running analyses directly from the command line.

We will be using the ``data.txt`` and ``map.txt`` files in the ``test/`` folder from
the GitHub repo for HyDe. If you don't have a clone of the repo, you can download
the files using the following commands:

.. code:: bash

  curl -O https://raw.githubusercontent.com/pblischak/HyDe/master/test/data.txt
  curl -O https://raw.githubusercontent.com/pblischak/HyDe/master/test/map.txt

Python Interface
----------------

Reading in Data
^^^^^^^^^^^^^^^

Reading data files (DNA sequences and taxon maps) into Python is done using the
``HydeData`` class. Making a new variable using the class requires passing six arguments
to the constructor (in order): (1) the name of the data file, (2) the name of the map file,
(3) the name of the outgroup taxon, (4) the number of individuals, (5) the number of taxa, and
(6) the number of sites. Names should be provided in quotes. The code below will read in the
``data.txt`` and ``map.txt`` files for us to analyze.

.. code:: py

  # Import the phyde module. For simplicity, we always import it as `hd`
  import phyde as hd

  # Read in the data using the HydeData class
  dat = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)

Conducting Individual Hypothesis Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With our data read in and stored using the variable ``dat``, we can begin to run hypothesis tests using
the methods provided in the ``HydeData`` class. The first of these is the ``test_triple()`` method, which will
conduct a hypothesis test at the population level for a specified triple of taxa (remember, the outgroup has
already been specified when we read in the data).

.. code:: py

  # Using the `dat` variable from the previous code section,
  # we'll run a hypothesis test on the triple (sp1, sp2, sp3).

  res = dat.test_triple("sp1", "sp2", "sp3")

``res`` is a variable that stores the results of our hypothesis test. More specifically, it is
a Python dictionary. To see what it contains we can type ``print res`` (or ``print(res)`` if you have
imported the ``print_function`` from the ``__future__`` module).

Using the same ``HydeData`` variable, we can test all of the individuals in the taxon "sp2" to see if
they are all hybrids. We do this using the ``test_individuals()`` method.

.. code:: py

  # The code here should look very similar to the previous code block
  # The only difference is that we are calling a different method

  res_ind = dat.test_individuals("sp1", "sp2", "sp3")

``res_ind`` stores the results of the hypothesis tests for each individual in population "sp2"
(individuals "i5" through "i10"). Since the results for each test is a dictionary, the ``res_ind``
variable is a dictionary of dictionaries with the individual names as the keys and the results of the
hypothesis test as the associated value. The code below shows a few examples of how to we can
work with these dictionary results in Python.

.. code:: py

  # Look at results for individual i5 in res_ind
  res_ind["i5"]

  # To look at a specific value, say "Gamma", we need two keys: one for each nested level of the dictionaries
  res_ind["i5"]["Gamma"]

  # To get all of the values of the test statistics ("Zscore"), we can use a dictionary comprehension
  zscores = {k:v["Zscore"] for k,v in res_ind.items()}
  print zscores

Testing All Triples
^^^^^^^^^^^^^^^^^^^

If you want to test all possible triples in the data set from within Python, we provide a wrapper function
that will call the ``hyde_cpp`` executable and will return the result back to Python for you to work with.
The ``run_hyde()`` function will return the main results as a HydeResult object. If you conduct bootstrapping,
it will also return a Bootstrap object. Using the ``run_hyde()`` function is similar to making a HydeData
object and the arguments are in the same order (data file, map file, outgroup,
number of individuals, number of taxa, and number of sites). To run bootstrapping,
add the ``bootReps=<#BOOTREPS>`` argument.

.. code:: py

  # Import the phyde module if you need to do so
  # import phyde as hd

  # If you aren't doing any bootstrapping, run_hyde() will only return a HydeResult object
  res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)

  # If you do run bootstrapping, run_hyde() will return a HydeResult and a Bootstrap object
  res, boot = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=100)

The results from using the ``run_hyde()`` function are also written to file
(``hyde-result.txt`` for the main results and ``hyde-boot.txt`` for the bootstrap replicates).
You can change the "hyde" prefix by supplying an extra argument, ``prefix="<PREFIX>"``, to ``run_hyde()``.

More detailed documentation on the HydeData, HydeResult, and Bootstrap classes can be found in the :ref:`API`.

Command Line Interface
----------------------

``run_hyde.py``
^^^^^^^^^^^^^^^

To run HyDe from the command line, we have provided a Python script (``run_hyde.py``) that will
test all triples in the same way that the ``run_hyde()`` function does from within Python.
The arguments for the script are passed using command line flags, all of which can be viewed
by typing ``run_hyde.py -h``.

.. code:: bash

  # Run HyDe with 100 bootstrap reps
  run_hyde.py -i data.txt -m map.txt -o out -n 16 -t 4 -s 50000 -b 100

The results will again be written to files with a prefix that can be supplied
using the ``--prefix`` flag (the default is 'hyde').

``hyde_cpp``
^^^^^^^^^^^^

If the Python module is not working for some reason (e.g., installation problems),
the ``hyde_cpp`` executable can actually function as a standalone program for
conducting hybridization detection analyses. It takes the same command line options
as the ``run_hyde.py`` script and outputs the same files. A makefile to compile the
executable is provided in the main HyDe folder.

.. code:: bash

  # View options using -h
  hyde_cpp -h

  # Run HyDe with 100 bootstrap reps with hyde_cpp
  hyde_cpp -i data.txt -m map.txt -o out -n 16 -t 4 -s 50000 -b 100
