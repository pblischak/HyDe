.. include:: links.rst

.. _Analyze:

Analyzing Data with HyDe
========================

Below we provide details on the three main scripts that are used to analyze data with HyDe.
The multithreaded versions of these scripts behave in the exact same way, but have an added
``--threads`` (``-j``) option to specify how many threads to use.

We will be using the ``data.txt`` and ``map.txt`` files in the ``test/`` folder from
the GitHub repo for HyDe. If you don't have a clone of the repo, you can download
the files using the following commands:

.. code:: bash

  curl -O https://raw.githubusercontent.com/pblischak/HyDe/master/test/data.txt
  curl -O https://raw.githubusercontent.com/pblischak/HyDe/master/test/map.txt


.. note:: **Recommended workflow**:

  When analyzing data with HyDe, we recomend the following workflow.

  #. Use the ``run_hyde.py`` script to analyze all possible triples. This will produce
     two output files, one with all results and one with only significant results
     (see **Note** below on filtered results).
  #. Next, if you want to see if certain individuals are hybrids, run the ``individual_hyde.py``
     script and use the filtered results file output from the previous step as the triples file.
  #. For the ``bootstrap_hyde.py`` script, we recommend using it when you don't have enough
     data (and therefore not enough power) to detect hybridization within a single individual.
     This will depend on how much hybridization has occurred, but it is typically difficult to
     detect hybridization with HyDe using less than 10,000 sites. Otherwise, we always recommend
     using the ``individual_hyde.py`` script.

Command Line Scripts
--------------------

``run_hyde.py``
^^^^^^^^^^^^^^^

To run HyDe from the command line, we have provided a Python script (``run_hyde.py``) that will
test all triples in all directions ("full" analysis).
It can also be used to test a predefined set of hypothesis tests using a triples file.
The arguments for the script are passed using command line flags, all of which can be viewed
by typing ``run_hyde.py -h``. Typing the name of the script with no arguments will print
out a docstring with additional details.

.. code:: bash

  # Run a full hybridization detection analysis
  run_hyde.py -i data.txt -m map.txt -o out -n 16 -t 4 -s 50000

The results will be written to file with a prefix that can be supplied
using the ``--prefix`` flag (``<prefix>-out.txt``; the default is 'hyde').

.. note:: **Filtered results**:

  We also write a file (``<prefix>-out-filtered.txt``) that filters
  the results from the hybridization detection
  analysis to only include significant results with sensible values of :math:`\gamma`
  (:math:`0 < \gamma < 1`). Some values of :math:`\gamma` in the original results
  file may be nonsensical because they will be either negative or greater than 1. However,
  :math:`\gamma` does not have any theoretical limits with regard to the hypothesis test in HyDe:
  it can range from :math:`-\infty` to :math:`\infty`. The reason that these values occur and may
  give a significant p-value is because they are testing a hypothesis that involves a hybrid but
  in the wrong direction (i.e., a hybrid is tested as one of the parental species). Testing all possible
  directions that hybridization can occur for three taxa is typically what causes these types of results to happen.


``individual_hyde.py``
^^^^^^^^^^^^^^^^^^^^^^

**Triples file required.**

If you want to test for hybridization at the individual level within the populations
that have significant levels of hybridization, you can use the filtered results file from
an analysis with ``run_hyde.py`` as the input triples file for the ``individual_hyde.py``
script. The options for ``individual_hyde.py`` are the same as for ``run_hyde.py``,
and typing the name of the script without arguments will again print out a
docstring with more details.

.. code:: bash

  # Test all individuals for the hybrid populations specified
  # in the file hyde-filtered-out.txt
  individual_hyde.py -i data.txt -m map.txt -tr hyde-out-filtered.txt -o out
                     -n 16 -t 4 -s 50000

The results of the individual level tests will be written to a file called ``<prefix>-ind.txt``,
where ``<prefix>`` can be set using the ``--prefix`` flag (default='hyde').

``bootstrap_hyde.py``
^^^^^^^^^^^^^^^^^^^^^

**Triples file required.**

The ``bootstrap_hyde.py`` conducts bootstrap resampling of individuals within
hybrid populations. The arguments are the same as the ``individual_hyde.py`` script
with the addition of specifying the number of bootstrap replicates (``--reps=<#reps>``; default=100).

.. code:: bash

  # Bootstrap resample individuals within hybrid populations
  # specified in hyde-out-filtered.txt
  bootstrap_hyde.py -i data.txt -m map.txt -tr hyde-out-filtered.txt -o out
                    -n 16 -t 4 -s 50000 --reps 200

The output file is named ``hyde-boot.txt``, but again this can be changed using the
``--prefix`` argument. Bootstrap replicates for each triple are separated by a line
with four pound symbols and a newline ("####\\n"; match this pattern to split results).

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

The final method of note implemented by the HydeData class is the ``bootstrap_triple()`` function,
which will randomly resample individuals in the hybrid population and run hypothesis tests for each
replicate. The number of bootstrap replicates to run is specified using the ``reps=#`` argument.

.. code:: py

  res_boot = dat.bootstrap_triple("sp1", "sp2", "sp3", reps=100)

``res_boot`` stores the results of this analysis as set of nested dictionaries.
The easiest way to see its structure is to print it using ``print(res_boot)``.

More detailed documentation on the HydeData class and other classes implemented in
the ``phyde`` module can be found in the :ref:`API`.
