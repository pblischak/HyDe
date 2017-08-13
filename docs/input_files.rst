.. include:: links.rst

.. _Input_Files:

HyDe Data Files
===============

The main data files for running HyDe are (1) the DNA sequences and (2) the
map of sampled individuals to their respective taxa/populations. Examples of these files
can be found in the ``test/`` and ``examples/`` directories within the main HyDe
folder. Below we describe each file in more detail. If you installed HyDe from PyPI
and do not have a clone of the GitHub repository, these files can be viewed `here <https://github.com/pblischak/HyDe>`__.

Sequence Data
-------------

The DNA sequence data should be in sequential Phylip format with individual names
and sequence data all on one line and separated by a tab. The header information
that is present in a typical Phylip file (# individuals and # sites) is not needed
for HyDe. Individual names do not have length restrictions (within reason) and do
not all need to be the same number of characters.

**Example:** ``data.txt``

.. code::

  Ind1  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  Ind2  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  Ind3  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  .
  .
  .
  Ind23 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  Ind24 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  Ind25 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  .
  .
  .
  Ind(N-2)  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  Ind(N-1)  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...
  IndN      AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA...

Taxon Map
---------

The taxon map file is organized similarly to the sequence data file with one individual
per row and a tab separating the individual's name from the name of the taxon/population it belongs
to. The individuals should be in the same order as the DNA sequence data file with
all individuals in a particular taxon grouped together sequentially.

**Example:** ``map.txt``

.. code::

  Ind1  Pop1
  Ind2  Pop1
  Ind3  Pop1
  .
  .
  .
  Ind23 Pop3
  Ind24 Pop3
  Ind25 Pop3
  .
  .
  .
  Ind(N-2)  PopM
  Ind(N-1)  PopM
  IndN  PopM

Triples file
------------

Each of the Python scripts that are part of HyDe take a file
of triples that are to be tested for hybridization. The format of this triples
file is a three column table where each row is a triple to be tested. Column one
of the table is parent one, column two is the putative hybrid,
and column three is parent two.
You can name the populations anything you like as long as the name doesn't have spaces.

.. code::

  Pop1  Pop2  Pop3
  Pop1  Pop3  Pop4
  Pop2  Pop4  Pop5
  .
  .
  .
  Pop4  Pop5  PopM

We have also written the scripts to take previous output files from HyDe as input.
It does this by ignoring the first-row header and using the triples specified in
the rows that follow.
