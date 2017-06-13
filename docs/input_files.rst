.. include:: links.rst

.. _Input_Files:

HyDe Data Files
===============

The main data files for running HyDe are (1) the DNA sequences and (2) the
map of sampled individuals to their respective taxa/OTUs. Examples of these files
can be found in the ``test/`` and ``examples/`` directories within the main HyDe
folder. Below we describe each file in more detail.

Sequence Data
-------------

The DNA sequence data should be in sequential Phylip format with individual names
and sequence data all on one line and separated by a tab. The header information
that is present in a typical Phylip file (# individuals and # sites) is not needed
for HyDe. Individual names do not have length restrictions (within reason) and do
not all need to be the same number of characters.

**Example:** ``data.txt``

.. code::

  Ind1  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  Ind2  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  Ind3  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  .
  .
  .
  Ind23 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  Ind24 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  Ind25 AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  .
  .
  .
  Ind(N-2)  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  Ind(N-1)  AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA
  IndN      AGATTGAGCTAGCAGACGTGACAGACAGATGACAGTGACGA

Taxon Map
---------

The taxon map file is organized similarly to the sequence data file with one individual
per row and a tab separating the individual's name from the name of the taxon it belongs
to. The individuals should be in the same order as the DNA sequence data file with
all individuals in a particular taxon grouped together sequentially.

**Example:** ``map.txt``

.. code::

  Ind1  Taxon1
  Ind2  Taxon1
  Ind3  Taxon1
  .
  .
  .
  Ind23 Taxon3
  Ind24 Taxon3
  Ind25 Taxon3
  .
  .
  .
  Ind(N-2)  TaxonM
  Ind(N-1)  TaxonM
  IndN  TaxonM
