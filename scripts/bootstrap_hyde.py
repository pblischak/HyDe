#!/usr/bin/env python
# -*- coding: utf-8 -*-

# bootstrap_hyde.py
# Written by PD Blischak

"""
Bootstrap resample individuals within a putative hybrid population for a
specified triple.

Specifying triples
------------------

The triples to test can be specified using a previous results
file (filtered for significant results):

```
P1  Hybrid  P2  Zscore  ...
sp1 sp2 sp3 2.321934    ...
sp1 sp2 sp4 14.38293    ...
```

or they can be tested using a user-made table with three
columns (P1, Hybrid, P2) and no headers:

```
sp1 sp2 sp3
sp1 sp2 sp4
```

Arguments
---------

    - infile <string>   : name of the DNA sequence data file.
    - mapfile <string>  : name of the taxon map file.
    - outgroup <string> : name of the outgroup.
    - triples <string>  : name of the file containing triples for testing.
    - reps <int>        : number of bootstrap replicates (default=100).
    - nind <int>        : number of sampled individuals.
    - nsites <int>      : number of sampled sites.
    - ntaxa <int>       : number of sampled taxa (populations, OTUs, etc.).

Output
------


"""

from __future__ import print_function
import phyde as hd
import argparse

def parse_triples(triples_file):
    """
    Parse a three column table or previous results file to get the names
    of the taxa that are going to be tested for hybridization.

    Returns a list of three-tuples of the form (P1, Hybrid, P2) for all triples
    that are to be tested.
    """
    triples = []
    with open(triples_file) as f:
        lines = f.splitlines()
        # remove header information if reading in a previous results file
        if lines[0].split()[0] == "P1" and lines[0].split()[1] == "Hybrid" and lines[0].split()[2] == "P2":
            lines = lines[1:]
        # catch the case where the last line in the file is blank
        if len(lines[-1]) == 0:
            triples = [(l.split()[0], l.split()[1], l.split()[2]) for l in lines[:-1]]
        else:
            triples = [(l.split()[0], l.split()[1], l.split()[2]) for l in lines]
    return triples

if __name__ == "__main__":
    """
    Runs the script.
    """
    print("Not yet implemented. Sorry!")
