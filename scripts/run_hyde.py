#!/usr/bin/env python
# -*- coding: utf-8 -*-

# run_hyde.py
# Written by PD Blischak

"""


Arguments
---------



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
    parser = argparse.ArgumentParser(description="Options for hyde.py",
                                     add_help=True)

    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--infile', action="store", type=str, required=True,
                          metavar='\b', help="name of the data input file")
    required.add_argument('-m', '--map', action="store", type=str, required=True,
                          metavar='\b', help="map of individuals to taxa")
    required.add_argument('-o', '--outgroup', action="store", type=str, required=True,
                          metavar='\b', help="name of the outgroup (only one accepted)")
    required.add_argument('-n', '--num_ind', action="store", type=str, required=True,
                          metavar='\b', help="number of individuals in data matrix")
    required.add_argument('-t', '--num_taxa', action="store", type=str, required=True,
                          metavar='\b', help="number of taxa (species, OTUs)")
    required.add_argument('-s', '--num_sites', action="store", type=str, required=True,
                          metavar='\b', help="number of sites in the data matrix")

    additional = parser.add_argument_group("additional arguments")
    additional.add_argument('-p', '--pvalue', action="store", type=float, default=0.05,
                            metavar='\b', help="p-value cutoff for test of significance [default=0.05]")
    additional.add_argument('--prefix', action="store", type=str, default="hyde",
                            metavar='\b', help="prefix appended to output files [default=hyde]")

    args     = parser.parse_args()
    infile   = args.infile
    mapfile  = args.map
    outgroup = args.outgroup
    nInd     = args.num_ind
    nTaxa    = args.num_taxa
    nSites   = args.num_sites
    pValue   = args.pvalue
    prefix   = args.prefix

    res = hd.run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue, 0, prefix)
