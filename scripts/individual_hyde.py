#!/usr/bin/env python
# -*- coding: utf-8 -*-

# individual_hyde.py
# Written by PD Blischak

"""
<< individual_hyde.py >>

Test all individuals within a putative hybrid lineage for a specified triple.

Arguments
---------

    - infile <string>   : name of the DNA sequence data file.
    - mapfile <string>  : name of the taxon map file.
    - outgroup <string> : name of the outgroup.
    - triples <string>  : name of the file containing triples for testing.
    - nind <int>        : number of sampled individuals.
    - nsites <int>      : number of sampled sites.
    - ntaxa <int>       : number of sampled taxa (populations, OTUs, etc.).

Output
------

    Writes a file ('hyde-ind.txt') with the results of the hybridization
    detection analysis for each individual tested against the putative parents
    specified in the triples file. P1 and P2 will be the name of the parents
    but the Hybrid column will contain the names of all of the individuals that
    were tested.
"""

from __future__ import print_function
import phyde as hd
import argparse
import sys
import os

def parse_triples(triples_file):
    """
    Parse a three column table or previous results file to get the names
    of the taxa that are going to be tested for hybridization.

    Returns a list of three-tuples of the form (P1, Hybrid, P2) for all triples
    that are to be tested.
    """
    triples = []
    with open(triples_file) as f:
        lines = f.read().splitlines()
        # remove header information if reading in a previous results file
        if lines[0].split()[0] == "P1" and lines[0].split()[1] == "Hybrid" and lines[0].split()[2] == "P2":
            lines = lines[1:]
        # catch the case where the last line in the file is blank
        if len(lines[-1]) == 0:
            triples = [(l.split()[0], l.split()[1], l.split()[2]) for l in lines[:-1]]
        else:
            triples = [(l.split()[0], l.split()[1], l.split()[2]) for l in lines]
    return triples

def write_ind(out, triple, outfile):
    """
    Write the current output dictionary from test_individuals()
    to the file passed as an argument to the function.
    """
    for k,v in out.items():
        print(triple[0], "\t", k, "\t", triple[2], "\t", sep='', end='', file=outfile)
        print(v['Zscore'], "\t", sep='', end='', file=outfile)
        print(v['Pvalue'], "\t", sep='', end='', file=outfile)
        print(v['Gamma'], "\t", sep='', end='', file=outfile)
        print(v['AAAA'], "\t", sep='', end='', file=outfile)
        print(v['AAAB'], "\t", sep='', end='', file=outfile)
        print(v['AABA'], "\t", sep='', end='', file=outfile)
        print(v['AABB'], "\t", sep='', end='', file=outfile)
        print(v['AABC'], "\t", sep='', end='', file=outfile)
        print(v['ABAA'], "\t", sep='', end='', file=outfile)
        print(v['ABAB'], "\t", sep='', end='', file=outfile)
        print(v['ABAC'], "\t", sep='', end='', file=outfile)
        print(v['ABBA'], "\t", sep='', end='', file=outfile)
        print(v['BAAA'], "\t", sep='', end='', file=outfile)
        print(v['ABBC'], "\t", sep='', end='', file=outfile)
        print(v['CABC'], "\t", sep='', end='', file=outfile)
        print(v['BACA'], "\t", sep='', end='', file=outfile)
        print(v['BCAA'], "\t", sep='', end='', file=outfile)
        print(v['ABCD'], "\n", sep='', end='', file=outfile)

if __name__ == "__main__":
    """
    Runs the scripts.
    """
    # print docstring if only the name of the script is given
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    parser = argparse.ArgumentParser(description="Options for individual_hyde.py",
                                     add_help=True)

    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--infile', action="store", type=str, required=True,
                          metavar='\b', help="name of the data input file")
    required.add_argument('-m', '--map', action="store", type=str, required=True,
                          metavar='\b', help="map of individuals to taxa")
    required.add_argument('-o', '--outgroup', action="store", type=str, required=True,
                          metavar='\b', help="name of the outgroup (only one accepted)")
    required.add_argument('-tr','--triples', action="store", type=str, required=True,
                          metavar='\b', help="table of triples to be analyzed using bootstrapping.")
    required.add_argument('-n', '--num_ind', action="store", type=int, required=True,
                          metavar='\b', help="number of individuals in data matrix")
    required.add_argument('-t', '--num_taxa', action="store", type=int, required=True,
                          metavar='\b', help="number of taxa (species, OTUs)")
    required.add_argument('-s', '--num_sites', action="store", type=int, required=True,
                          metavar='\b', help="number of sites in the data matrix")

    additional = parser.add_argument_group("additional arguments")
    additional.add_argument('--prefix', action="store", type=str, default="hyde",
                            metavar='\b', help="prefix appended to output files [default=hyde]")

    args     = parser.parse_args()
    infile   = args.infile
    mapfile  = args.map
    outgroup = args.outgroup
    triples  = parse_triples(args.triples)
    nind     = args.num_ind
    ntaxa    = args.num_taxa
    nsites   = args.num_sites
    prefix   = args.prefix

    # Read data into a HydeData object
    data = hd.HydeData(infile, mapfile, outgroup, nind, ntaxa, nsites)

    if os.path.exists(prefix+"-ind.txt"):
        print("\n**  Warning: File '"+prefix+"-ind.txt' already exists. **")
        print("**  Renaming to 'old-"+prefix+"-ind.txt'. **\n")
        os.rename(prefix+"-ind.txt", "old-"+prefix+"-ind.txt")
        outfile = open(prefix+"-ind.txt", 'wa')
    else:
        outfile = open(prefix+"-ind.txt", 'wa')

    print("P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\t\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n", end='', file=outfile)

    for t in triples:
        res = data.test_individuals(t[0], t[1], t[2])
        write_ind(res, t, outfile)
