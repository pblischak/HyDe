#!/usr/bin/env python
# -*- coding: utf-8 -*-

# hyde.py
# Written by PD Blischak

import phyde
import argparse

if __name__ == "__main__":
    """
    Run the main script.
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
    additional.add_argument('-b', '--bootstrap', action="store", type=int, default=0,
                            metavar='\b', help="number of bootstrap replicates [default=0]")
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
    bootReps = args.bootstrap
    prefix   = args.prefix

    if bootReps > 0:
        res, boot = phyde.run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue, bootReps, prefix)
    else:
        res = phyde.run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue, bootReps, prefix)
#    hyde.make_cf_table(prefix+"-out.txt", outgroup, prefix)
