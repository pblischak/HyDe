#!/usr/bin/env python
# -*- coding: utf-8 -*-

# hyde.py
# Written by PD Blischak

from __future__ import print_function
import numpy as np
import pandas as pd
import subprocess as sps
import scipy as sci
import sys
import os
import argparse
from collections import Counter

def run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue=0.05, bootReps=0, threads=1, prefix="hyde"):
    """
    Wrapper for running HyDe C++ executable.
    """

    if os.path.exists(prefix+"-out.txt"):
        print("\n**  Warning: File '"+prefix+"-out.txt' already exists.")
        print("**  Renaming to 'old-"+prefix+"-out.txt'.\n")
        os.rename(prefix+"-out.txt", "old-"+prefix+"-out.txt")
    else:
        pass

    # Check that hyde executable works
    test_hyde = sps.Popen(['hyde'], stdout=sps.PIPE, stderr=sps.PIPE, shell=True)
    (out, err) = test_hyde.communicate()
    #print(repr(str(err)))
    if not str(err).startswith('\n** ERROR'):
        sys.exit("\n\n** ERROR: could not execute hyde **\n\n")
    else:
        pass

    hyde_cmd = [
        "hyde",
        "-n", str(nInd),
        "-t", str(nTaxa),
        "-s", str(nSites),
        "-i", infile,
        "-m", mapfile,
        "-o", outgroup,
        "-j", str(threads),
        "-p", str(pValue),
        "-b", str(bootReps),
        "--prefix", prefix
    ]

    proc = sps.call(hyde_cmd)

#def read_hyde_boot():
#    """
#    Read in bootstrap replicates file from HyDe.
#    """

def make_cf_table(hyde_out, prefix="hyde", outgroup):
    """
    Process HyDe output.
    """
    #print("TESTING: Making concordance factor table...", end='')
    df = pd.read_csv(prefix+"-out.txt", sep='\t')
    cf_file = open(prefix+"-cf-table.txt", 'w')
    print("t1", "t2", "t3", "t4", "CF12_34", "CF13_24", "CF14_23", sep=',', file=cf_file)
    df2 = df[['P1', 'Hybrid', 'P2', 'gamma']]
    unique_taxa = np.union1d(df2.P1.unique(), df2.Hybrid.unique())
    prev_triplet = []
    for r in range(df2.shape[0]):
        if Counter([df2.iloc[r][0],df2.iloc[r][1],df2.iloc[r][2]]) == Counter(prev_triplet):
            pass
        else:
            print(df2.iloc[r][0], df2.iloc[r][1], df2.iloc[r][2], \
                  outgroup, abs(df2.iloc[r][3]), 0.0, 1.0-abs(df2.iloc[r][3]), sep=',', file=cf_file)
            prev_triplet = [df2.iloc[r][0],df2.iloc[r][1],df2.iloc[r][2]]

    #print("Done.\n")

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
    additional.add_argument('-j', '--threads', action="store", type=int, default=1,
                            metavar='\b', help="number of threads [default=1]")
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
    threads  = args.threads
    prefix   = args.prefix

    run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue, bootReps, threads, prefix)
    make_cf_table(prefix+"-out.txt", prefix, outgroup)
