# -*- coding: utf-8 -*-

# hyde.py
# Written by PD Blischak

from __future__ import print_function
import numpy as np
import pandas as pd
try:
    import subprocess32 as sps
except ImportError:
    import subprocess as sps
import sys
import os
from collections import Counter
from ..core.result import HydeResult
from ..core.bootstrap import Bootstrap

def run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue=0.05, bootReps=0, prefix="hyde"):
    """
    Wrapper for running HyDe C++ executable.

    Arguments
    ---------

        - infile: data file with DNA sequence data in sequential Phylip format (no header).
        - mapfile: two-column file assigning individuals to taxa.
        - outgroup: name of the outgroup taxon.
        - nInd: the number of individuals.
        - nTaxa: the number of taxa.
        - nSites: the number of sites.
        - pValue: p-value for the hypothesis test (default=0.05).
        - bootReps: number of bootstrap replicates to perform (default=0).
        - prefix: name added to output files (default="hyde").

    Returns
    -------

        If no bootstrapping is done (i.e., bootReps=0), run_hyde() returns a HydeResult
        object that contains the results for all of the hypothesis tests conduced.

        If bootstrapping is completed, the function returns a HydeResult object for the Making
        results (all hypothesis tests) and a Bootstrap object containing all of the bootstrap reps.

    Examples
    --------

        import phyde as hd

        # no bootstrapping
        res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)

        # w/ bootstrapping
        res, boot = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=1000)
    """

    if os.path.exists(prefix+"-out.txt"):
        print("\n**  Warning: File '"+prefix+"-out.txt' already exists. **")
        print("**  Renaming to 'old-"+prefix+"-out.txt'. **\n")
        os.rename(prefix+"-out.txt", "old-"+prefix+"-out.txt")
    else:
        pass

    if os.path.exists(prefix+"-boot.txt"):
        print("\n**  Warning: File '"+prefix+"-boot.txt' already exists. **")
        print("**  Renaming to 'old-"+prefix+"-boot.txt'. **\n")
        os.rename(prefix+"-boot.txt", "old-"+prefix+"-boot.txt")
    else:
        pass

    # Check that hyde executable works
    test_hyde = sps.Popen(['hyde_cpp'], stdout=sps.PIPE, stderr=sps.PIPE, shell=True)
    (out, err) = test_hyde.communicate()
    #print(repr(str(err)))
    if not str(err).startswith('\n** ERROR'):
        sys.exit("\n\n** ERROR: could not execute hyde_cpp **\n\n")
    else:
        pass

    hyde_cmd = [
        'hyde_cpp',
        "-n", str(nInd),
        "-t", str(nTaxa),
        "-s", str(nSites),
        "-i", infile,
        "-m", mapfile,
        "-o", outgroup,
        "-p", str(pValue),
        "-b", str(bootReps),
        "--prefix", prefix
    ]

    proc = sps.call(hyde_cmd)

    if bootReps > 0:
        res  = HydeResult(prefix+"-out.txt")
        boot = Bootstrap(prefix+"-boot.txt")
        return (res, boot)
    else:
        res = HydeResult(prefix+"-out.txt")
        return res

#def read_hyde_boot():
#    """
#    Read in bootstrap replicates file from HyDe.
#    """

#def make_cf_table(hyde_out, outgroup, prefix="hyde"):
#    """
#    Process HyDe output.
#    """
#    #print("TESTING: Making concordance factor table...", end='')
#    df = pd.read_csv(prefix+"-out.txt", sep='\t')
#    cf_file = open(prefix+"-cf-table.txt", 'w')
#    print("t1", "t2", "t3", "t4", "CF12_34", "CF13_24", "CF14_23", sep=',', file=cf_file)
#    df2 = df[['P1', 'Hybrid', 'P2', 'gamma']]
#    unique_taxa = np.union1d(df2.P1.unique(), df2.Hybrid.unique())
#    prev_triplet = []
#    for r in range(df2.shape[0]):
#        if Counter([df2.iloc[r][0],df2.iloc[r][1],df2.iloc[r][2]]) == Counter(prev_triplet):
#            pass
#        else:
#            print(df2.iloc[r][0], df2.iloc[r][1], df2.iloc[r][2], \
#                  outgroup, abs(df2.iloc[r][3]), 0.0, 1.0-abs(df2.iloc[r][3]), sep=',', file=cf_file)
#            prev_triplet = [df2.iloc[r][0],df2.iloc[r][1],df2.iloc[r][2]]
#
#    #print("Done.\n")
