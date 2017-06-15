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
