from __future__ import print_function
import numpy as np
import pandas as pd
try:
    import subprocess32 as sps
except ImportError:
    import subprocess as sps
import sys
import os
#from collections import Counter
from ..core.result import HydeResult
#from ..core.bootstrap import Bootstrap

def run_hyde(infile, mapfile, outgroup, nInd, nTaxa, nSites, pValue=0.05, prefix="hyde"):
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
        - prefix: name added to output files (default="hyde").

    Returns
    -------

        run_hyde() returns a HydeResult object that contains the results for
        all of the hypothesis tests conduced.

    Examples
    --------

        import phyde as hd
        res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)
    """

    if os.path.exists(prefix+"-out.txt"):
        print("\n**  Warning: File '"+prefix+"-out.txt' already exists. **")
        print("**  Renaming to 'old-"+prefix+"-out.txt'. **\n")
        os.rename(prefix+"-out.txt", "old-"+prefix+"-out.txt")
    else:
        pass

    #if os.path.exists(prefix+"-boot.txt"):
    #    print("\n**  Warning: File '"+prefix+"-boot.txt' already exists. **")
    #    print("**  Renaming to 'old-"+prefix+"-boot.txt'. **\n")
    #    os.rename(prefix+"-boot.txt", "old-"+prefix+"-boot.txt")
    #else:
    #    pass

    # Check that hyde executable works
    test_hyde = sps.Popen(['hyde_cpp'], stdout=sps.PIPE, stderr=sps.PIPE, shell=True)
    (out, err) = test_hyde.communicate()
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
        "--prefix", prefix
    ]

    proc = sps.call(hyde_cmd)

    #if bootReps > 0:
    #    res  = HydeResult(prefix+"-out.txt")
    #    boot = Bootstrap(prefix+"-boot.txt")
    #    return (res, boot)
    #else:
    res = HydeResult(prefix+"-out.txt")
    return res
