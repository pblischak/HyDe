""" Class for processing HyDe output. """

from __future__ import print_function
import numpy as np

class HydeResult:
    """
    A class for reading in and working with results from a HyDe analysis
    using the C++ interface. It is mostly a simple container for storing
    triples as tuples and using them as keys in a dictionary with values
    that are also dictionaries containing the results from the HyDe analysis.
    """
    def __init__(self, infile):
        """

        """
        self.infile = infile
        self.res = {}
        self.triples = []
        self._read_hyde_results(infile)

    def _read_hyde_results(self, file):
        """

        """
        with open(file) as f:
            results = f.read().splitlines()[1:-1]
            for r in results:
                rs    = r.split('\t')
                print(rs)
                tripl = (rs[0], rs[1], rs[2])
                if tripl not in self.triples:
                    self.triples.append(tripl)
                    self.res[tripl] = self._hyde_info(rs[3:])
                else:
                    print("\nERROR:")
                    print("  The triple ", tripl, " was tested more than once.\n")

    def _hyde_info(self, b):
        """

        """
        if len(b) != 18:
            raise ValueError("** Warning: length of hyde entry is incorrect. **")
        else:
            pass
        hyde_info = {
            "Zscore" : float(b[0]),
            "Pvalue" : float(b[1]),
            "Gamma"  : float(b[2]),
            "AAAA"   : float(b[3]),
            "AAAB"   : float(b[4]),
            "AABA"   : float(b[5]),
            "AABB"   : float(b[6]),
            "AABC"   : float(b[7]),
            "ABAA"   : float(b[8]),
            "ABAB"   : float(b[9]),
            "ABAC"   : float(b[10]),
            "ABBA"   : float(b[11]),
            "BAAA"   : float(b[12]),
            "ABBC"   : float(b[13]),
            "CABC"   : float(b[14]),
            "BACA"   : float(b[15]),
            "BCAA"   : float(b[16]),
            "ABCD"   : float(b[17])
        }
        return hyde_info
