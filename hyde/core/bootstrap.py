""" Class for processing HyDe bootsrapping output. """

from __future__ import print_function
import numpy as np
import pandas as pd

class Bootstrap:
    """
    A class representing a dictionary of
    species triplets with info for each
    bootstrap replicate stored as a dictionary
    containing parameter values.
    """
    def __init__(self, bootfile):
        self.bootfile = bootfile
        self.brs = {}
        self.tps = []
        self._read_bootstraps(bootfile)

    def _read_bootstraps(self, bootfile):
        with open(bootfile) as b:
            boot_reps = b.read().split("####\n")
            split_bootreps = [r.splitlines() for r in boot_reps]
            for s in split_bootreps:
                lines = [line.split('\t') for line in s[1:]]
                for l in lines:
                    tripl = (l[0],l[1],l[2])
                    if tripl not in self.tps:
                        self.tps.append(tripl)
                        self.brs[tripl] = []
                    else:
                        pass
                    self.brs[tripl].append(self._boot_info(l[3:]))

    def _boot_info(self, b):
        if len(b) != 18:
            print(len(b), "** Warning: length of bootrep entry is incorrect. **")
        else:
            pass
        boot_info = {
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
        return boot_info

    def summarize(self, tripl):
        """

        """
        summary = {
            "Zscore" : [],
            "Pvalue" : [],
            "Gamma"  : []
        }
        for t in self.brs[tripl]:
            summary['Zscore'].append(t['Zscore'])
            summary['Pvalue'].append(t['Pvalue'])
            summary['Gamma'].append(t['Gamma'])
        return summary
