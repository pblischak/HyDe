""" Class for processing HyDe bootsrapping output. """

from __future__ import print_function
import numpy as np

class Bootstrap:
    """
    A class representing a dictionary of
    species triplets with info for each
    bootstrap replicate stored as a dictionary
    containing parameter values.
    """
    def __init__(self, bootfile):
        self.bootfile = bootfile
        self.breps = {}
        self.triples    = []
        self._read_bootstraps(bootfile)

    def _read_bootstraps(self, bootfile):
        with open(bootfile) as b:
            boot_reps = b.read()[:-1].split("####\n")
            print("Number of boot reps:", len(boot_reps))
            split_bootreps = [r.splitlines() for r in boot_reps]
            for s in split_bootreps:
                lines = [line.split('\t') for line in s[1:]]
                for l in lines:
                    tripl = (l[0],l[1],l[2])
                    if tripl not in self.triples:
                        self.triples.append(tripl)
                        self.breps[tripl] = []
                    else:
                        pass
                    self.breps[tripl].append(self._boot_info(l[3:]))

    def _boot_info(self, b):
        if len(b) != 18:
            raise ValueError("** Warning: length of bootrep entry is incorrect. **")
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

    def summarize(self):
        """

        """
        summaries = {}
        for t in self.triples:
            summaries[t] = {}
            summaries[t]["Zscore"] = [np.mean(self.zscore(t)), np.std(self.zscore(t))]
            summaries[t]["Pvalue"] = [np.mean(self.pvalue(t)), np.std(self.pvalue(t))]
            summaries[t]["Gamma"]  = [np.mean(self.gamma(t)), np.std(self.gamma(t))]
        return summaries

    def gamma(self, tripl):
        return [t['Gamma'] for t in self.breps[tripl]]

    def zscore(self, tripl):
        return [t['Zscore'] for t in self.breps[tripl]]

    def pvalue(self, tripl):
        return [t['Pvalue'] for t in self.breps[tripl]]

    def site_patterns(self, tripl):
        return {
            "AAAA" : [t['AAAA'] for t in self.breps[tripl]],
            "AAAB" : [t['AAAB'] for t in self.breps[tripl]],
            "AABA" : [t['AABA'] for t in self.breps[tripl]],
            "AABB" : [t['AABB'] for t in self.breps[tripl]],
            "AABC" : [t['AABC'] for t in self.breps[tripl]],
            "ABAA" : [t['ABAA'] for t in self.breps[tripl]],
            "ABAB" : [t['ABAB'] for t in self.breps[tripl]],
            "ABAC" : [t['ABAC'] for t in self.breps[tripl]],
            "ABBA" : [t['ABBA'] for t in self.breps[tripl]],
            "BAAA" : [t['BAAA'] for t in self.breps[tripl]],
            "ABBC" : [t['ABBC'] for t in self.breps[tripl]],
            "CABC" : [t['CABC'] for t in self.breps[tripl]],
            "BACA" : [t['BACA'] for t in self.breps[tripl]],
            "BCAA" : [t['BCAA'] for t in self.breps[tripl]],
            "ABCD" : [t['ABCD'] for t in self.breps[tripl]],
        }

    def write_summary(self, summary_file):
        print("Writing summary to file:", summary_file)
        summ = self.summarize()
        with open(summary_file, 'w') as f:
            print("P1","Hybrid","P2","Zscore_Mean","Zscore_StdDev",
                  "Pvalue_Mean","Pvalue_StdDev","Gamma_Mean","Gamma_StdDev", sep='\t', file=f)
            for k in summ.keys():
                for i in list(k):
                    print(i,"\t", sep='', end='', file=f)
                print(summ[k]['Zscore'][0], summ[k]['Zscore'][1], summ[k]['Pvalue'][0],
                      summ[k]['Pvalue'][1], summ[k]['Gamma'][0], summ[k]['Gamma'][1],
                      sep='\t', file=f)
