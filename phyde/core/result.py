""" Class for processing HyDe output. """

from __future__ import print_function
import numpy as np

class HydeResult:
    """
    A class for reading in and working with results from a HyDe analysis.
    It is mostly a simple container for storing triples as tuples and using
    them as keys in a dictionary with values that are also dictionaries
    containing the results from the HyDe analysis.

    :param str infile: name of results file.

    Example:

    .. code:: py

        import phyde as hd
        res = hd.HydeResult("hyde-out.txt")
    """
    def __init__(self, infile):
        """
        HydeResult constructor.
        """
        self.infile = infile
        self.res = {}
        self.triples = []
        self._read_hyde_results(infile)

    def __call__(self, attr, p1, hyb, p2):
        """
        A callable method (the object can be called as a function) for
        accessing information in a ``phyde.HydeResult`` object.

        :param str attr: name of hypothesis test attribute to plot (e.g., "Gamma", "Zscore", "Pvalue", etc.)
        :param str p1: parent one.
        :param str hyb: putative hybrid.
        :param str p2: parent two.

        Example:

        .. code:: py

            import phyde as hd
            res = hd.HydeResult("hyde-out.txt")
            res("Zscore", "sp1", "sp2", "sp3")
        """
        return self.res[(p1, hyb, p2)][attr]

    def _read_hyde_results(self, file):
        # Fxn for reading in results file from a hyde_cpp analysis.
        with open(file) as f:
            results_read_in = f.read()
            if results_read_in[-1] == "\n":
                results = results_read_in[:-1].splitlines()[1:]
            else:
                results = results_read_in.splitlines()[1:]
            for r in results:
                rs    = r.split('\t')
                tripl = (rs[0], rs[1], rs[2])
                if tripl not in self.triples:
                    self.triples.append(tripl)
                    self.res[tripl] = self._hyde_info(rs[3:])
                else:
                    print("\nERROR:")
                    print("  The triple ", tripl, " was tested more than once.\n")

    def _hyde_info(self, b):
        # Store information for HyDe hypothesis test.
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

    def abba_baba(self, p1, hyb, p2):
        """
        Calculate Patterson's D-Statistic for the triple (p1, hyb, p2).

        :param str p1: parent 1.
        :param str hyb: putative hybrid.
        :param str p2: parent 2.

        Example:

        .. code:: py

            import phyde as hd
            res = hd.HydeResult("hyde-out.txt")
            res.abba_baba("sp1", "sp2", "sp3")
        """
        return ((np.array(self.res[(p1,hyb,p2)]['ABBA']) - np.array(self.res[(p1,hyb,p2)]['ABAB'])) /
                (np.array(self.res[(p1,hyb,p2)]['ABBA']) + np.array(self.res[(p1,hyb,p2)]['ABAB'])))
