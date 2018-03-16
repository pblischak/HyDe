# cython: profile = True
# cython: lineprofile = True
# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True
# cython: embedsignature = True
# distutils: language = c++

from __future__ import print_function
import numpy as np
cimport numpy as np
cimport cython
import sys

from libc.math cimport fabs, sqrt, pow, exp
from libcpp.vector cimport vector

# types and typedefs for DNA bases and indices
ctypedef np.uint8_t DNA_t
DNA = np.uint8
ctypedef np.uint64_t INDEX_t
INDEX = np.uint64

# Hardcoded switch to test resolving missing and/or ambiguous bases
_ignore_amb_sites = False
if _ignore_amb_sites:
    print("\nWARNING: Ignoring sites with missing/ambiguous bases.")

cdef dict _BASE_TO_UINT8 = {
    "A": 0,  "a": 0,
    "G": 1,  "g": 1,
    "C": 2,  "c": 2,
    "T": 3,  "t": 3,
    "U": 3,  "u": 3,
    "M": 5,  "m": 5,
    "R": 6,  "r": 6,
    "W": 7,  "w": 7,
    "S": 8,  "s": 8,
    "Y": 9,  "y": 9,
    "K": 10, "k": 10,
    "B": 11, "b": 11,
    "D": 12, "d": 12,
    "H": 13, "h": 13,
    "V": 14, "v": 14,
    "N": 15, "n": 15, "?": 15,
    "-": 4
}

cdef dict _BASE_LOOKUP = {
    0: [0], 1: [1], 2: [2], 3: [3], 4: [4], 5: [0,2],
    6: [0,1], 7: [0,3], 8: [1,2], 9: [2,3], 10: [1,3],
    11: [1,2,3], 12: [0,1,3], 13: [0,2,3], 14: [0,1,2],
    15: [1,2,3,4]
}


cdef vector[vector[int]] _baseLookup = [ #/* {A,G,C,T} == {0,1,2,3} */
    [0],
    [1],
    [2],
    [3],
    [4],          #/* - == ignore */
    [0, 2],       #/* M == A or C */
    [0, 1],       #/* R == A or G */
    [0, 3],       #/* W == A ot T */
    [1, 2],       #/* S == G or C */
    [2, 3],       #/* Y == C or T */
    [1, 3],       #/* K == G or T */
    [1, 2, 3],    #/* B == C or G or T */
    [0, 1, 3],    #/* D == A or G or T */
    [0, 2, 3],    #/* H == A or C or T */
    [0, 1, 2],    #/* V == A or G or C */
    [0, 1, 2, 3]  #/* N == A or G or C or T */
]

cdef class HydeData(object):
    """
    Class for storing (1) a matrix of DNA bases as unsigned, 8-bit
    integers and (2) mapping of individuals to taxa.

    :param str infile: the name of the input DNA sequence data file.
    :param str mapfile: the name of the individual mapping file.
    :param str outgroup: the name of the outgroup.
    :param int nind: the number of individuals.
    :param int ntaxa: the number of populations/taxa.
    :param int nsites: the number of sites.
    :param bool quiet: suppress printing output.

    Example:

    .. code:: py

        import phyde as hd
        data = hd.HydeData("infile.txt", "mapfile.txt", "outgroup", 100, 6, 10000)
    """
    cdef DNA_t[:, ::1] dnaMat
    cdef INDEX_t[::1] outIndex # outgroup sequence indices
    cdef int nind, nsites
    cdef dict taxonMap
    cdef dict taxonMap_cp
    cdef double counts[16][16]
    cdef double site_pattern_probs[15]
    cdef double ind_nucl_probs[4][15]
    cdef str outgroup
    cdef bint quiet

    def __init__(self, infile=None, mapfile=None, str outgroup=None,
                 int nind=-1, int ntaxa=-1, int nsites=-1, bint quiet=False):
        """
        HydeData class constructor.
        """
        self.nind = nind
        self.nsites = nsites
        self.dnaMat = np.zeros((self.nind, self.nsites), dtype=DNA)
        self.taxonMap = {}
        self.taxonMap_cp = {}
        self.outgroup = outgroup
        self.quiet = quiet
        if not self.quiet:
            print("\nReading input file",end='')
        self._read_infile(infile)
        if not self.quiet:
            print("Done.")
        if not self.quiet:
            print("Reading map file  ",end='')
        self._read_mapfile(mapfile)
        if not self.quiet:
            print("Done.")
        self.outIndex = np.array([i[0] for i in self.taxonMap[outgroup]], dtype=INDEX)

    def _read_infile(self, infile):
        counter = 0
        first_line = 1
        valid_header = 1
        nind_not_int = 0
        with open(infile) as f:
            for line in f:
                if first_line == 1:
                    first_line = 0
                    nind_phylip = line.split()[0]
                    nsites_phylip = line.split()[1]

                    try:
                        int(nind_phylip)
                    except ValueError:
                        nind_not_int = 1

                    if nind_not_int:
                        pass
                    else:
                        continue
                        #print("\nERROR:")
                        #print("  First line of data file does not contain the correct header information.")
                        #print("  User input:")
                        #if len(nsites_phylip) > 50:
                        #    print("\n    ", nind_phylip, " ", nsites_phylip[0:49], "...", sep='')
                        #else:
                        #    print(line)
                        #sys.exit(-1)

                #try:
                #    bases = line.split()[1]
                #except IndexError:
                #    break
                if len(line) > self.nsites:
                    bases = line.split()[1]
                else:
                    break
                if len(bases) != self.nsites:
                    print("\nERROR:")
                    print("  Number of sites specified (", self.nsites, ") is not equal", sep='')
                    print("  to the number of sites in the data file (", len(bases), ").\n", sep='')
                    sys.exit(-1)
                self._convert(counter, bases)
                if not self.quiet:
                    print(".", end='')
                    sys.stdout.flush()
                #print("**",counter,"**",sep='')
                counter += 1
                if counter > self.nind:
                    print("\nERROR:")
                    print("  Number of individuals specified (", self.nind, ") is not equal", sep='')
                    print("  to the number of individuals in the data file (>= ", counter, ").\n", sep='')
                    sys.exit(-1)

    def _read_mapfile(self, mapfile):
        with open(mapfile) as f:
            lines = f.read().splitlines()
            taxa = []
            for i, l in enumerate(lines):
                if l.split()[1] not in taxa:
                    taxa.append(l.split()[1])
                    self.taxonMap[l.split()[1]] = []
                    self.taxonMap_cp[l.split()[1]] = []
                self.taxonMap[l.split()[1]].append((i, l.split()[0]))
                self.taxonMap_cp[l.split()[1]].append((i, l.split()[0]))
                if not self.quiet:
                    print(".", end='')

    def resetOutgroup(self, newOut):
        """
        Reset outgroup population.

        :param str newOut: Name of the new outgroup.
        """
        self.outgroup = newOut
        self.outIndex = np.array([i[0] for i in self.taxonMap[newOut]], dtype=INDEX)

    cdef void _convert(self, int row, str d):
        cdef unsigned long s
        for s in range(len(d)):
            self.dnaMat[row,s] = _BASE_TO_UINT8[d[s]]

    cpdef dict test_triple(self, str p1, str hyb, str p2):
        """
        Main method for testing a hypothesis on a specified triple.
        ((P1,Hyb),P2)::math:`\gamma` and (P1,(Hyb,P2))::math:`1-\gamma`.
        It is a wrapper for the C++ based methods `_test_triple_c()`,
        which is not accessible from Python.

        :param str p1:  parent one.
        :param str hyb: the putative hybrid.
        :param str p2:  parent two.
        :rtype: dict

        Returns a dictionary with the values for the Z-score, P-value, estimate
        of gamma, and all of the site pattern counts.

        Example:

        .. code:: py

          import phyde as hd
          data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
          res = data.test_triple("sp1", "sp2", "sp3")
        """
        cdef np.ndarray[INDEX_t, ndim=1] p1_rows  = np.array([i[0] for i in self.taxonMap[p1]], dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] hyb_rows = np.array([i[0] for i in self.taxonMap[hyb]],  dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] p2_rows  = np.array([i[0] for i in self.taxonMap[p2]],  dtype=INDEX)
        cdef dict res
        res = self._test_triple_c(p1_rows, hyb_rows, p2_rows)
        return res

    cpdef dict test_individuals(self, str p1, str hyb, str p2):
        """
        Test all individuals in a given putative hybrid lineage (``hyb``).

        :param str p1:  parent one.
        :param str hyb: the putative hybrid.
        :param str p2:  parent two.
        :rtype: dict

        Example:

        .. code:: py

          import phyde as hd
          data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
          res = data.test_individuals("sp1", "sp2", "sp3")
        """
        cdef:
            np.ndarray[INDEX_t, ndim=1] p1_rows  = np.array([i[0] for i in self.taxonMap[p1]], dtype=INDEX)
            np.ndarray[INDEX_t, ndim=1] hyb_rows = np.array([i[0] for i in self.taxonMap[hyb]],  dtype=INDEX)
            np.ndarray[INDEX_t, ndim=1] p2_rows  = np.array([i[0] for i in self.taxonMap[p2]],  dtype=INDEX)
            np.ndarray[INDEX_t, ndim=1] curr_ind = np.array([0], dtype=INDEX)
            int n_hyb = hyb_rows.shape[0], t
            dict res = {}

        for t in range(n_hyb):
            curr_ind[0] = hyb_rows[t]
            res[self.taxonMap[hyb][t][1]] = {}
            res[self.taxonMap[hyb][t][1]] = self._test_triple_c(p1_rows, curr_ind, p2_rows)
            #tmp_res = self._test_triple_c(p1_rows, hyb_rows[t], p2_rows)
            #res[self.taxonMap[hyb][t]] = tmp_res
        return res

    cpdef dict bootstrap_triple(self, str p1, str hyb, str p2, int reps=100):
        """
        Performs bootstrap resampling of individuals within the specified hybrid
        population.

        :param str p1:  parent one.
        :param str hyb: the putative hybrid.
        :param str p2:  parent two.
        :param int reps: number of boostrap replicates (default=100).
        :rtype: dict

        Example:

        .. code:: py

          import phyde as hd
          data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
          res = data.bootstrap_triple("sp1", "sp2", "sp3")
        """
        cdef np.ndarray[INDEX_t, ndim=1] p1_rows  = np.array([i[0] for i in self.taxonMap[p1]], dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] hyb_rows = np.array([i[0] for i in self.taxonMap[hyb]],  dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] p2_rows  = np.array([i[0] for i in self.taxonMap[p2]],  dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] hyb_resampled = hyb_rows
        cdef dict res = {}
        cdef int r
        for r in range(reps):
            hyb_resampled = np.random.choice(hyb_rows, hyb_rows.shape[0], replace=True)
            res[r+1] = {}
            res[r+1] = self._test_triple_c(p1_rows, hyb_resampled, p2_rows)
        return res

    cpdef list list_triples(self):
        """
        List all possible triples based on the populations in the given
        map file.

        :rtype: list

        Returns a list of 3-tuples: ``(p1, hyb, p2)``
        """
        res = []
        ks = list(self.taxonMap)
        for a in range(len(ks)):
            if ks[a] == self.outgroup:
                del ks[a]
        for i in range(len(ks)-2):
            for j in range(i+1, len(ks)-1):
                for k in range(j+1, len(ks)):
                    res.append((ks[i],ks[j], ks[k]))
                    res.append((ks[i],ks[k], ks[j]))
                    res.append((ks[j],ks[i], ks[k]))
        return res

    cdef dict _test_triple_c(self, np.ndarray[INDEX_t, ndim=1] p1,
                             np.ndarray[INDEX_t, ndim=1] hyb,
                             np.ndarray[INDEX_t, ndim=1] p2):
        """

        """
        cdef:
            int i, j, k, l, c1, c2
            int n_out = self.outIndex.shape[0], n_p1 = p1.shape[0]
            int n_hyb = hyb.shape[0], n_p2 = p2.shape[0]
            double num_obs = 0.0, avg_obs = 0.0, z_val = 0.0, p_val = 0.0
            double _c_num = 0.0, _c_denom = 0.0, _c = 0.0, gamma = 0.0

        for c1 in range(16):
            for c2 in range(16):
                self.counts[c1][c2] = 0.0

        for i in range(n_out):
            for j in range(n_p1):
                for k in range(n_hyb):
                    for l in range(n_p2):
                        num_obs += self._get_counts(self.outIndex[i], p1[j], hyb[k], p2[l])

        avg_obs  = num_obs / (self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
        z_val    = self._calc_gh(num_obs, avg_obs)
        p_val    = self._calc_p_value(z_val)
        _c_num   = avg_obs * (self.site_pattern_probs[8] - self.site_pattern_probs[6])
        _c_denom = avg_obs * (self.site_pattern_probs[3] - self.site_pattern_probs[6])
        _c       = _c_num / _c_denom
        gamma = _c / (1 + _c)

        return {
            "Zscore": z_val,
            "Pvalue": p_val,
            "Gamma" : gamma,
            "AAAA"  : self.site_pattern_probs[0],
            "AAAB"  : self.site_pattern_probs[1],
            "AABA"  : self.site_pattern_probs[2],
            "AABB"  : self.site_pattern_probs[3],
            "AABC"  : self.site_pattern_probs[4],
            "ABAA"  : self.site_pattern_probs[5],
            "ABAB"  : self.site_pattern_probs[6],
            "ABAC"  : self.site_pattern_probs[7],
            "ABBA"  : self.site_pattern_probs[8],
            "BAAA"  : self.site_pattern_probs[9],
            "ABBC"  : self.site_pattern_probs[10],
            "CABC"  : self.site_pattern_probs[11],
            "BACA"  : self.site_pattern_probs[12],
            "BCAA"  : self.site_pattern_probs[13],
            "ABCD"  : self.site_pattern_probs[14]
        }

    @cython.nonecheck(False)
    cdef double _get_counts(self, int out, int p1, int hyb, int p2):
        """

        """
        cdef:
            double nn = 0.0, resolved = 0.0
            DNA_t i, j, k, l
            int s, sites = self.dnaMat.shape[1]
        for s in range(sites):
            i = self.dnaMat[out,s]
            j = self.dnaMat[p1,s]
            k = self.dnaMat[hyb,s]
            l = self.dnaMat[p2,s]
            if i < 4 and j < 4 and k < 4 and l < 4:
                nn += 1.0
                self.counts[i * 4 + j][k * 4 + l] += 1.0
            elif not _ignore_amb_sites:
                resolved = self._resolve_ambiguity_cpp(i, j, k, l)
                nn += resolved
            else:
              pass

        return nn

    cdef double _calc_gh(self, double nobs, double avobs):
        """

        """
        # AAAA
        self.site_pattern_probs[0] = self.counts[0][0] + self.counts[5][5] + self.counts[10][10] + self.counts[15][15]
        # AAAB
        self.site_pattern_probs[1] = (self.counts[0][1] + self.counts[0][2] + self.counts[0][3] + self.counts[5][4]
                                    + self.counts[5][6] + self.counts[5][7] + self.counts[10][8] + self.counts[10][9]
                                    + self.counts[10][11] + self.counts[15][12] + self.counts[15][13] + self.counts[15][14])
        # AABA
        self.site_pattern_probs[2] = (self.counts[0][4] + self.counts[0][8] + self.counts[0][12] + self.counts[5][1]
                                    + self.counts[5][9] + self.counts[5][13] + self.counts[10][2] + self.counts[10][6]
                                    + self.counts[10][14] + self.counts[15][3] + self.counts[15][7] + self.counts[15][11])
        # AABB
        self.site_pattern_probs[3] = (self.counts[0][5] + self.counts[0][10] + self.counts[0][15] + self.counts[5][0]
                                    + self.counts[5][10] + self.counts[5][15] + self.counts[10][0] + self.counts[10][5]
                                    + self.counts[10][15] + self.counts[15][0] + self.counts[15][5] + self.counts[15][10])
        # AABC
        self.site_pattern_probs[4] = (self.counts[0][6] + self.counts[0][7] + self.counts[0][9] + self.counts[0][11]
                                    + self.counts[0][13] + self.counts[0][14] + self.counts[5][2] + self.counts[5][3]
                                    + self.counts[5][8] + self.counts[5][11] + self.counts[5][12] + self.counts[5][14]
                                    + self.counts[10][1] + self.counts[10][3] + self.counts[10][4] + self.counts[10][7]
                                    + self.counts[10][12] + self.counts[10][13] + self.counts[15][1] + self.counts[15][2]
                                    + self.counts[15][4] + self.counts[15][6] + self.counts[15][8] + self.counts[15][9])
        # ABAA
        self.site_pattern_probs[5] = (self.counts[1][0] + self.counts[2][0] + self.counts[3][0] + self.counts[4][5]
                                    + self.counts[6][5] + self.counts[7][5] + self.counts[8][10] + self.counts[9][10]
                                    + self.counts[11][10] + self.counts[12][15] + self.counts[13][15] + self.counts[14][15])
        # ABAB
        self.site_pattern_probs[6] = (self.counts[1][1] + self.counts[2][2] + self.counts[3][3] + self.counts[4][4]
                                    + self.counts[6][6] + self.counts[7][7] + self.counts[8][8] + self.counts[9][9]
                                    + self.counts[11][11] + self.counts[12][12] + self.counts[13][13] + self.counts[14][14])
        # ABAC
        self.site_pattern_probs[7] = (self.counts[1][2] + self.counts[1][3] + self.counts[2][1] + self.counts[2][3]
                                    + self.counts[3][1] + self.counts[3][2] + self.counts[4][6] + self.counts[4][7]
                                    + self.counts[6][4] + self.counts[6][7] + self.counts[7][4] + self.counts[7][6]
                                    + self.counts[8][9] + self.counts[8][11] + self.counts[9][8] + self.counts[9][11]
                                    + self.counts[11][8] + self.counts[11][9] + self.counts[12][13] + self.counts[12][14]
                                    + self.counts[13][12] + self.counts[13][14] + self.counts[14][12] + self.counts[14][13])
        # ABBA
        self.site_pattern_probs[8] = (self.counts[1][4] + self.counts[2][8] + self.counts[3][12] + self.counts[4][1]
                                    + self.counts[6][9] + self.counts[7][13] + self.counts[8][2] + self.counts[9][6]
                                    + self.counts[11][14] + self.counts[12][3] + self.counts[13][7] + self.counts[14][11])
        # BAAA
        self.site_pattern_probs[9] = (self.counts[4][0] + self.counts[8][0] + self.counts[12][0] + self.counts[1][5]
                                    + self.counts[9][5] + self.counts[13][5] + self.counts[2][10] + self.counts[6][10]
                                    + self.counts[14][10] + self.counts[3][15] + self.counts[7][15] + self.counts[11][15])
        # ABBC
        self.site_pattern_probs[10] = (self.counts[1][6] + self.counts[1][7] + self.counts[2][9] + self.counts[2][11]
                                     + self.counts[3][13] + self.counts[3][14] + self.counts[4][2] + self.counts[4][3]
                                     + self.counts[6][8] + self.counts[6][11] + self.counts[7][12] + self.counts[7][14]
                                     + self.counts[8][1] + self.counts[8][3] + self.counts[9][4] + self.counts[9][7]
                                     + self.counts[11][12] + self.counts[11][13] + self.counts[12][1] + self.counts[12][2]
                                     + self.counts[13][4] + self.counts[13][6] + self.counts[14][8] + self.counts[14][9])
        # CABC
        self.site_pattern_probs[11] = (self.counts[8][6] + self.counts[12][7] + self.counts[4][9] + self.counts[12][11]
                                     + self.counts[4][13] + self.counts[8][14] + self.counts[9][2] + self.counts[13][3]
                                     + self.counts[1][8] + self.counts[13][11] + self.counts[1][12] + self.counts[9][14]
                                     + self.counts[6][1] + self.counts[14][3] + self.counts[2][4] + self.counts[14][7]
                                     + self.counts[2][12] + self.counts[6][13] + self.counts[7][1] + self.counts[11][2]
                                     + self.counts[3][4] + self.counts[11][6] + self.counts[3][8] + self.counts[7][9])
        # BACA
        self.site_pattern_probs[12] = (self.counts[4][8] + self.counts[4][12] + self.counts[8][4] + self.counts[8][12]
                                     + self.counts[12][4] + self.counts[12][8] + self.counts[1][9] + self.counts[1][13]
                                     + self.counts[9][1] + self.counts[9][13] + self.counts[13][1] + self.counts[13][9]
                                     + self.counts[2][6] + self.counts[2][14] + self.counts[6][2] + self.counts[6][14]
                                     + self.counts[14][2] + self.counts[14][6] + self.counts[3][7] + self.counts[3][11]
                                     + self.counts[7][3] + self.counts[7][11] + self.counts[11][3] + self.counts[11][7])
        # BCAA
        self.site_pattern_probs[13] = (self.counts[6][0] + self.counts[7][0] + self.counts[9][0] + self.counts[11][0]
                                     + self.counts[13][0] + self.counts[14][0] + self.counts[2][5] + self.counts[3][5]
                                     + self.counts[8][5] + self.counts[11][5] + self.counts[12][5] + self.counts[14][5]
                                     + self.counts[1][10] + self.counts[3][10] + self.counts[4][10] + self.counts[7][10]
                                     + self.counts[12][10] + self.counts[13][10] + self.counts[1][15] + self.counts[2][15]
                                     + self.counts[4][15] + self.counts[6][15] + self.counts[8][15] + self.counts[9][15])
        # ABCD
        self.site_pattern_probs[14] = (self.counts[1][11] + self.counts[1][14] + self.counts[2][7] + self.counts[2][13]
                                     + self.counts[3][6] + self.counts[3][9] + self.counts[4][11] + self.counts[4][14]
                                     + self.counts[6][3] + self.counts[6][12] + self.counts[7][2] + self.counts[7][8]
                                     + self.counts[8][7] + self.counts[8][13] + self.counts[9][3] + self.counts[9][12]
                                     + self.counts[11][1] + self.counts[11][4] + self.counts[12][6] + self.counts[12][9]
                                     + self.counts[13][2] + self.counts[13][8] + self.counts[14][1] + self.counts[14][4])

        if (fabs((1.0 / nobs) * (self.site_pattern_probs[0] + self.site_pattern_probs[1] + self.site_pattern_probs[2] + self.site_pattern_probs[3]
                              + self.site_pattern_probs[4]  + self.site_pattern_probs[5] + self.site_pattern_probs[6] + self.site_pattern_probs[7]
                              + self.site_pattern_probs[8]  + self.site_pattern_probs[9] + self.site_pattern_probs[10] + self.site_pattern_probs[11]
                              + self.site_pattern_probs[12] + self.site_pattern_probs[13] + self.site_pattern_probs[14]) - 1.0) > 0.05):
            if not self.quiet: print("** WARNING: There was a problem counting site patterns. **")
            return -99999.9

        cdef:
            double p9 = (self.site_pattern_probs[8] + 0.05) / nobs
            double p7 = (self.site_pattern_probs[6] + 0.05) / nobs
            double p4 = (self.site_pattern_probs[3] + 0.05) / nobs
            double obs_invp1 = avobs * (p9 - p7)
            double obs_invp2 = avobs * (p4 - p7)
        if obs_invp1 == 0:
            obs_invp1 += 1.0
            obs_invp2 += 1.0
        cdef:
            double obs_var_invp1 = avobs * p9 * (1 - p9) + avobs * p7 * (1 - p7) + 2 * avobs * p9 * p7
            double obs_var_invp2 = avobs * p4 * (1 - p4) + avobs * p7 * (1 - p7) + 2 * avobs * p4 * p7
            double obs_cov_invp1_invp2 = -1 * avobs * p9 * p4 + avobs * p9 * p7 + avobs * p7 * p4 + avobs * p7 * (1 - p7)
            double ratio = obs_invp2 / obs_invp1;
            double GH_ts = ((obs_invp1) * (ratio) / sqrt(obs_var_invp1 * (pow(ratio, 2.0)) - 2.0
                            * obs_cov_invp1_invp2 * ratio + obs_var_invp2))

        cdef double temp = -99999.9
        if p7 > p9 and p7 < p4:
            return temp
        elif GH_ts > -99999.9 and GH_ts < 99999.9:
            return GH_ts
        else:
            return temp

    cdef double _calc_p_value(self, double my_z):
        """

        """
        cdef:
            double a1 =  0.254829592;
            double a2 = -0.284496736;
            double a3 =  1.421413741;
            double a4 = -1.453152027;
            double a5 =  1.061405429;
            double p  =  0.3275911;
            int sign = 1;
        if my_z < 0:
            sign = -1;
        cdef:
            double z = fabs(my_z) / sqrt(2.0);
            double t = 1.0 / (1.0 + p * z);
            double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-z * z);

        return 1.0 - (0.5 * (1.0 + sign * y))

    cdef double _resolve_ambiguity(self, int out, int p1, int hyb, int p2):
        """

        """
        cdef double denom = 0.0
        cdef unsigned i, j, k, l
        if (out >= 4 and p1 >= 4 and hyb >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            p1 >= 4 and hyb >= 4 and p2 >= 4):
            return 0.0
        elif (out == 4 or p1 == 4 or hyb ==4 or p2 == 4):
            return 0.0
        else:
            denom = len(_BASE_LOOKUP[out]) * len(_BASE_LOOKUP[p1]) * len(_BASE_LOOKUP[hyb]) * len(_BASE_LOOKUP[p2])
            for i in range(len(_BASE_LOOKUP[out])):
                for j in range(len(_BASE_LOOKUP[p1])):
                    for k in range(len(_BASE_LOOKUP[hyb])):
                        for l in range(len(_BASE_LOOKUP[p2])):
                            self.counts[_BASE_LOOKUP[out][i] * 4 + _BASE_LOOKUP[p1][j]][_BASE_LOOKUP[hyb][k] * 4 + _BASE_LOOKUP[p2][l]] += 1.0 / denom
            return 1.0

    cdef double _resolve_ambiguity_cpp(self, int out, int p1, int hyb, int p2):
        """

        """
        cdef double denom = 0.0
        cdef unsigned i, j, k, l
        if (out >= 4 and p1 >= 4 and hyb >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            p1 >= 4 and hyb >= 4 and p2 >= 4):
            return 0.0
        elif (out == 4 or p1 == 4 or hyb ==4 or p2 == 4):
            return 0.0
        else:
            denom = _baseLookup[out].size() * _baseLookup[p1].size() * _baseLookup[hyb].size() * _baseLookup[p2].size()
            for i in range(_baseLookup[out].size()):
                for j in range(_baseLookup[p1].size()):
                    for k in range(_baseLookup[hyb].size()):
                        for l in range(_baseLookup[p2].size()):
                            self.counts[_baseLookup[out][i] * 4 + _baseLookup[p1][j]][_baseLookup[hyb][k] * 4 + _baseLookup[p2][l]] += (1.0 / denom)
            return 1.0
