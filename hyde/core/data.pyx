# cython: profile=True
# cython: line_profile:True

from __future__ import print_function, division
import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport fabs, sqrt, pow, exp

# types and typedefs for DNA bases and indices
ctypedef np.uint8_t DNA_t
DNA = np.uint8
ctypedef np.uint64_t INDEX_t
INDEX = np.uint64

cdef dict _BASE_TO_UINT8 = {
    "A": 0,  "G": 1,  "C": 2,  "T": 3,  "U": 3,
    "M": 5,  "R": 6,  "W": 7,  "S": 8,  "Y": 9,
    "K": 10, "B": 11, "D": 12, "H": 13, "V": 14,
    "N": 15, "?": 15, "-": 4
}

cdef dict _BASE_LOOKUP = {
    0: [0], 1: [1], 2: [2], 3: [3], 4: [4], 5: [0,2],
    6: [0,1], 7: [0,3], 8: [1,2], 9: [2,3], 10: [1,3],
    11: [1,2,3], 12: [0,1,3], 13: [0,2,3], 14: [0,1,2],
    15: [1,2,3,4]
}

cdef class HydeData:
    """
    Class for storing (1) a matrix of DNA bases as unsigned, 8-bit
    integers; (2) mapping of individuals to taxa (optional); and (3)
    partition information for multilocus data (optional).
    """
    cdef np.ndarray dnaMat
    cdef np.ndarray outIndex # outgroup sequence indices
    cdef dict taxonMap
    cdef double counts[16][16]
    cdef double site_pattern_probs[15]
    cdef double ind_nucl_probs[4][15]
    cdef bytes outgroup

    def __init__(self, infile, mapfile, outgroup, int nind, int nsites, int ntaxa, partition=None):
        """
        Constructor:
            Read infile with DNA characters. Parse mapfile and partition
            file if they are given. Convert DNA bases to integer codes.

        Example:
            >>> data = HydeData("infile.txt", "mapfile.txt", "outgroup", 100, 10000, 6)
            >>> data.dnaMat

        Arguments
        =========

            - infile <str>: name of the input file

            - mapfile <str>: name of the individual-to-OTU map

            - outgroup <str>: name of the outgroup (can be
              reset using the `resetOutgroup()` method)

            - nind <int>: total number of individuals sampled

            - nsites <int>: number of sites

            - ntaxa <int>: number of taxa
        """
        self.dnaMat = np.zeros((nind, nsites), dtype=DNA)
        self.taxonMap = {}
        self.outgroup = outgroup
        if partition is not None:
            self.partitions = {}
            self._read_partfile(partition)
        else:
            pass
        print("Reading input file",end='')
        self._read_infile(infile)
        print("Done.")
        print("Reading map file",end='')
        self._read_mapfile(mapfile)
        print("Done.")
        self.outIndex = np.array([i[0] for i in self.taxonMap[outgroup]], dtype=INDEX)

    def __call__(self, p1, hyb, p2, fast=False):
        cdef np.ndarray[INDEX_t, ndim=1] p1_rows  = np.array([i[0] for i in self.taxonMap[p1]], dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] hyb_rows = np.array([i[0] for i in self.taxonMap[hyb]],  dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] p2_rows  = np.array([i[0] for i in self.taxonMap[p2]],  dtype=INDEX)
        if fast:
            res = self._test_triple_c_fast(p1_rows, hyb_rows, p2_rows)
        elif not fast:
            res = self._test_triple_c(p1_rows, hyb_rows, p2_rows)
        else:
            pass

    def _read_infile(self, infile):
        counter = 0
        with open(infile) as f:
            for line in f:
                try:
                    bases = line.split()[1]
                except IndexError:
                    break
                self._convert(counter, bases)
                print(".", end='')
                counter += 1

    def _read_mapfile(self, mapfile):
        with open(mapfile) as f:
            lines = f.read().splitlines()
            taxa = []
            for i, l in enumerate(lines):
                if l.split()[1] not in taxa:
                    taxa.append(l.split()[1])
                    self.taxonMap[l.split()[1]] = []
                self.taxonMap[l.split()[1]].append((i, l.split()[0]))
                print(".", end='')

    def _read_partfile(self, partfile):
        with open(partfile) as f:
            lines = f.read().splitlines()
            for l in lines:
                entry = l.split("=")
                start_stop = entry[1].split("-")
                self.partitions[entry[0]] = start_stop

    def resetOutgroup(self, newOut):
        self.outgroup = newOut
        self.outIndex = np.array([i[0] for i in self.taxonMap[newOut]], dtype=INDEX)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _convert(self, int row, bytes d):
        cdef unsigned long s
        for s in range(len(d)):
            self.dnaMat[row,s] = _BASE_TO_UINT8[d[s]]

    cpdef dict test_triple(self, bytes p1, bytes hyb, bytes p2, bint fast=False):
        """
        Main method for testing a hypothesis on a specified triple.
        ((P1,Hyb),P2):gamma and (P1,(Hyb,P2)):1-gamma.

        It is a wrapper for the C based methods `_test_triple_c()` and
        `_test_triple_c_fast()`, which are not accessible from Python.

        Arguments
        =========

            - p1 <str>: parent one
            - hyb <str>: the putative hybrid
            - p2 <str>: parent two
            - fast: logical variable to indicate whether a fast
              analysis is run or not (default=False).

        A note on "fast" analyses
        -------------------------

        Fast analyses are fast because they precalculate the frequency of each site pattern
        in each taxon first and then use them in the calculation of the overall site patterns.
        The alternative is to run all permutations of individuals within taxa, which is more
        thorough but is also much more time consuming.
        """
        cdef np.ndarray[INDEX_t, ndim=1] p1_rows  = np.array([i[0] for i in self.taxonMap[p1]], dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] hyb_rows = np.array([i[0] for i in self.taxonMap[hyb]],  dtype=INDEX)
        cdef np.ndarray[INDEX_t, ndim=1] p2_rows  = np.array([i[0] for i in self.taxonMap[p2]],  dtype=INDEX)
        cdef dict res = {}
        if fast:
            res = self._test_triple_c_fast(p1_rows, hyb_rows, p2_rows)
        elif not fast:
            res = self._test_triple_c(p1_rows, hyb_rows, p2_rows)
        else:
            pass

        return res
    #def test_triple_py(self, p1, hyb, p2):
    #    for i in self.outIndex:
    #        for j in [i[0] for i in self.taxonMap[p1]]:
    #            for k in [i[0] for i in self.taxonMap[hyb]]:
    #                for l in [i[0] for i in self.taxonMap[p2]]:
    #                    for s in range(self.dnaMat.shape[1]):
    #                        i * j * k * l * s

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef dict _test_triple_c(self, np.ndarray[INDEX_t, ndim=1] p1,
                             np.ndarray[INDEX_t, ndim=1] hyb,
                             np.ndarray[INDEX_t, ndim=1] p2):
        """

        """
        cdef int i, j, k, l, s, sites = self.dnaMat.shape[1]
        cdef double num_obs = 0.0, avg_obs = 0.0, z_val = 0.0, p_val = 0.0, gamma = 0.0
        for i in self.outIndex:
            for j in p1:
                for k in hyb:
                    for l in p2:
                        num_obs += self._get_counts(i, j, k, l)

        avg_obs = num_obs / (self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
        z_val   = self._calc_gh(num_obs, avg_obs)
        p_val   = self._calc_p_value(z_val)

        return {
            "GH_statistic": z_val,
            "P_value": p_val,
            "Gamma": self.site_pattern_probs[3] - self.site_pattern_probs[6]
        }


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef dict _test_triple_c_fast(self,  np.ndarray[INDEX_t, ndim=1] p1,
                                  np.ndarray[INDEX_t, ndim=1] hyb,
                                  np.ndarray[INDEX_t, ndim=1] p2):
        """

        """
        cdef int i, j, k, l
        cdef int s, sites = self.dnaMat.shape[1]
        for i in self.outIndex:
            for j in p1:
                for k in hyb:
                    for l in p2:
                        for s in range(sites):
                            i * j * k * l * s

    cdef double _get_counts(self, int i, int j, int k, int l):
        """

        """
        cdef double nn = 0.0, resolved = 0.0
        cdef int s, sites = self.dnaMat.shape[1]
        for s in range(sites):
            if self.dnaMat[i][s] < 4 and self.dnaMat[j][s] < 4 and self.dnaMat[k][s] < 4 and self.dnaMat[l][s] < 4:
                nn += 1.0
                self.counts[self.dnaMat[i][s] * 4 + self.dnaMat[j][s]][self.dnaMat[k][s] * 4 + self.dnaMat[l][s]] += 1.0
            else:
                resolved = self._resolve_ambiguity(self.dnaMat[i][s], self.dnaMat[j][s],
                                              self.dnaMat[k][s], self.dnaMat[l][s])
                nn += resolved

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
                              + self.site_pattern_probs[4] + self.site_pattern_probs[5] + self.site_pattern_probs[6] + self.site_pattern_probs[7]
                              + self.site_pattern_probs[8] + self.site_pattern_probs[9] + self.site_pattern_probs[10] + self.site_pattern_probs[11]
                              + self.site_pattern_probs[12] + self.site_pattern_probs[13] + self.site_pattern_probs[14]) - 1.0) > 0.05):
            print("** WARNING: There was a problem counting site patterns. **")
            return -99999.9

        cdef double p9 = (self.site_pattern_probs[8] + 0.05) / nobs
        cdef double p7 = (self.site_pattern_probs[6] + 0.05) / nobs
        cdef double p4 = (self.site_pattern_probs[3] + 0.05) / nobs
        cdef double obs_invp1 = avobs * (p9 - p7)
        cdef double obs_invp2 = avobs * (p4 - p7)
        if obs_invp1 == 0:
            obs_invp1 += 1.0
            obs_invp2 += 1.0
        cdef double obs_var_invp1 = avobs * p9 * (1 - p9) + avobs * p7 * (1 - p7) + 2 * avobs * p9 * p7
        cdef double obs_var_invp2 = avobs * p4 * (1 - p4) + avobs * p7 * (1 - p7) + 2 * avobs * p4 * p7
        cdef double obs_cov_invp1_invp2 = -1 * avobs * p9 * p4 + avobs * p9 * p7 + avobs * p7 * p4 + avobs * p7 * (1 - p7)
        cdef double ratio = obs_invp2 / obs_invp1;
        cdef double GH_ts = ((obs_invp1) * (ratio) / sqrt(obs_var_invp1 * (pow(ratio, 2.0)) - 2.0
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
        cdef double a1 =  0.254829592;
        cdef double a2 = -0.284496736;
        cdef double a3 =  1.421413741;
        cdef double a4 = -1.453152027;
        cdef double a5 =  1.061405429;
        cdef double p  =  0.3275911;
        cdef int sign = 1;
        if my_z < 0:
            sign = -1;
        cdef double z = fabs(my_z) / sqrt(2.0);
        cdef double t = 1.0 / (1.0 + p * z);
        cdef double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-z * z);

        return 1.0 - (0.5 * (1.0 + sign * y))

    cdef double _resolve_ambiguity(self, int out, int p1, int hyb, int p2):
        """

        """
        cdef double denom = 0.0
        cdef int i, j, k, l
        if (out >= 4 and p1 >= 4 and hyb >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            out >= 4 and hyb >= 4 and p2 >= 4 or
            p1 >= 4 and hyb >= 4 and p2 >= 4):
            return 0.0
        elif (out == 4 or p1 == 4 or hyb ==4 or p2 == 4):
            return 0.0
        else:
            denom = len(_BASE_LOOKUP[out]) * len(_BASE_LOOKUP[p1]) * len(_BASE_LOOKUP[hyb]) * len(_BASE_LOOKUP[p1])
            for i in range(len(_BASE_LOOKUP[out])):
                for j in range(len(_BASE_LOOKUP[p1])):
                    for k in range(len(_BASE_LOOKUP[hyb])):
                        for l in range(len(_BASE_LOOKUP[p2])):
                            self.counts[_BASE_LOOKUP[out][i] * 4 + _BASE_LOOKUP[p1][j]][_BASE_LOOKUP[hyb][k] * 4 + _BASE_LOOKUP[p2][l]] += 1.0 / denom
            return 1.0
