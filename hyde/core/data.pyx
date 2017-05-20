from __future__ import print_function, division
import numpy as np
cimport numpy as np
cimport cython

# types and typedefs for DNA bases and indices
ctypedef np.uint8_t DNA_t
DNA = np.uint8
ctypedef np.uint64_t INDEX_t
INDEX = np.uint64

_BASE_TO_UINT8 = {
    "A": 0,  "G": 1,  "C": 2,  "T": 3,  "U": 3,
    "M": 5,  "R": 6,  "W": 7,  "S": 8,  "Y": 9,
    "K": 10, "B": 11, "D": 12, "H": 13, "V": 14,
    "N": 15, "?": 15, "-": 4
}

_BASE_LOOKUP = {

}

cdef class HydeData:
    """
    Class for storing (1) a matrix of DNA bases as unsigned, 8-bit
    integers; (2) mapping of individuals to taxa (optional); and (3)
    partition information for multilocus data (optional).
    """
    cdef public np.ndarray dnaMat
    cdef public np.ndarray outIndex # outgroup sequence indices
    cdef public dict taxonMap
    cdef double counts[16][16]
    cdef double probs[4][15]

    def __init__(self, infile, mapfile, outgroup, int nind, int nsites, int ntaxa, partition=None):
        """
        Constructor:
            Read infile with DNA characters. Parse mapfile and partition
            file if they are given. Convert DNA bases to integer codes.
        """
        self.dnaMat = np.zeros((nind, nsites), dtype=DNA)
        self.taxonMap = {}
        if partition is not None:
            self.partitions = {}
            self._read_partfile(partition)
        else:
            pass
        self._read_infile(infile)
        self._read_mapfile(mapfile)
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
        with open(infile) as f:
            lines = f.read().splitlines()
            dna = [lines[i].split()[1] for i in range(len(lines))]
            self._convert(dna)

    def _read_mapfile(self, mapfile):
        with open(mapfile) as f:
            lines = f.read().splitlines()
            taxa = []
            for i, l in enumerate(lines):
                if l.split()[1] not in taxa:
                    taxa.append(l.split()[1])
                    self.taxonMap[l.split()[1]] = []
                self.taxonMap[l.split()[1]].append((i, l.split()[0]))

    def _read_partfile(self, partfile):
        with open(partfile) as f:
            lines = f.read().splitlines()
            for l in lines:
                entry = l.split("=")
                start_stop = entry[1].split("-")
                self.partitions[entry[0]] = start_stop

    def _convert(self, d):
        for i in range(len(d)):
            for s in range(len(d[i])):
                self.dnaMat[i,s] = _BASE_TO_UINT8[d[i][s]]

    def test_triple(self, p1, hyb, p2, fast=False):
        """
        Main method for testing a hypothesis on a specified triple.
        ((P1,Hyb),P2):gamma and (P1,(Hyb,P2)):1-gamma.

        It is a wrapper for the C based methods `_test_triple_c()` and
        `_test_triple_c_fast()`, which are not accessible from Python.

        Arguments
        =========

            - P1: parent one
            - Hyb: the putative hybrid
            - P2: parent two
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
        if fast:
            res = self._test_triple_c_fast(p1_rows, hyb_rows, p2_rows)
        elif not fast:
            res = self._test_triple_c(p1_rows, hyb_rows, p2_rows)
        else:
            pass

    #def test_triple_py(self, p1, hyb, p2):
    #    for i in self.outIndex:
    #        for j in [i[0] for i in self.taxonMap[p1]]:
    #            for k in [i[0] for i in self.taxonMap[hyb]]:
    #                for l in [i[0] for i in self.taxonMap[p2]]:
    #                    for s in range(self.dnaMat.shape[1]):
    #                        i * j * k * l * s

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _test_triple_c(self, np.ndarray[INDEX_t, ndim=1] p1,
                             np.ndarray[INDEX_t, ndim=1] hyb,
                             np.ndarray[INDEX_t, ndim=1] p2):
        """

        """
        cdef int i, j, k, l, s, sites = self.dnaMat.shape[1]
        for i in self.outIndex:
            for j in p1:
                for k in hyb:
                    for l in p2:
                        for s in range(sites):
                            i * j * k * l * s

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _test_triple_c_fast(self,  np.ndarray[INDEX_t, ndim=1] p1,
                                  np.ndarray[INDEX_t, ndim=1] hyb,
                                  np.ndarray[INDEX_t, ndim=1] p2):
        """

        """
        cdef INDEX_t i, j, k, l
        cdef int s, sites = self.dnaMat.shape[1]
        for i in self.outIndex:
            for j in p1:
                for k in hyb:
                    for l in p2:
                        for s in range(sites):
                            i * j * k * l * s
