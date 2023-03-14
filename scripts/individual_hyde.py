#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# individual_hyde.py
# Written by PD Blischak

"""
<< individual_hyde.py >>

Test all individuals within a putative hybrid lineage for a specified triple.

Arguments
---------

For more details on script arguments, type: individual_hyde.py -h

    - infile         <string> : name of the DNA sequence data file.
    - mapfile        <string> : name of the taxon map file.
    - outgroup       <string> : name of the outgroup.
    - triples        <string> : name of the file containing triples for testing.
    - nind              <int> : number of sampled individuals.
    - nsites            <int> : number of sampled sites.
    - ntaxa             <int> : number of sampled taxa (populations, OTUs, etc.).
    - prefix         <string> : name added to the beginning of output file.
    - quiet            <flag> : suppress printing to stdout.
    - ignore_amb_sites <flag> : ignore missing/ambiguous sites.

Output
------

    Writes a file ('hyde-ind.txt') with the results of the hybridization
    detection analysis for each individual tested against the putative parents
    specified in the triples file. P1 and P2 will be the name of the parents
    but the Hybrid column will contain the names of all of the individuals that
    were tested.
"""

import phyde as hd
import argparse
import sys
import os


def parse_triples(triples_file):
    """
    Parse a three column table or previous results file to get the names
    of the taxa that are going to be tested for hybridization.

    Returns a list of three-tuples of the form (P1, Hybrid, P2) for all triples
    that are to be tested.
    """
    triples = []
    with open(triples_file) as f:
        lines = f.read().splitlines()
        # remove header information if reading in a previous results file
        if (
            lines[0].split()[0] == "P1"
            and lines[0].split()[1] == "Hybrid"
            and lines[0].split()[2] == "P2"
        ):
            lines = lines[1:]
        # catch the case where the last line in the file is blank
        if len(lines[-1]) == 0:
            triples = [
                (l.split()[0], l.split()[1], l.split()[2]) for l in lines[:-1]
            ]
        else:
            triples = [
                (l.split()[0], l.split()[1], l.split()[2]) for l in lines
            ]
    return triples


def write_ind(out, triple, outfile):
    """
    Write the current output dictionary from test_individuals()
    to the file passed as an argument to the function.
    """
    for k, v in out.items():
        print(
            triple[0],
            "\t",
            k,
            "\t",
            triple[2],
            "\t",
            sep="",
            end="",
            file=outfile,
        )
        print(v["Zscore"], "\t", sep="", end="", file=outfile)
        print(v["Pvalue"], "\t", sep="", end="", file=outfile)
        print(v["Gamma"], "\t", sep="", end="", file=outfile)
        print(v["AAAA"], "\t", sep="", end="", file=outfile)
        print(v["AAAB"], "\t", sep="", end="", file=outfile)
        print(v["AABA"], "\t", sep="", end="", file=outfile)
        print(v["AABB"], "\t", sep="", end="", file=outfile)
        print(v["AABC"], "\t", sep="", end="", file=outfile)
        print(v["ABAA"], "\t", sep="", end="", file=outfile)
        print(v["ABAB"], "\t", sep="", end="", file=outfile)
        print(v["ABAC"], "\t", sep="", end="", file=outfile)
        print(v["ABBA"], "\t", sep="", end="", file=outfile)
        print(v["BAAA"], "\t", sep="", end="", file=outfile)
        print(v["ABBC"], "\t", sep="", end="", file=outfile)
        print(v["CABC"], "\t", sep="", end="", file=outfile)
        print(v["BACA"], "\t", sep="", end="", file=outfile)
        print(v["BCAA"], "\t", sep="", end="", file=outfile)
        print(v["ABCD"], "\n", sep="", end="", file=outfile)


if __name__ == "__main__":
    """
    Runs the scripts.
    """
    # print docstring if only the name of the script is given
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Options for individual_hyde.py", add_help=True
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        "--infile",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="name of the data input file",
    )
    required.add_argument(
        "-m",
        "--map",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="map of individuals to taxa",
    )
    required.add_argument(
        "-o",
        "--outgroup",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="name of the outgroup (only one accepted)",
    )
    required.add_argument(
        "-tr",
        "--triples",
        action="store",
        type=str,
        required=True,
        metavar="\b",
        help="table of triples to be analyzed.",
    )
    required.add_argument(
        "-n",
        "--num_ind",
        action="store",
        type=int,
        required=True,
        metavar="\b",
        help="number of individuals in data matrix",
    )
    required.add_argument(
        "-t",
        "--num_taxa",
        action="store",
        type=int,
        required=True,
        metavar="\b",
        help="number of taxa (species, OTUs)",
    )
    required.add_argument(
        "-s",
        "--num_sites",
        action="store",
        type=int,
        required=True,
        metavar="\b",
        help="number of sites in the data matrix",
    )

    additional = parser.add_argument_group("additional arguments")
    additional.add_argument(
        "--prefix",
        action="store",
        type=str,
        default="hyde",
        metavar="\b",
        help="prefix appended to output files [default=hyde]",
    )
    additional.add_argument(
        "-q", "--quiet", action="store_true", help="supress printing to stdout"
    )
    additional.add_argument(
        "--ignore_amb_sites",
        action="store_true",
        help="ignore missing/ambiguous sites",
    )

    args = parser.parse_args()
    infile = args.infile
    mapfile = args.map
    outgroup = args.outgroup
    triples = parse_triples(args.triples)
    nind = args.num_ind
    ntaxa = args.num_taxa
    nsites = args.num_sites
    prefix = args.prefix
    quiet = args.quiet
    ignore_amb_sites = args.ignore_amb_sites

    if not quiet:
        print("\nRunning individual_hyde.py")

    data = hd.HydeData(
        infile, mapfile, outgroup, nind, ntaxa, nsites, quiet, ignore_amb_sites
    )

    # Read data into a HydeData object
    if not quiet:
        print("\nAnalyzing", len(triples), "triple(s).", sep=" ")

    outpath = hd.expand_prefix(prefix)
    if os.path.exists(outpath + "-ind.txt"):
        if not quiet:
            print(
                "\n**  Warning: File '"
                + outpath
                + "-ind.txt' already exists. **"
            )
        if not quiet:
            print("**  Renaming to '" + outpath + "-ind-old.txt'. **\n")
        os.rename(outpath + "-ind.txt", outpath + "-ind-old.txt")
        outfile = open(outpath + "-ind.txt", "a")
    else:
        outfile = open(outpath + "-ind.txt", "a")

    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=outfile,
    )

    for t in triples:
        res = data.test_individuals(t[0], t[1], t[2])
        write_ind(res, t, outfile)
