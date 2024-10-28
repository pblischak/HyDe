#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# bootstrap_hyde.py
# Written by PD Blischak

"""
<< bootstrap_hyde_mp.py >>

Bootstrap resample individuals within a putative hybrid population for a
specified triple -- multithreaded version.

Arguments
---------

For more details on script arguments, type: bootstrap_hyde_mp.py -h

    - infile         <string> : name of the DNA sequence data file.
    - mapfile        <string> : name of the taxon map file.
    - outgroup       <string> : name of the outgroup.
    - triples        <string> : name of the file containing triples for testing.
    - reps              <int> : number of bootstrap replicates [default=100].
    - nind              <int> : number of sampled individuals.
    - nsites            <int> : number of sampled sites.
    - ntaxa             <int> : number of sampled taxa/populations.
    - threads           <int> : number of threads [default=all available]
    - prefix         <string> : name added to the beginning of output file.
    - quiet            <flag> : suppress printing to stdout.
    - ignore_amb_sites <flag> : ignore missing/ambiguous bases.

Output
------

    Writes a file ('hyde-boot.txt') with the bootstrap results for each
    triple. These are Written in the same format as the normal hyde results
    file except that the bootstrap replicates for each triple are separated
    by four pound symbols and a new line "####\\n" (can be used to split the
    bootstrap replicates with the split() function)."""

import phyde as hd
import argparse
import sys
import os
import multiprocess as mp


def parse_triples(triples_file):
    """
    Parse a three column table or previous results file to get the names
    of the taxa that are going to be tested for hybridization.

    Returns a list of three-tuples of the form (P1, Hybrid, P2) for all triples
    that are to be tested.
    """
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
        return triples
    triples = [
        (l.split()[0], l.split()[1], l.split()[2]) for l in lines
    ]
    return triples


def write_boot(boot, triple, outfile):
    """
    Write the current dictionary of bootstrapping output from bootstrap_triple()
    to the file passed as an argument to the function.
    """
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC"
        "\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=outfile,
    )
    for _, v in boot.items():
        print(
            triple[0],
            "\t",
            triple[1],
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
    Runs the script.
    """
    # print docstring if only the name of the script is given
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Options for bootstrap_hyde.py", add_help=True
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
        "-r",
        "--reps",
        action="store",
        type=int,
        default=100,
        metavar="\b",
        help="number of bootstrap replicates [default=100]",
    )
    additional.add_argument(
        "-j",
        "--threads",
        action="store",
        type=int,
        default=mp.cpu_count(),
        metavar="\b",
        help="number of threads [default=all available]",
    )
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
        help="ignore missing/ambiguous bases.",
    )

    args = parser.parse_args()
    infile = args.infile
    mapfile = args.map
    outgroup = args.outgroup
    triples = parse_triples(args.triples)
    nind = args.num_ind
    ntaxa = args.num_taxa
    nsites = args.num_sites
    reps = args.reps
    prefix = args.prefix
    quiet = args.quiet
    ignore_amb_sites = args.ignore_amb_sites
    threads = args.threads

    if not quiet:
        print("\nRunning bootstrap_hyde_mp.py")

    # Read data into a HydeData object
    data = hd.HydeData(
        infile, mapfile, outgroup, nind, ntaxa, nsites, quiet, ignore_amb_sites
    )

    if not quiet:
        print(
            "\nAnalyzing ",
            len(triples),
            " triple(s) (",
            reps,
            " bootstrap replicates each) using ",
            threads,
            " thread(s).",
            sep="",
        )

    outpath = hd.expand_prefix(prefix)
    if os.path.exists(outpath + "-boot.txt"):
        if not quiet:
            print(
                "\n**  Warning: File '"
                + outpath
                + "-boot.txt' already exists. **"
            )
        if not quiet:
            print("**  Renaming to '" + outpath + "-boot-old.txt'. **\n")
        os.rename(prefix + "-boot.txt", outpath + "-boot-old.txt")
        outfile = open(outpath + "-boot.txt", "a")
    else:
        outfile = open(outpath + "-boot.txt", "a")

    def wrap_test(tr):
        """
        Wrapper function for running a hypothesis
        test on a given triple.
        """
        res = {(tr[0], tr[1], tr[2]): data.bootstrap_triple(
            tr[0], tr[1], tr[2], reps
        )}
        return res

    def mp_run():
        """
        Run tests on multiple threads.
        """
        p = mp.Pool(threads)
        res = p.map(wrap_test, triples)
        return res

    # run the bootsrap replicates
    out = mp_run()
    counter = 0
    for o in out:
        key = list(o.keys())[0]
        value = list(o.values())[0]
        write_boot(value, key, outfile)
        counter += 1
        if counter != len(triples):
            print("####\n", end="", file=outfile)
