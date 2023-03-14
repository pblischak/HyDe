#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# hyde_gui.py
# Written by PD Blischak

import phyde as hd
from tkinter import *
from tkinter import ttk
import sys, os

try:
    from progress.bar import Bar
except ImportError:
    print("Couldn't import progress module...")
    print("Please use the following command to install it:")
    print("\n  pip install progress")
    sys.exit(0)

long_description = """\
Hybridization detection using phylogenetic invariants. \
For details on each analysis please consult the full documentation at \
https://hybridization-detection.rtfd.io/."""

quiet = (
    0  # Since we're using a GUI, don't need to worry about suppressing output
)


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


def write_out(out, triple, outfile):
    """
    Take the output from test_triple() and write it to file.
    """
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
    print(out["Zscore"], "\t", sep="", end="", file=outfile)
    print(out["Pvalue"], "\t", sep="", end="", file=outfile)
    print(out["Gamma"], "\t", sep="", end="", file=outfile)
    print(out["AAAA"], "\t", sep="", end="", file=outfile)
    print(out["AAAB"], "\t", sep="", end="", file=outfile)
    print(out["AABA"], "\t", sep="", end="", file=outfile)
    print(out["AABB"], "\t", sep="", end="", file=outfile)
    print(out["AABC"], "\t", sep="", end="", file=outfile)
    print(out["ABAA"], "\t", sep="", end="", file=outfile)
    print(out["ABAB"], "\t", sep="", end="", file=outfile)
    print(out["ABAC"], "\t", sep="", end="", file=outfile)
    print(out["ABBA"], "\t", sep="", end="", file=outfile)
    print(out["BAAA"], "\t", sep="", end="", file=outfile)
    print(out["ABBC"], "\t", sep="", end="", file=outfile)
    print(out["CABC"], "\t", sep="", end="", file=outfile)
    print(out["BACA"], "\t", sep="", end="", file=outfile)
    print(out["BCAA"], "\t", sep="", end="", file=outfile)
    print(out["ABCD"], "\n", sep="", end="", file=outfile)


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


def write_boot(boot, triple, outfile):
    """
    Write the current dictionary of bootstrapping output from bootstrap_triple()
    to the file passed as an argument to the function.
    """
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
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


def do_run_hyde(
    infile,
    mapfile,
    outgroup,
    nind,
    ntaxa,
    nsites,
    outfile,
    tripletsfile,
    pvalue,
    ignore_amb_sites,
):
    # Read in data as HydeData object
    data = hd.HydeData(
        infile, mapfile, outgroup, nind, ntaxa, nsites, quiet, ignore_amb_sites
    )

    # Get triples
    if tripletsfile != "none":
        triples = parse_triples(tripletsfile)
    else:
        triples = data.list_triples()

    print("\nAnalyzing", len(triples), "triple(s).", sep=" ")
    if tripletsfile != "none":
        print("Using triplets in file ", tripletsfile, ".", sep="")

    # Check if output file exists
    # and exit if it does
    if os.path.exists(outfile):
        print("Warning: Output file already exists...")
        print("Please choose different name.")
        return 0
    else:
        pass

    # print outfile header
    outf = open(outfile, "a")
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=outf,
    )

    filtered_outfile = open("filtered-" + outfile, "a")
    # print filtered outfile header
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=filtered_outfile,
    )

    print("")
    bar = Bar("Progress", max=len(triples))
    for t in triples:
        res = data.test_triple(t[0], t[1], t[2])
        write_out(res, t, outf)
        if (
            res["Pvalue"] < (pvalue / len(triples))
            and abs(res["Zscore"]) != 99999.9
            and res["Gamma"] > 0.0
            and res["Gamma"] < 1.0
        ):
            write_out(res, t, filtered_outfile)
        bar.next()

    print("\n\nDone!")


def do_ind_hyde(
    infile,
    mapfile,
    outgroup,
    nind,
    ntaxa,
    nsites,
    outfile,
    tripletsfile,
    pvalue,
    ignore_amb_sites,
):
    # Read in data as HydeData object
    data = hd.HydeData(
        infile, mapfile, outgroup, nind, ntaxa, nsites, quiet, ignore_amb_sites
    )

    # Get triples
    if tripletsfile != "none":
        triples = parse_triples(tripletsfile)
    else:
        print("Could not open triplets file...")
        return -1

    print("\nAnalyzing", len(triples), "triple(s).", sep=" ")
    print("Using triplets in file ", tripletsfile, ".", sep="")

    # Check if output file exists
    # and exit if it does.
    if os.path.exists(outfile):
        print("Warning: Output file already exists...")
        print("Please choose different name.")
        return 0
    else:
        pass

    # print outfile header
    outf = open(outfile, "a")
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=outf,
    )

    print("")
    bar = Bar("Progress", max=len(triples))
    for t in triples:
        res = data.test_individuals(t[0], t[1], t[2])
        write_ind(res, t, outf)
        bar.next()

    print("\n\nDone!")


def do_boot_hyde(
    infile,
    mapfile,
    outgroup,
    nind,
    ntaxa,
    nsites,
    outfile,
    tripletsfile,
    pvalue,
    ignore_amb_sites,
    reps,
):
    # Read in data as HydeData object
    data = hd.HydeData(
        infile, mapfile, outgroup, nind, ntaxa, nsites, quiet, ignore_amb_sites
    )

    # Get triples
    if tripletsfile != "none":
        triples = parse_triples(tripletsfile)
    else:
        triples = data.list_triples()

    print(
        "\nAnalyzing ",
        len(triples),
        " triple(s) (",
        reps,
        " bootstrap replicates each).",
        sep="",
    )
    if tripletsfile != "none":
        print("Using triplets in file ", tripletsfile, ".", sep="")

    # Check if output file exists
    # and exit if it does
    if os.path.exists(outfile):
        print("Warning: Output file already exists...")
        print("Please choose different name.")
        return 0
    else:
        pass

    # print outfile header
    outf = open(outfile, "a")
    print(
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\tABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n",
        end="",
        file=outf,
    )

    counter = 0
    print("")
    bar = Bar("Progress", max=len(triples))
    for t in triples:
        res = data.bootstrap_triple(t[0], t[1], t[2], reps)
        write_boot(res, t, outf)
        counter += 1
        if counter != len(triples):
            print("####\n", end="", file=outf)
        bar.next()

    print("\n\nDone!")


def check_input(
    infile,
    mapfile,
    outgroup,
    nind,
    ntaxa,
    nsites,
    outfile,
    tripletsfile,
    pvalue,
    ignore_amb_sites,
    reps,
    analysis,
):
    all_good = True
    if analysis == 0:
        print("Please select one of the available analyses...")
        all_good = False
    elif infile == "":
        print("No input file specified...")
        all_good = False
    elif mapfile == "":
        print("No map file specified...")
        all_good = False
    elif outgroup == "":
        print("No outgroup specified...")
        all_good = False
    elif nind == 0:
        print("Number of individuals is not specified...")
        all_good = False
    elif ntaxa == 0:
        print("Number of taxa is not specified...")
        all_good = False
    elif nsites == 0:
        print("Number of sites is not specified...")
        all_good = False
    else:
        pass

    return all_good


def main(
    infile,
    mapfile,
    outgroup,
    nind,
    ntaxa,
    nsites,
    outfile,
    tripletsfile,
    pvalue,
    ignore_amb_sites,
    reps,
    analysis,
):
    # Check the input
    good_to_go = check_input(
        infile,
        mapfile,
        outgroup,
        nind,
        ntaxa,
        nsites,
        outfile,
        tripletsfile,
        pvalue,
        ignore_amb_sites,
        reps,
        analysis,
    )
    if good_to_go:
        pass
    else:
        return 0

    if analysis == 1:
        do_run_hyde(
            infile,
            mapfile,
            outgroup,
            nind,
            ntaxa,
            nsites,
            outfile,
            tripletsfile,
            pvalue,
            ignore_amb_sites,
        )
    elif analysis == 2:
        do_ind_hyde(
            infile,
            mapfile,
            outgroup,
            nind,
            ntaxa,
            nsites,
            outfile,
            tripletsfile,
            pvalue,
            ignore_amb_sites,
        )
    elif analysis == 3:
        do_boot_hyde(
            infile,
            mapfile,
            outgroup,
            nind,
            ntaxa,
            nsites,
            outfile,
            tripletsfile,
            pvalue,
            ignore_amb_sites,
            reps,
        )
    else:
        print("Warning: Unrecognized option...")
        # sys.exit(-1)


if __name__ == "__main__":
    root = Tk()
    root.title("HyDe GUI v{}".format(hd.__version__))

    mainframe = ttk.Frame(root, padding="25 25 25 25")
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)

    infile = StringVar()
    mapfile = StringVar()
    outgroup = StringVar()
    nind = IntVar()
    ntaxa = IntVar()
    nsites = IntVar()
    outfile = StringVar()
    tripletsfile = StringVar(value="none")
    pvalue = DoubleVar(value=0.05)
    ignore_amb_sites = IntVar()
    reps = IntVar(value=100)
    analysis = IntVar()
    ttk.Label(
        mainframe, text=long_description, wraplength=600, justify=LEFT
    ).grid(column=1, row=1, columnspan=4, sticky=(W, E))
    infile_entry = ttk.Entry(mainframe, width=15, textvariable=infile)
    infile_entry.grid(column=2, row=2, sticky=(W, E))
    mapfile_entry = ttk.Entry(mainframe, width=15, textvariable=mapfile)
    mapfile_entry.grid(column=2, row=3, sticky=(W, E))
    outgroup_entry = ttk.Entry(mainframe, width=15, textvariable=outgroup)
    outgroup_entry.grid(column=2, row=4, sticky=(W, E))
    nind_entry = ttk.Entry(mainframe, width=15, textvariable=nind)
    nind_entry.grid(column=2, row=5, sticky=(W, E))
    ntaxa_entry = ttk.Entry(mainframe, width=15, textvariable=ntaxa)
    ntaxa_entry.grid(column=2, row=6, sticky=(W, E))
    nsites_entry = ttk.Entry(mainframe, width=15, textvariable=nsites)
    nsites_entry.grid(column=2, row=7, sticky=(W, E))
    outfile_entry = ttk.Entry(mainframe, width=15, textvariable=outfile)
    outfile_entry.grid(column=2, row=8, sticky=(W, E))
    triplets_entry = ttk.Entry(mainframe, width=15, textvariable=tripletsfile)
    triplets_entry.grid(column=2, row=9, sticky=(W, E))
    pvalue_entry = ttk.Entry(mainframe, width=15, textvariable=pvalue)
    pvalue_entry.grid(column=2, row=10, sticky=(W, E))
    ignore_box = ttk.Checkbutton(
        mainframe,
        text="Ignore missing/ambiguous sites?",
        variable=ignore_amb_sites,
    )
    ignore_box.grid(column=2, row=11, sticky=(W, E))
    run_hyde = ttk.Radiobutton(
        mainframe, text="Full HyDe analysis", variable=analysis, value=1
    )
    run_hyde.grid(column=2, row=13, sticky=(W, E))
    ind_hyde = ttk.Radiobutton(
        mainframe,
        text="Individual-based HyDe analysis",
        variable=analysis,
        value=2,
    )
    ind_hyde.grid(column=2, row=14, sticky=(W, E))
    boot_hyde = ttk.Radiobutton(
        mainframe,
        text="Individual-based bootstrapping analysis",
        variable=analysis,
        value=3,
    )
    boot_hyde.grid(column=2, row=15, sticky=(W, E))
    reps_entry = ttk.Entry(mainframe, width=4, textvariable=reps)
    reps_entry.grid(column=3, row=15, sticky=W)
    ttk.Separator(mainframe, orient=HORIZONTAL).grid(
        column=1, row=16, columnspan=4, sticky=(W, E)
    )

    ttk.Button(
        mainframe,
        text="Run",
        command=lambda: main(
            infile.get(),
            mapfile.get(),
            outgroup.get(),
            nind.get(),
            ntaxa.get(),
            nsites.get(),
            outfile.get(),
            tripletsfile.get(),
            pvalue.get(),
            ignore_amb_sites.get(),
            reps.get(),
            analysis.get(),
        ),
    ).grid(column=4, row=17, sticky=W)

    ttk.Button(mainframe, text="Quit", command=lambda: sys.exit(0)).grid(
        column=1, row=17, sticky=W
    )

    ttk.Label(mainframe, text="Input file:").grid(column=1, row=2, sticky=W)
    ttk.Label(mainframe, text="Map file:").grid(column=1, row=3, sticky=W)
    ttk.Label(mainframe, text="Outgroup:").grid(column=1, row=4, sticky=W)
    ttk.Label(mainframe, text="Number of individuals:").grid(
        column=1, row=5, sticky=W
    )
    ttk.Label(mainframe, text="Number of taxa:").grid(column=1, row=6, sticky=W)
    ttk.Label(mainframe, text="Number of sites:").grid(
        column=1, row=7, sticky=W
    )
    ttk.Label(mainframe, text="Output file:").grid(column=1, row=8, sticky=W)
    ttk.Label(mainframe, text="Triplets file:").grid(column=1, row=9, sticky=W)
    ttk.Label(mainframe, text="P-value:").grid(column=1, row=10, sticky=W)
    ttk.Label(mainframe, text="Analysis options:").grid(
        column=1, row=12, sticky=W
    )
    ttk.Label(mainframe, text="reps").grid(column=4, row=15, sticky=W)

    for child in mainframe.winfo_children():
        child.grid_configure(padx=5, pady=3)

    root.bind(
        "<Return>",
        lambda event: main(
            infile.get(),
            mapfile.get(),
            outgroup.get(),
            nind.get(),
            ntaxa.get(),
            nsites.get(),
            outfile.get(),
            tripletsfile.get(),
            pvalue.get(),
            ignore_amb_sites.get(),
            reps.get(),
            analysis.get(),
        ),
    )

    windowWidth = root.winfo_reqwidth()
    windowHeight = root.winfo_reqheight()
    positionRight = int(root.winfo_screenwidth() / 3 - windowWidth / 2)
    positionDown = int(root.winfo_screenheight() / 4 - windowHeight / 2)
    root.geometry("+{}+{}".format(positionRight, positionDown))

    root.mainloop()
