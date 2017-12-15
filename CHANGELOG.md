# Change log

## v0.4.0

**Removing dependence on `hyde_cpp`**

 - The `hyde_cpp` executable has been deprecated in favor of using the
   `*_hyde.py` scripts, which are just as fast. It will be available
   as its own GitHub repository for legacy purposes.

 - The `phyde.analyze` submodule has been removed, and all functions calling the
   `hyde_cpp` program have been either modified or removed.

 - We anticipate that this will make the installation process simpler, and will
   allow us to distribute the software as purely a Python package (w/ Cython for speed).

## v0.3.3

**New multithreaded scripts**

 - All of the `*_hyde.py` scripts now have multithreaded versions:
   `run_hyde_mp.py`, `individual_hyde_mp.py`, and `bootstrap_hyde_mp.py`.
   Parallelization with these scripts works best with data sets that have a lot of species.
   We use the `multiprocess` module, which will need to be installed: `pip install multiprocess`.

 - We also added a quiet option (`-q` or `--quiet`) to suppress printing to stdout
   while the program is running.

## v0.3.2

**Input file format update**

 - All of the `*_hyde.py` scripts now accept DNA sequence data files
   in sequential Phylip format (ie, you don't need to remove the header line).
   `hyde_cpp` also accepts data in this format.
 - We maintain backwards compatibility with the previous input format
   without the header information as well.

## v0.3.1

**Bug fixes**

 - Corrected bootstrap_hyde.py output.
 - Bootstrap class now correctly calculates # of bootstrap reps.
 - Better handling of files with extra empty line at the end.

## v0.3.0

**Changes introducing incompatibilities with previous versions:**

 - The ability to do bootstrapping was removed from `hyde_cpp`. This
   functionality is now provided by the `bootstrap_hyde.py` script. All
   scripts relying on bootstrapping with `hyde_cpp` have been changed
   accordingly.

**New workflow scripts:**

We introduce a new workflow for HyDe using three main Python scripts that can be
called from the command line. They are automatically installed with the `phyde`
module. A new input file is also introduced for specifying specific triples that
are to be tested using the methods implemented in the scripts (see below).

 - `run_hyde.py`: this script remains largely the same as in previous versions.
   We have added the ability to test specific hypotheses rather than running a
   full analysis (all triples in all directions) as well.
 - `individual_hyde.py`: this new script takes a table of specified triples
   and tests all individuals in the hybrid population separately.
 - `bootstrap_hyde.py`: this new script will bootstrap resample individuals
   within the hybrid populations that are specified by the table of triples.

**Input format for specifying triples:**

The three scripts listed above can take as input a three column table where each
row specifies the names of taxa that are to be tested for hybridization.
The order is (Parent 1, Hybrid, Parent 2) as below:

```
sp1 sp2 sp3
sp1 sp2 sp4
sp2 sp3 sp4
.
.
.
```

The `individual_hyde.py` and `bootstrap_hyde.py` files can also take output files
from previous runs of HyDe. These contain header information that the scripts will
automatically check for and will then parse the first three columns of the file.

## v0.2.3

### **Final release of v0.2.* series.**

This release can be used as a legacy version for the old workflow for HyDe.
Ideally users should switch to the newest version (see above).

## v0.2.0

First release of HyDe.

**Features:**

`hyde_cpp` executable for performing full scale hybridization detection analyses on phylogenomic data sets.

`phyde` Python module for performing individual hypothesis tests, parsing and analyzing results from hyde_cpp, plotting, and calculating D-Statistics.
