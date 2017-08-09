# Change log

## v0.3.0 (forthcoming)

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

## v0.2.0

First release of HyDe.

**Features:**

`hyde_cpp` executable for performing full scale hybridization detection analyses on phylogenomic data sets.

`phyde` Python module for performing individual hypothesis tests, parsing and analyzing results from hyde_cpp, plotting, and calculating D-Statistics.
