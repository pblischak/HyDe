## HyDe: hybridization detection using phylogenetic invariants

## Installation

```bash
git clone https://github.com/pblischak/HyDe HyDe
cd HyDe/src

# If you don't have OpenMP for multithreading
make

# If you do have OpenMP for multithreading
make OPENMP=yes

# Copy to hyde and hyde.py to /usr/local/bin
sudo make install
```

## Running HyDe

Two main ways to run HyDe.

```
hyde.py -n <num-ind> -t <num-taxa> -s <num-sites> -i <infile> -m <map-file> -j <num-threads>
```
