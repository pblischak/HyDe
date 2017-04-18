## HyDe: hybridization detection using phylogenetic invariants

## Installation

```bash
git clone https://github.com/pblischak/HyDe.git
cd HyDe/src

# If you don't have OpenMP for multithreading
make

# If you do have OpenMP for multithreading
make OPENMP=yes

# Copy hyde and hyde.py to /usr/local/bin
sudo make install
```

## Running HyDe

Type `hyde.py -h` for options.

```
hyde.py -i <infile> -m <map-file> -o <outgroup> \
        -n <num-ind> -t <num-taxa> -s <num-sites> \
        -j <num-threads> --prefix <prefix>
```
