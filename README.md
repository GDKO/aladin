ALADIN
===============
### Overview

ALADIN (mitochondriAL circulAr Dna reconstItutioN) use DNA long reads (Nanopore or Pacbio) and a mitochondrial reference of the same species or the close on, and get the mitochondrial genome of the reads species

Requirements
------------
* Programs
  - [minimap2](https://github.com/lh3/minimap2) - v2.17 or later
  - [canu](https://github.com/marbl/canu) - v1.8 or later
  - [mummer](https://github.com/mummer4/mummer) - v4 or later
  - [aragorn](http://mbio-serv2.mbioekol.lu.se/ARAGORN) - v1.2.38
* Python Libraries
  - Biopython
  - docopt

Quick usage
------------
```aladin.py -r reference -i reads -d N```

Full usage
-----------
```
Usage:
./aladin.py -r <FILE> -i <FILE> [-o <DIR>] [-d <STR>] [-m <STR>] [-l <INT>] [-f <FLOAT>] [-s <STR>] [-t <INT>] [--cleanup]

Options:
    -h, --help                     show options
    -r, --reference <FILE>         reference input
    -i, --reads <FILE>             Reads input
    -o, --dir_output <DIR>         creates a directory for all output files [default: results]
    -d, --data_format <STR>        (N) Nanopore or (P) Pacbio [default: N]
    -m, --mode <STR>               (M) Mitochondion or (C) Chloroplast [default: M]
    -l, --length <INT>             break reads into chunks of this length [default: 4000]
    -f, --fraction <FLOAT>         fraction of length at which the end of the sequence gets split into a new sequence [default: 0.1]
    -s, --genome_size <INT>        set expected mitichondrion size in kb [default: 20]
    -t, --threads <INT>            minimap2 threads [default: 3]
    --cleanup                      remove intermediate files
    --version                      print version
```
