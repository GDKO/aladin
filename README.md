ALADIN
===============
### Overview

ALADIN (mitochondriAL circulAr Dna reconstItutioN) uses DNA long reads (Nanopore or Pacbio) and a reference sequence to assemble the mitochondrial genome.

Requirements
------------
* Programs
  - [minimap2](https://github.com/lh3/minimap2) - v2.17 or later
  - [miniasm](https://github.com/lh3/miniasm) - r179 or later
  - [Minipolish](https://github.com/rrwick/Minipolish) - v0.1.3 or later
  - [aragorn](http://mbio-serv2.mbioekol.lu.se/ARAGORN) - v1.2.38
  - [seqtk](https://github.com/lh3/seqtk) - v1.3-r106
  - [racon](https://github.com/isovic/racon) - v.1.4.3 or later
* Python Libraries
  - Biopython
  - docopt

Quick usage
------------
```aladin.py -r reference -i reads```

Full usage
-----------
```
    Usage:
    ./aladin.py -r <FILE> -i <FILE> [-o <DIR>] [-d <STR>] [-m <STR>] [-l <INT>] 
                                    [-f <FLOAT>] [-s <STR>] [-t <INT>] [--cleanup]

    Options:
        -h, --help                     show options
        -r, --reference <FILE>         reference input
        -i, --reads <FILE>             Reads input
        -o, --dir_output <DIR>         creates a directory for all output files [default: results]
        -d, --data_format <STR>        (N) Nanopore or (P) Pacbio [default: N]
        -m, --mode <STR>               (M) Mitochondion or (C) Chloroplast [default: M]
        -l, --length <INT>             break reads into chunks of this length [default: 4000]
        -f, --fraction <FLOAT>         fraction of length at which the end of the sequence gets split into a new sequence [default: 0.1]
        -s, --genome_size <INT>        set expected size in kb [default: 20]
        -t, --threads <INT>            number of threads [default: 3]
        --cleanup                      remove intermediate files
        --version                      print version
```
