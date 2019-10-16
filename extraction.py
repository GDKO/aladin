#!/usr/bin/env python3

"""
extraction.py
    Usage:
    ./extraction.py -r <FILE> -i <FILE>

    Options:
        -h, --help                     show options
        -r, --reads <FILE>             reads file
        -i, --ids <FILE>               read ids
        --version                      print version
"""

from Bio import SeqIO
from docopt import docopt
from depot.AladinFun import open_files, return_format


def extract_reads(reads, read_ids, outfile):

    file_in = open_files(reads)
    file_out = open(outfile, 'w')

    for record in SeqIO.parse(file_in, return_format(reads)):
        if record.id in read_ids:
            file_out.write(">" + str(record.id) + "\n")
            file_out.write(str(record.seq) + "\n")

    file_in.close()
    file_out.close()

############### Main ###############
def main():
    args = docopt(__doc__, version='1.00')
    reads = str(args['--reads'])
    ids = str(args['--ids'])

    file_list = open_files(ids)

    read_ids = []
    for line in file_list:
        line = line.rstrip("\n")
        read_ids.append(line)

    file_list.close()
    outfile = "read_ids.txt.fasta"
    extract_reads(reads, read_ids, outfile)

if __name__ == "__main__":
    main()
