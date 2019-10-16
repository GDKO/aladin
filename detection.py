#!/usr/bin/env python3

"""
detection.py
    Usage:
    ./detection.py -g <FILE> -f <FILE> -o <FILE> -i <STR>

    Options:
        -h, --help                     show options
        -g, --gfa_file <FILE>          gfa file
        -f, --fasta_file <FILE>        fasta file
        -o, --output_file <FILE>       output file
        -i, --contig_id <STR>          contig id
        --version                      print version
"""

import re
from docopt import docopt


def forma_cigar(cigar):
    cigar = cigar.rstrip()
    type_counts = re.split("\D+", cigar)
    types = re.split("\d+", cigar)
    types.pop(0)
    length_cigar = 0
    for c, t in zip(type_counts, types):
        if t != "I":
            length_cigar += int(c)
    return length_cigar

def forma_cigar_minus(cigar):
    cigar = cigar.rstrip()
    type_counts = re.split("\D+", cigar)
    types = re.split("\d+", cigar)
    types.pop(0)
    length_cigar = 0
    for c, t in zip(type_counts, types):
        if t != "D":
            length_cigar += int(c)
    return length_cigar

def return_pos(file_gfa, contig):

    file_opened = open(file_gfa, "r")
    cigar_plus = ""
    length_cigar = 0

    for line in file_opened:
        if line.startswith("L"):
            elem = line.split("\t")
            if elem[1] == contig and elem[1] == elem[3] and elem[2] == "+" and elem[4] == "+":
                cigar_plus = elem[5]
                length_cigar = forma_cigar(cigar_plus)
    file_opened.close()

    file_opened = open(file_gfa, "r")
    if length_cigar == 0:
        for line in file_opened:
            if line.startswith("L"):
                elem = line.split("\t")
                if elem[1] == contig and elem[1] == elem[3] and elem[2] == "-" and elem[4] == "-":
                    cigar_minus = elem[5]
                    length_cigar = forma_cigar_minus(cigar_minus)

    file_opened.close()
    length_plus = length_cigar
    return length_plus

############### Main ###############

def main():
    args = docopt(__doc__, version='1.00')
    file_gfa = str(args['--gfa_file'])
    file_fasta = str(args['--fasta_file'])
    file_output = str(args['--output_file'])
    contig = str(args['--contig_id'])

    length_plus = return_pos(file_gfa, contig)
    sequence = ""

    file_seq = open(file_fasta, "r")

    for line in file_seq:
        line = line.rstrip("\n")
        if line.startswith(">"):
            header = line
        else:
            sequence += line

    file_seq.close()

    file_new_seq = open(file_output, 'w')
    file_new_seq.write(header + "\n")
    file_new_seq.write(sequence[length_plus:] + "\n")
    file_new_seq.close()


if __name__ == "__main__":
    main()
