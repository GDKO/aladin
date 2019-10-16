#!/usr/bin/env python3

"""
fragmenter.py
    Usage:
    ./fragmenter.py -i <FILE> [-l <INT>] [-f <FLOAT>]

    Options:
        -h, --help                     show options
        -i, --input_file <FILE>        file to fragment
        -l, --length <INT>             split reads into chunks of this length [default: 4000]
        -f, --fraction <FLOAT>         fraction of length at which the end of the sequence gets split into a new sequence [default: 0.1]
        --version                      print version
"""

from Bio import SeqIO
from docopt import docopt
from depot.AladinFun import open_files, return_format

def fragment_fasta(input_file, length_frag, fraction_select):
    file_in = open_files(input_file)

    outfile = input_file+".frag.fasta"
    file_out = open(outfile, 'w')

    for record in SeqIO.parse(file_in, return_format(input_file)):
        name = record.id
        seq = str(record.seq)
        if length_frag < 1:
            file_out.write(">" + name + "\n")
            file_out.write(seq + "\n")
        else:
            subseqs = [seq[i:i+length_frag] for i in range(0, len(seq), length_frag)]
            if len(subseqs) > 1 and len(subseqs[-1]) < (fraction_select*length_frag):
                toadd = subseqs.pop()
                subseqs[-1] += toadd
            i = 0
            for subseq in subseqs:
                file_out.write(">" + name + "_" + str(i) + "\n")
                file_out.write(subseq + "\n")
                i += 1

    file_in.close()
    file_out.close()

############### Main ###############
def main():
    args = docopt(__doc__, version='1.00')
    input_file = str(args['--input_file'])
    length_frag = int(args['--length'])
    fraction_select = float(args['--fraction'])

    fragment_fasta(input_file, length_frag, fraction_select)

if __name__ == "__main__":
    main()
