#!/usr/bin/env python3

"""
    Usage:
    ./reorder.py -r <FILE> -f <FILE>

    Options:
    -h, --help                     show options
    -r, --reference <FILE>         reference sequence input
    -f, --fasta <FILE>             fasta to reorder
    --version                      print version
"""
import os
import subprocess
from Bio import SeqIO, SearchIO, Seq
from docopt import docopt
from depot.AladinFun import fichier_vide

def reorder(reference,fasta):
    FNULL = open(os.devnull, 'w')
    blat_file = fasta + ".psl"
    subprocess.call(["blat", fasta, reference, blat_file],stdout=FNULL,stderr=FNULL)
    fichier_vide(blat_file)

    qresult = SearchIO.read(blat_file, 'blat-psl')
    hit_start = qresult[0][0].hit_range[0]
    hit_end = qresult[0][0].hit_range[1]
    hit_id = qresult[0].id
    orientation = qresult[0][0].query_strand

    if orientation > 0:
        orientation = "+"
    else:
        orientation = "-"

    assembly = {}
    for record in SeqIO.parse(fasta, "fasta"):
        elem = str(record.id).split(" ")
        if orientation=="+":
            assembly[elem[0]] = str(record.seq)
        else:
            assembly[elem[0]] = str(record.seq.reverse_complement())

    circ_file = open(fasta+".reordered", 'w')
    circ_file.write(">" + hit_id + "\n")
    if orientation=="+":
        pos = hit_start
        circ_file.write(assembly[hit_id][pos:] + assembly[hit_id][:pos] + "\n")
    else:
        pos = len(assembly[hit_id])-hit_end
        circ_file.write(assembly[hit_id][pos:] + assembly[hit_id][:pos] + "\n")
    circ_file.close()

    return(pos,orientation)

def main():
    args = docopt(__doc__, version='1.1')
    reference = os.path.abspath(str(args['--reference']))
    fasta = os.path.abspath(str(args['--fasta']))

    reorder_position, orientation = reorder(reference,fasta)

    print("Reordered fasta file at position " + str(reorder_position) + " with " + str(orientation) + " orientation.")

if __name__ == "__main__":
    main()
