#!/usr/bin/env python3

"""
aladin.py
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
        -s, --genome_size <INT>        set expected size in kb [default: 20]
        -t, --threads <INT>            minimap2 threads [default: 3]
        --cleanup                      remove intermediate files
        --version                      print version
"""

import os
import re
import shutil
import subprocess
from Bio import SeqIO
from docopt import docopt
from depot.AladinFun import *
from detection import return_pos
from extraction import extract_reads
from fragmenter import fragment_fasta

def run_nucmer(assembly_file, contig):
    FERR = open("ALADIN.err", 'a')
    coords_file = "nucmer/nucmer.coords"
    delta_file = "nucmer/nucmer.delta"
    coords_out = open(coords_file, 'w')
    subprocess.call(["nucmer", "--maxmatch", "-c", "100", "--delta",
                     delta_file, assembly_file, assembly_file],
                     stdout=FERR, stderr=FERR)
    subprocess.call(["show-coords", "-B", "-r",
                     delta_file],
                     stdout=coords_out, stderr=FERR)
    coords_file = open_files(coords_file)
    length_coords = []
    contigs = []
    trim_list = []
    for line in coords_file:
        elems = line.split("\t")
        q_contig_name = elems[0]
        contig_length = int(elems[2])
        s_contig_name = elems[5]
        q_start = int(elems[6])
        q_end = int(elems[7])
        s_start = int(elems[8])
        s_end = int(elems[9])
        strand = elems[17]
        #check that the other hit is close to the end
        if q_contig_name  == s_contig_name and strand == "Plus" and q_contig_name == contig:
            if q_start < 5 and s_start > 1000 and (s_end + 100) > contig_length:
                if q_end not in trim_list:
                    length_coords.append(q_end)
                    contigs.append(q_contig_name)
                    trim_list.append(q_end)
            elif s_start < 5 and q_start > 1000 and (q_end + 100) > contig_length:
                if s_end not in trim_list:
                    length_coords.append(s_end)
                    contigs.append(q_contig_name)
                    trim_list.append(s_end)
    return (contigs, length_coords, len(length_coords))

def run_pipeline(data_format, mode, reference, reads, genome_size, length, fraction, threads):
    # Directory structure
    get_outdir("map")
    get_outdir("nucmer")

    # File naming
    FERR = open("ALADIN.err", 'w')
    sam_file = "map/aln.sam"
    extract_reads_file = "map/read_ids.txt.fasta"
    extract_reads_frag_file = "map/read_ids.txt.fasta.frag.fasta"
    extract_reads_frag_file_100x = "map/read_ids.txt.fasta.frag.100x.fasta"
    assembly_file = "canu/extract.contigs.fasta"
    aragorn_file = "canu/aragorn.txt"
    gfa_file = "canu/extract.contigs.gfa"

    #Set data format
    if data_format == "N":
        format_mini = "map-ont"
        format_canu = "-nanopore-raw"
    elif data_format == "P":
        format_mini = "map-pb"
        format_canu = "-pacbio-raw"

    # Map reads
    print(status_d['0'])
    sam_out = open(sam_file, 'w')
    subprocess.call(["minimap2", "--sam-hit-only", "-t", threads, "--secondary=no", "-ax",
                     format_mini, reference, reads],
                     stdout=sam_out, stderr=FERR)
    sam_out.close()
    fichier_vide(sam_file)

    # Extract and fragment reads
    print(status_d['1'])
    read_ids = []
    sam_input = open_files(sam_file)
    for line in sam_input:
        elems = line.split("\t")
        if not line.startswith("@") and str(elems[4]) == "60":
            read_ids.append(elems[0])

    read_ids = list(set(read_ids))
    extract_reads(reads, read_ids, extract_reads_file)
    fragment_fasta(extract_reads_file, length, fraction)

    # Select 100x of reads
    num_frags = len([1 for line in open(extract_reads_frag_file) if line.startswith(">")])
    num_100x = (genome_size*1000*200)/length
    fraction_select = min(num_100x/num_frags,1)
    seqtk_out = open(extract_reads_frag_file_100x,'w')
    subprocess.call(["seqtk", "seq", "-f", str(fraction_select), extract_reads_frag_file],
                     stdout=seqtk_out,stderr=FERR)
    seqtk_out.close()

    # Run canu
    print(status_d['2'])
    genome_size = "genomeSize=" + str(genome_size) + "k"
    subprocess.call(["canu", "-d", "canu", "-p", "extract", genome_size, "corOutCoverage=200", "correctedErrorRate=0.15",
                     format_canu, extract_reads_frag_file_100x],
                     stdout=FERR, stderr=FERR)
    fichier_vide(assembly_file)

    # Run aragorn to select contig
    print(status_d['3'])
    aragorn_out = open(aragorn_file, 'w')
    if mode == "M":
        subprocess.call(["aragorn", "-l", "-gc5", "-mt", "-w", assembly_file],
                        stdout=aragorn_out, stderr=FERR)
    else:
        subprocess.call(["aragorn", "-l", "-gc1", "-m", "-t", "-w", assembly_file],
                        stdout=aragorn_out, stderr=FERR)
    aragorn_input = open_files(aragorn_file)
    trna_contigs = []
    trna_genes = []
    for line in aragorn_input:
        if line.startswith(">") and not line.startswith(">end"):
            elem = line.split(" ")
            header = elem[0].replace(">","")
            trna_contigs.append(header)
        elif "genes found" in line:
            elem = line.split(" ")
            num_trnas = int(elem[0])
            trna_genes.append(num_trnas)
    aragorn_contig = trna_contigs[trna_genes.index(max(trna_genes))]

    print("[!] Found contig " + str(aragorn_contig) + " with " +  str(max(trna_genes)) + " mtRNAs")

    num_circ = 0
    length_canu = 0
    trim_len = 0
    assembly = {}

    # Check canu output
    for record in SeqIO.parse(assembly_file, "fasta"):
        elem = str(record.id).split(" ")
        if re.search(aragorn_contig,str(record.id)):
            if re.search("suggestCircular=yes", str(record.description)):
                assembly[elem[0]] = str(record.seq)
                num_circ = 1

    contigs, length_coords, nucmer_cases = run_nucmer(assembly_file, aragorn_contig)
    mito_contig = aragorn_contig

    if num_circ == 1:
        length_canu = return_pos(gfa_file, aragorn_contig)
        trim_len = length_canu
        if str(nucmer_cases) == "0":
            print(warn_d['2'] % trim_len)
        elif str(nucmer_cases) == "1":
            abs_trim = abs(int(length_canu) - int(length_coords[0]))
            if abs_trim < 11:
                print(status_d['4'] % (trim_len, length_coords[0]))
            else:
                print(warn_d['1'] % abs_trim)
        else:
            print(warn_d['3'])
    else:
        print(warn_d['0'])
        if str(nucmer_cases) == "0":
            warn_d('4')
        elif str(nucmer_cases) == "1":
            trim_len = length_coords[0]
            print(status_d['5'] % trim_len)
        else:
            warn_d('5')

    circ_file = open("Canu_trimmed.fasta", 'w')
    circ_file.write(">" + mito_contig + "\n")
    circ_file.write(assembly[mito_contig][int(trim_len):] + "\n")
    circ_file.close()

    FERR.close()


############### Main ###############

def main():
    args = docopt(__doc__, version='1.01')
    reference = os.path.abspath(str(args['--reference']))
    reads = os.path.abspath(str(args['--reads']))
    dir_output = get_outdir(args['--dir_output'])
    length = int(args['--length'])
    fraction = float(args['--fraction'])
    threads = str(args['--threads'])
    cleanup = args['--cleanup']
    data_format = str(args['--data_format'])
    mode = str(args['--mode'])
    genome_size = int(args['--genome_size'])

    os.chdir(dir_output)

    check_programs("minimap2", "canu", "nucmer", "show-coords", "aragorn", "seqtk")

    run_pipeline(data_format, mode, reference, reads, genome_size, length, fraction, threads)

    if cleanup:
        shutil.rmtree("canu")
        shutil.rmtree("map")
        shutil.rmtree("nucmer")

    print(status_d['6'])
    resultats_dir = os.getcwd()
    print("Results be in "+resultats_dir+"/Canu_trimmed.fasta")

if __name__ == "__main__":
    main()
