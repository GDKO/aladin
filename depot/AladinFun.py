import sys
import os
import shutil
import gzip
from Bio import SeqIO

# Reporting
def error(message, *argv):
    if argv is None:
        sys.exit(error_d[message])
    else:
        sys.exit(error_d[message] % (argv))

error_d = {
    '0': '[ERROR:0]\t: File %s does not exist.',
    '1': '[ERROR:1]\t: File %s is empty.',
    '2': '[ERROR:2]\t: %s not found! Please install and add to the PATH variable.'
}

status_d = {
    '0': '[+] Mapping reads with minimap2.',
    '1': '[+] Extracting reads.',
    '2': '[+] Assembling the reads.',
    '3': '[+] Finding the mitochondrion.',
    '4': '[+] Pipeline finished.'
}


# Check if empty or does not exist
def fichier_vide(fichier):
    if os.path.isfile(fichier) is True:
        if os.path.getsize(fichier) <= 200:
            error('1', fichier)
    else:
        error('0', fichier)

# Check if program is in path

def check_programs(*arg):
    error_list = []
    for program in arg:
        if which(program) is False:
            error_list.append(program)
    if error_list:
        programs = ', '.join(error_list)
        error('2', programs)

def which(program):
    if shutil.which(program):
        return program
    else:
        return False

# Directory checking

def get_outdir(out_directory, add_dir=""):
    """generates output directory in case it does not exist."""
    if type(out_directory) != str:
        print("\t[!] {} is NOT a directory! Please specify an output directory".format(out_directory))
        sys.exit()
    elif os.path.isfile(out_directory):
        print("\t[!] {} is a File! Please specify an output directory".format(out_directory))
        sys.exit()
    elif not os.path.exists(os.path.join(out_directory, add_dir)):
        os.mkdir(os.path.join(out_directory, add_dir))
        return os.path.abspath(os.path.join(out_directory, add_dir))
    else:
        return os.path.abspath(os.path.join(out_directory, add_dir))

def check_indir(input_dir):
    if not os.path.exists(input_dir):
        print("\t[!] FATAL ERROR: '{}' directory not found".format(input_dir))
        sys.exit()
    elif not os.path.isdir(input_dir):
        print("\t[!] FATAL ERROR: '{}' is not a directory".format(input_dir))
        sys.exit()
    else:
        return input_dir

# Can open gzip files

def open_files(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname, 'rt', encoding="latin-1")
    else:
        return open(fname, 'rt', encoding="latin-1")

# Check if fasta

def is_fasta(filename):
    with open_files(filename) as handle:
        if handle.readline().startswith(">"):
            return True
        else:
            return False

def return_format(filename):
    if is_fasta(filename):
        return "fasta"
    else:
        return "fastq"
