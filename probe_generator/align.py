#!/usr/bin/env python3
"""This script performs multiple sequence alignment of an input FASTA file.
The user can specify the alignment software from the following options:
ClustalO (default), MUSCLE, MAFFT, and T-Coffee.
The user must have the alignment software installed and available.
"""

import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description='Perform sequence alignment using ClustalO, MUSCLE, MAFFT or T-Coffee.')

parser.add_argument('-i', '--input', required=True, help='Path to the file containing the sequences to be aligned.') # Defines an argument, the input file is required.
parser.add_argument('-o', '--output', default='out.fasta', help='Name of the output file.') # output file, default is out.fasta
parser.add_argument('-p', '--program', default='clustalo', help='Program to be used for alignment: clustalo, muscle, mafft, tcoffee.')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to be used.')

args = parser.parse_args() # Retrieve the arguments.

input_file = args.input # Use the values of the arguments.
aligned_file = args.output # `aligned_file` contains the name of the output file.
program = args.program.lower() # The program name is converted to lowercase for convenience.
threads = int(args.threads) # contains the numbers of threads to be used.

allowed_programs = ['clustalo', 'muscle', 'mafft', 'tcoffee'] # List of supported programs.

if threads < 1: # The numbers of threads is 0 or lower (imposible)
    print("The number of threads must be greater than or equal to 1.")
elif program in allowed_programs:
    try:
        if(program == 'muscle'):
            if threads:
                print("WARNING: MUSCLE does not support multiple threads.  The program will be executed with a single thread.")
            subprocess.check_call(["muscle", "-in", input_file, "-out", aligned_file])
        elif(program == 'mafft'):
            subprocess.check_call(['mafft', '--auto', '--thread', str(threads), input_file], stdout=open(aligned_file, 'w'))
        elif program == 'tcoffee':
            # n_core is deprecated; use -thread <N> in newer versions of T-Coffee
            subprocess.check_call(["t_coffee", "-in", input_file, "-output", "fasta", "-outfile", aligned_file, "-mode", "procoffee", "-n_core", str(threads)])
        else: # clustalo (forcing overwriting)
            subprocess.check_call(["clustalo", "-i", input_file, "-o", aligned_file, "--threads", str(threads), "--force"])
        print(program.upper(), "completed. The aligned sequences were saved in", aligned_file)
    except OSError as e:
        print("ERROR: Make sure that", program.upper(), "is installed and in your PATH.")
        sys.exit(1)
    except Exception as e:
        print("Error code:", e)
        sys.exit(2)
else:
    print("The specified program is not supported. Supported programs are: ClustalO, MUSCLE, MAFFT, and T-Coffee.")