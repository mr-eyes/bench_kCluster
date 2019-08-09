import sys
import os

input_file = sys.argv[1]
names_file = os.path.basename(input_file) + ".names"
names = open(names_file,"w")

with open (input_file, 'r') as fasta_sequences, open(names_file,"w") as names:
    for fasta in fasta_sequences:
        if fasta[0] == ">":
            header = fasta.strip()[1:]
            names.write(header + "\t" + header + "\n")