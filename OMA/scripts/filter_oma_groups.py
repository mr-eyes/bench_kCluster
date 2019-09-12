import sys
import os
import pickle
import re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

"""
1- Parse system argyments for species of interest.
2- Load groups_to_species.pickle
3- Load and Parse Eukorayotes fasta file, Generate new folder with fasta + names
"""

BLOSUM = {
    "A" : "A",
    "C" : "C",
    "G" : "G",
    "T" : "T",
    "U" : "T",
    "W" : "A",
    "S" : "G",
    "M" : "A",
    "K" : "G",
    "R" : "A",
    "Y" : "C",
    "B" : "G",
    "D" : "A",
    "H" : "A",
    "V" : "A",
    "N" : "A",
    "Z" : "",
    "X" : "A",
}

###########
#   (1)   #
###########

SPECIES = list()
oma_groups_path = str()

if len(sys.argv) > 1:
    SPECIES = set([x.upper() for x in sys.argv[1:]])
    SPECIES_NO = len(sys.argv) - 1

else:
    print("Failed in parsing arguments..")
    sys.exit(1)

###########
#   (2)   #
###########

ABS_PATH = os.path.abspath(os.path.split(sys.argv[0])[0])
pickle_path = os.path.join(ABS_PATH, "..", "oma_data", "maps", "oma_group_to_species.pickle")
print("loading oma_group_to_species.pickle")
with open(pickle_path, 'rb') as handle:
    all_group_to_species = pickle.load(handle)

protein_to_oma_group = dict()
for oma_group, species in all_group_to_species.items():
    temp_proteins = list()    

    for entry in species:
        if re.findall("[A-Z]*", entry)[0] in SPECIES:
            temp_proteins.append(entry)
    
    else:
        if len(temp_proteins): # == SPECIES_NO:
            for prot in temp_proteins:
                protein_to_oma_group[prot] = oma_group

            temp_proteins.clear()

protein_targets_keys = list(protein_to_oma_group.keys())

###########
#   (3)   #
###########

print("Processing...")

output_dir = os.path.join(ABS_PATH, ".." , "oma_seqs") # , "merged_" + "_".join([x.lower() for x in SPECIES]))
exp_id = sum(os.path.isdir(os.path.join(output_dir, i)) for i in os.listdir(output_dir)) + 1
exp_id = str(exp_id)
exp_dir = os.path.join(output_dir, "exp_" + exp_id)

if not os.path.exists(exp_dir):
    os.makedirs(exp_dir)
else:
    print("unable to create output dir")
    sys.exit(1)

with open(os.path.join(exp_dir, "info.txt"), 'w') as info:
    for entry in SPECIES:
        info.write(entry + "\n")


cdna_fasta = os.path.join(ABS_PATH, ".." , "oma_data", "eukaryotes.cdna.fa.gz")
fasta_sequences = SeqIO.parse(gzip.open(cdna_fasta,'rt'),'fasta')

fasta_output = os.path.join(exp_dir, "exp" + exp_id + ".fa")
names_output = fasta_output + ".names"

group_to_names = dict()

with open(fasta_output, 'w') as fastaFile, open(names_output, 'w') as namesFile:
    for seq in fasta_sequences:
        protein_name = str(seq.description).strip()

        if protein_name in protein_to_oma_group:
            seq.description = seq.id = seq.name = protein_name + "|" + str(protein_to_oma_group[protein_name])
        else:
            continue

        if re.findall("[A-Z]*", protein_name)[0] in SPECIES:
            new_seq = str()
            for ch in seq.seq:
                new_seq += BLOSUM[ch]
            seq.seq = Seq(new_seq)
            SeqIO.write(seq, fastaFile, "fasta")
            namesFile.write(seq.description + '\t' + seq.description + '\n')

            if protein_name in protein_targets_keys:
                group_id = protein_to_oma_group[protein_name]

                if group_id in group_to_names:
                    group_to_names[group_id].append(protein_name)
                else:
                    group_to_names[group_id] = [protein_name]

# with open(names_output, 'w') as namesFile:
#     for group, proteins in group_to_names.items():
#         if len(proteins) == SPECIES_NO:
#             for protein_name in proteins:
#                 namesFile.write(protein_name + "\t" + protein_name + "\n")