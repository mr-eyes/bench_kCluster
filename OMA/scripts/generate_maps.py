"""
This script for generating a pickle file of {refseqID : taxID}
"""

import gzip
import sys
import pickle
import re
import os
import os,sys


class MAPPER:

    ABS_PATH = os.path.abspath(os.path.split(sys.argv[0])[0])

    #######################################
    #  OMA-group-id : [all-species-names] #
    #######################################

    @staticmethod
    def groups_to_species():
        # OPTION : groups_to_species
        print("Mapping oma_group to species names...", file=sys.stderr)
        group_to_species = dict()
        filename= os.path.join(MAPPER.ABS_PATH, "..", "oma_data", "oma-groups.txt.gz")
        with gzip.open(filename, 'rt') as oma_groups:
            for i in range(3): next(oma_groups)
            
            for line in oma_groups:
                line = line.strip().split()
                oma_group = int(line[0])
                # proteins = [re.findall("[A-Z]*", entry)[0] for entry in line[1:] if entry != "n/a"]
                proteins = [entry for entry in line[1:] if entry != "n/a"]
                group_to_species[oma_group] = proteins

        output_pickle = os.path.join(MAPPER.ABS_PATH, "..", "oma_data", "maps", "oma_group_to_species.pickle")
        with open(output_pickle, 'wb') as handle:
            pickle.dump(group_to_species, handle, protocol=pickle.HIGHEST_PROTOCOL)

        group_to_species.clear()


    
    @staticmethod
    def oma_to_refseq():
        # OPTION : groups_to_refseq
        print("Mapping oma_group to refseq...", file=sys.stderr)
        groups_to_refseq = dict()
        filename= os.path.join(MAPPER.ABS_PATH, "..", "oma_data", "oma-refseq.txt.gz")
        with gzip.open(filename, 'rt') as oma_refseq:
            next(oma_refseq); next(oma_refseq)

            for line in oma_refseq:
                line = line.strip().split()
                groups_to_refseq[line[0]] =  line[1]


        output_pickle = os.path.join(MAPPER.ABS_PATH, "..", "oma_data", "maps", "oma_to_refseq.pickle")
        with open(output_pickle, 'wb') as handle:
            pickle.dump(groups_to_refseq, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        groups_to_refseq.clear()





#######################################|
#######################################|



# mapped = dict()
# with gzip.open("gene2refseq.gz", 'rt') as gene2refseq:
#     next(gene2refseq)
#     for line in gene2refseq:
#         line = line.strip().split()
#         if line[5] != "-":
#             mapped[line[5]] =  int(line[1])

# print("writing pickle...")

# with open('refseq_to_gene.pickle', 'wb') as handle:
#     pickle.dump(mapped, handle, protocol=pickle.HIGHEST_PROTOCOL)

# mapped.clear()


######################################
#  OMA-group : [all-species-names]   #
######################################





if __name__ == "__main__":
    MAP = {
        "groups_to_species" : MAPPER.groups_to_species,
        "groups_to_refseq" : MAPPER.oma_to_refseq
    }
    
    if len(sys.argv) > 1:
        if sys.argv[1] in MAP:
            MAP[sys.argv[1]]()

    else:
        print("Please select from: ", str(", ".join(list(MAP.keys()))))