"""
This script for generating a pickle file of {refseqID : taxID}
"""

import gzip
import sys
import pickle

mapped = dict()

with gzip.open("gene2refseq.gz", 'rt') as gene2refseq:
    next(gene2refseq)
    for line in gene2refseq:
        line = line.strip().split()
        if line[5] != "-":
            mapped[line[5]] =  int(line[0])


print("writing pickle...")

with open('refseq_to_tax.pickle', 'wb') as handle:
    pickle.dump(mapped, handle, protocol=pickle.HIGHEST_PROTOCOL)

