import pickle
import sys
from Bio import SeqIO
import re
import gzip

TAX_IDS = [9606, 9544] # Homo Sapien, Macaca mulatta
SPECIES = ["HUMAN", "MACMU"]

print("loading refseq_to_tax.pickle", file = sys.stderr)
with open('refseq_to_tax.pickle', 'rb') as handle:
    refseq_to_tax = pickle.load(handle)


print("loading oma_to_refseq.pickle", file = sys.stderr)
with open('oma_to_refseq.pickle', 'rb') as handle:
    oma_to_refseq = pickle.load(handle)


cdna_fasta = "eukaryotes.cdna.fa.gz"
fasta_sequences = SeqIO.parse(gzip.open(cdna_fasta,'rt'),'fasta')

output_file = "merged_" + "-".join(list(map(str, TAX_IDS))) + ".fa"

skip1 = 0
skip2 = 0
non_acgts = 0

with open(output_file, 'w') as merged_fasta, open(output_file + ".names", 'w') as namesFile:
    for seq in fasta_sequences:
        is_DNA = re.match("^[ACGT]*$", str(seq.seq)) is not None

        oma_id = seq.description.replace(' ','')
        
        if not is_DNA:
            # print(f"{seq.description.strip()} contains non-ACGT characters")
            non_acgts += 1
            continue
        
        for species in SPECIES:
            if species in oma_id:
                if oma_id not in oma_to_refseq:
                    skip1 += 1
                    continue

                refseq_id = oma_to_refseq[oma_id]

                if refseq_id not in refseq_to_tax:
                    skip2+=1
                    continue

                namesFile.write(refseq_id + "\t" + refseq_id + "\n")
                SeqIO.write(seq, merged_fasta, "fasta")

