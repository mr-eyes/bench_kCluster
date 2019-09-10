from __future__ import division
import tqdm
import re
import sys
import json
import math
import os


def transcripts_to_loci(transcipts_ids):
    loci = {}  # {locus:counted}
    for transcript in transcipts_ids:
        locus = transcript_locus[transcript]

        if locus in loci:
            loci[locus] += 1
        else:
            loci[locus] = 1

    return loci


def transcripts_to_genes(transcipts_ids):
    genes = {}  # {locus:counted}
    for transcript in transcipts_ids:
        gene = transcript_gene[transcript]

        if gene in genes:
            genes[gene] += 1
        else:
            genes[gene] = 1

    return genes


def how_many_loci(cluster):
    return len(transcripts_to_loci(cluster))


def how_many_genes(cluster):
    return len(transcripts_to_genes(cluster))


def how_many_complete_loci(cluster):
    loci = transcripts_to_loci(cluster)
    complete = 0
    for key, value in loci.items():
        if len(locus_transcripts[key]) == value:
            complete += 1

    return complete


def how_many_complete_genes(cluster):
    genes = transcripts_to_genes(cluster)
    complete = 0
    for key, value in genes.items():
        if len(gene_transcripts[key]) == value:
            complete += 1

    return complete


def Q1(cluster):
    cluster_genes = transcripts_to_genes(cluster)
    for key, value in cluster_genes.items():
        if len(gene_transcripts[key]) != value:
            return False

    return True


def Q2(cluster):
    if len(transcripts_to_loci(cluster)) > 1:
        return True
    else:
        return False


def build_stats(cluster_type, no_loci, no_genes, no_complete_genes, no_complete_loci):
    stats["loci"][cluster_type].append(no_loci)
    stats["genes"][cluster_type].append(no_genes)
    stats["complete-genes"][cluster_type].append(no_complete_genes)
    stats["complete-loci"][cluster_type].append(no_complete_loci)


fasta_file_path = ""
clstr_file_path = ""
output_file = ""

if len(sys.argv) < 2:
    sys.exit(
        "Kindly pass positional arguments, ex: python clusters_assessment.py [clstr_file]")

else:
    clstr_file_path = sys.argv[1]
    output_file = os.path.dirname(clstr_file_path) + '/' + "assessement_" + os.path.basename(clstr_file_path)

_complete_mixed = 0
_complete_clean = 0
_incomplete_mixed = 0
_incomplete_clean = 0

stats = {"loci": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "genes": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "complete-genes": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "complete-loci": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []}
         }


locus_transcripts = {}  # locus_0 : gene1,gene2,...
transcript_locus = {}  # gene1:locus1, gene2:locus1, gene3:locus3

gene_transcripts = {}  # gene_id: transcript1,transcript2,.....
transcript_gene = {}   # transcript1: gene9, transcript2: gene1, ....


with open(clstr_file_path) as fa:
    next(fa) # skip header
    for line in fa:
        fields = line.strip().split('\t')[1].split(',')
        for field in fields:
            transcript_id = field.split("|")[0]
            gene = field.split("|")[1]
            locus = gene
            transcript_locus[transcript_id] = locus
            transcript_gene[transcript_id] = gene

            if gene in gene_transcripts:
                gene_transcripts[gene].append(transcript_id)
            else:
                gene_transcripts[gene] = [transcript_id]

            if locus in locus_transcripts:
                locus_transcripts[locus].append(transcript_id)
            else:
                locus_transcripts[locus] = [transcript_id]

clusters_transcripts_ids = {}

with open(clstr_file_path, "r") as clusters:
    next(clusters)
    for line in clusters:
        fields = line.strip().split('\t')
        cluster_id = int(fields[0])
        _transcripts = fields[1].split(",")
        transcripts = []
        for tr in _transcripts:
            transcripts.append(tr.split("|")[0])

        clusters_transcripts_ids[cluster_id] = transcripts

res = open(output_file, "w")
res.write("cluster_id\tQ1\tQ2\tOMA_groups\tcomplete_OMA_groups\n")


_cc_transcripts_count = 0
_ic_transcripts_count = 0
_cm_transcripts_count = 0
_im_transcripts_count = 0

for cluster_id, transcripts_ids in sorted(clusters_transcripts_ids.items()):
    q1 = Q1(transcripts_ids)
    q2 = Q2(transcripts_ids)
    no_loci = how_many_loci(transcripts_ids)
    no_genes = how_many_genes(transcripts_ids)
    no_complete_loci = how_many_complete_loci(transcripts_ids)
    no_complete_genes = how_many_complete_genes(transcripts_ids)

    ans1 = ""
    ans2 = ""

    if q1 == True and q2 == True:
        ans1, ans2 = "Complete", "Mixed"
        _complete_mixed += 1
        build_stats("_complete_mixed",no_loci, no_genes, no_complete_genes, no_complete_loci)
        _cm_transcripts_count += len(transcripts_ids)
    if q1 == True and q2 == False:
        ans1, ans2 = "Complete", "Clean"
        _complete_clean += 1
        build_stats("_complete_clean",no_loci, no_genes, no_complete_genes, no_complete_loci)
        _cc_transcripts_count += len(transcripts_ids)
    if q1 == False and q2 == True:
        ans1, ans2 = "InComplete", "Mixed"
        _incomplete_mixed += 1
        build_stats("_incomplete_mixed",no_loci, no_genes, no_complete_genes, no_complete_loci)
        _im_transcripts_count += len(transcripts_ids)
    if q1 == False and q2 == False:
        ans1, ans2 = "InComplete", "Clean"
        _incomplete_clean += 1
        build_stats("_incomplete_clean",no_loci, no_genes, no_complete_genes, no_complete_loci)
        _ic_transcripts_count += len(transcripts_ids)

    line = str(cluster_id) + "\t" + ans1 + "\t" + ans2 + "\t" + str(no_loci) + "\t" + \
        str(no_complete_loci) + "\t" + str(no_genes) + \
        "\t" + str(no_complete_genes) + "\n"

    line = f"{cluster_id}\t{ans1}\t{ans2}\t{no_loci}\t{no_complete_loci}\n"

    res.write(line)

res.close()

# Writing summary file of counts ________________________________
summary = open(output_file.replace(".tsv","") + "_summary.txt", "w")
summary.write("seqs\tclstrs\ttype\n")
summary.write(("%d\t%d\tcm\n") % (_cm_transcripts_count, _complete_mixed))
summary.write(("%d\t%d\tcc\n") % (_cc_transcripts_count, _complete_clean))
summary.write(("%d\t%d\tim\n") % (_im_transcripts_count, _incomplete_mixed))
summary.write(("%d\t%d\tic\n") % (_ic_transcripts_count, _incomplete_clean))
summary.close()


# Disabled
# Writing statistics json file ________________________________

# json_output = {}

# for cluster_type in ["_incomplete_mixed", "_incomplete_clean", "_complete_mixed", "_complete_clean"]:

#     result = {
#             'no_genes': stats["genes"][cluster_type],
#             'complete_genes':  stats["complete-genes"][cluster_type],
#             'no_loci':  stats["loci"][cluster_type],
#             'complete_loci':  stats["complete-loci"][cluster_type]
#         }
#     json_output[cluster_type] = result


# json_file = open(output_file.split(".")[0] + "_stats.json", "w")

# json_file.write(json.dumps(json_output, sort_keys=True,
#                            indent=4, separators=(',', ': ')))

# json_file.close()
