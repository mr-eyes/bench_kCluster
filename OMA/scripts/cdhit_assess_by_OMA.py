from __future__ import division
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


def _mean(lst):
    return round(sum(lst) / len(lst),2)


def _std(lst):
    mean = sum(lst) / len(lst)   # mean
    var = sum(pow(x-mean, 2) for x in lst) / len(lst)  # variance
    std = math.sqrt(var)  # standard deviation
    return round(std, 2)

def _median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

fasta_file_path = ""
clstr_file_path = ""
output_file = ""

if len(sys.argv) < 3:
    sys.exit(
        "Kindly pass positional arguments, ex: python clusters_assessment.py <fasta_file> <clstr_file> <output_file>")

else:
    fasta_file_path = sys.argv[1]
    clstr_file_path = sys.argv[2]

output_dir = os.path.dirname(clstr_file_path)
output_prefix = os.path.basename(clstr_file_path).replace(".cdhit.clstr","")


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


with open(fasta_file_path) as fa:
    for line in fa:
        if line[0] != ">":
            continue
        fields = line.split("|")
        transcript_id = fields[0][1:]
        gene = fields[1]
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


clstr_file = open(clstr_file_path, "r")
clstr_data = clstr_file.read()
clstr_file.close()

rep = {"\t": ",", "at +/": "", "at -/": "",
       "...": ",", "nt": "", "%": "", " ": ""}
rep = dict((re.escape(k), v) for k, v in rep.items())
pattern = re.compile("|".join(rep.keys()))
clstr_data = pattern.sub(lambda m: rep[re.escape(m.group(0))], clstr_data)
all_clusters = clstr_data.split(">Cluster")

clusters_transcripts_ids = {}

for i in range(1, len(all_clusters), 1):
    cluster = all_clusters[i]
    cluster = cluster.split("\n")
    cluster_id = int(cluster[0])

    for item in cluster[1:-1]:
        item = item.replace(">", "").split(",")
        transcript_id = item[2].split("|")[0]
        gene_id = item[2].split("|")[1]

        if cluster_id in clusters_transcripts_ids:
            clusters_transcripts_ids[cluster_id].append(transcript_id)

        else:
            clusters_transcripts_ids[cluster_id] = [transcript_id]

output_file = os.path.join(output_dir, "details", "assessement_" +output_prefix + ".tsv")
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
        build_stats("_complete_mixed", no_loci, no_genes,
                    no_complete_genes, no_complete_loci)
        _cm_transcripts_count += len(transcripts_ids)
    if q1 == True and q2 == False:
        ans1, ans2 = "Complete", "Clean"
        _complete_clean += 1
        build_stats("_complete_clean", no_loci, no_genes,
                    no_complete_genes, no_complete_loci)
        _cc_transcripts_count += len(transcripts_ids)
    if q1 == False and q2 == True:
        ans1, ans2 = "InComplete", "Mixed"
        _incomplete_mixed += 1
        build_stats("_incomplete_mixed", no_loci, no_genes,
                    no_complete_genes, no_complete_loci)
        _im_transcripts_count += len(transcripts_ids)
    if q1 == False and q2 == False:
        ans1, ans2 = "InComplete", "Clean"
        _incomplete_clean += 1
        build_stats("_incomplete_clean", no_loci, no_genes,
                    no_complete_genes, no_complete_loci)
        _ic_transcripts_count += len(transcripts_ids)

    line = str(cluster_id) + "\t" + ans1 + "\t" + ans2 + "\t" + str(no_loci) + "\t" + str(no_complete_loci) + "\n"
    res.write(line)

res.close()


# Writing summary file of counts _______________________________
tmp_prefix = output_prefix.split("_")
tmp_prefix = tmp_prefix[0] + "_0." + tmp_prefix[1] + "%"
output_file = os.path.join(output_dir, "summaries", "assessement_" + tmp_prefix + "_summary.txt")
summary = open(output_file, 'w')
summary.write("seqs\tclstrs\ttype\n")
summary.write(("%d\t%d\tcm\n") % (_cm_transcripts_count, _complete_mixed))
summary.write(("%d\t%d\tcc\n") % (_cc_transcripts_count, _complete_clean))
summary.write(("%d\t%d\tim\n") % (_im_transcripts_count, _incomplete_mixed))
summary.write(("%d\t%d\tic\n") % (_ic_transcripts_count, _incomplete_clean))
summary.close()



##### DISABLED

# Writing statistics json file ________________________________

# json_output = {}
# for cluster_type in ["_complete_mixed","_complete_clean","_incomplete_mixed","_incomplete_clean"]:
#     result =  {
#         'mean': {
#             'no_genes': _mean(stats["genes"][cluster_type]),
#             'complete_genes':  _mean(stats["complete-genes"][cluster_type]),
#             'no_loci':  _mean(stats["loci"][cluster_type]),
#             'complete_loci':  _mean(stats["complete-loci"][cluster_type])
#         },
#         'std': {
#             'no_genes': _std(stats["genes"][cluster_type]),
#             'complete_genes':  _std(stats["complete-genes"][cluster_type]),
#             'no_loci':  _std(stats["loci"][cluster_type]),
#             'complete_loci':  _std(stats["complete-loci"][cluster_type])
#         },
#         'min': {"no_genes": min(stats["genes"][cluster_type]),
#                 "complete_genes": min(stats["complete-genes"][cluster_type]),
#                 "no_loci": min(stats["loci"][cluster_type]),
#                 "complete_loci": min(stats["complete-loci"][cluster_type])},
#         'max': {"no_genes": max(stats["genes"][cluster_type]),
#                 "complete_genes": max(stats["complete-genes"][cluster_type]),
#                 "no_loci": max(stats["loci"][cluster_type]),
#                 "complete_loci": max(stats["complete-loci"][cluster_type])},
#         'median': {"no_genes": _median(stats["genes"][cluster_type]),
#                    "complete_genes": _median(stats["complete-genes"][cluster_type]),
#                    "no_loci": _median(stats["loci"][cluster_type]),
#                    "complete_loci": _median(stats["complete-loci"][cluster_type])}
#     }
#     json_output[cluster_type] = result


# json_file = open("".join(output_file.split(".")[:-2])[1:] + "_stats.json", "w")
# json_file.write(json.dumps(json_output, sort_keys=True,
#                            indent=4, separators=(',', ': ')))
# json_file.close()
