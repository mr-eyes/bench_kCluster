"""
input : directory contains summary files |seqs\tclusters\ttype|
summary*txt files must named with t{similarity_theshold} with trailing zeros to be sorted correctly
"""


import re
import sys
import os
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline
from glob import glob

summary_dir = sys.argv[1]
output_html = os.path.join(os.path.dirname(summary_dir),"visualization.html")

thresholds = []
types = {}
counts = {}
#kmers_count = os.popen('grep "" kmers_clustering/full_human_transcriptome_results/*/*summary*').read()

#for line in kmers_count.split("\n"):

for summary_file in sorted(glob(summary_dir + "/*sum*txt")):
    f = open(summary_file, 'r')
    next(f)  # skip headers
    _threshold = re.findall(r"_(\d+\.\d+)%", summary_file)[0].replace(".","")

    thresholds.append(_threshold)

    for line in f:
        #print line
        line = line.split()
        seq_count = line[0]
        clstrs_count = line[1]
        _type = line[2]
        #print _type

        if _type not in types:
            types[_type] = [clstrs_count]
            counts[_type] = [seq_count]
        else:
            types[_type].append(clstrs_count)
            counts[_type].append(seq_count)
    f.close()

# print thresholds

# print types
# print counts

complete_clean = go.Scatter(
    x=thresholds,
    y=counts["cc"],
    mode='lines+markers',
    name='Complete Clean'
)

complete_mixed = go.Scatter(
    x=thresholds,
    y=counts["cm"],
    mode='lines+markers',
    name='Complete Mixed'
)

incomplete_clean = go.Scatter(
    x=thresholds,
    y=counts["ic"],
    mode='lines+markers',
    name='Incomplete Clean'
)

incomplete_mixed = go.Scatter(
    x=thresholds,
    y=counts["im"],
    mode='lines+markers',
    name='Incomplete Mixed'
)

data = [complete_mixed, complete_clean, incomplete_mixed, incomplete_clean]

layout = dict(title='<b>CD-HIT Clustering Assessment</b>',
              xaxis=dict(title='<b>Threshold</b>', type='category'),
              yaxis=dict(title='<b>Number Of Transcripts</b>', type="log"),
              )

fig = dict(data=data, layout=layout)

plotly.offline.plot(fig, filename=output_html, auto_open=False)
