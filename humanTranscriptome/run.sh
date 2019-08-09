#!/usr/bin/env bash

# Stop on errors.
set -euo pipefail

########################
#      VARIABLES       #
#----------------------#

curr_pwd=$(pwd)
RUNLOG=${curr_pwd}/humanTranscriptome.log
data_dir=${curr_pwd}/data

KSIZE=31

MIN_Q=5
MAX_Q=31
STEP_Q=1

echo "Run by `whoami` on `date`" > $RUNLOG


########################
#      DOWNLOAD        #
#----------------------#


mkdir data
echo "Downloading..."
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz

echo "Preparing Data..."
gunzip gencode.v31.transcripts.fa.gz
python generate_names.py gencode.v31.transcripts.fa
mv gencode.v31.transcripts* data/

########################
#      Processing      #
#----------------------#


FASTA=${data_dir}/gencode.v31.transcripts.fa
NAMES="${FASTA}.names"
IDX="idx_gencode.v31.transcripts"

kCluster index_kmers -n ${NAMES} -f ${FASTA} -k ${KSIZE}

kCluster pairwise --min-q ${MIN_Q} --max-q ${MAX_Q} --step-q ${STEP_Q} -i ${IDX}