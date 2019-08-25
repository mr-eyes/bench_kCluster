#!/usr/bin/env bash
set -euo pipefail

if ! [ -x "$(command -v jellyfish)" ]; then
  echo 'Error: Jellyfish is not installed.' >&2
  exit 1
fi

if [ "$#" -ne 3 ]; then
    echo "run: ./shared_kmers.sh 31 seq1.fa seq2.fa "
    exit 2
fi

KSIZE=$1
SEQ1=$2
SEQ2=$3

echo "Processing ..."
jellyfish count -m ${KSIZE} -s 100M -t 10 -C ${SEQ1} &> /dev/null
jellyfish dump -t -c mer_counts.jf | sort | uniq | sed '/^[[:blank:]]*$/d' | awk -F '\t' '{print $1}' > seq1.kmers
jellyfish count -m ${KSIZE} -s 100M -t 10 -C ${SEQ2} &> /dev/null
jellyfish dump -t -c mer_counts.jf | sort | uniq | sed '/^[[:blank:]]*$/d' | awk -F '\t' '{print $1}' > seq2.kmers

SHARED=$(awk 'a[$0]++' seq1.kmers seq2.kmers | wc -l)
# rm mer_counts.jf seq1.kmers seq2.kmers

echo "${SEQ1} & ${SEQ2} have ${SHARED} shared unique kmers"