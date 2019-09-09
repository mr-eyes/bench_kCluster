#!/usr/bin/env bash
set -euo pipefail

kCluster_pairwise=$(pwd)/../kCluster2/build/kCluster_pairwise

if ! [ -x "$(command -v kCluster)" ]; then
  echo 'Error: kCluster is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v ${kCluster_pairwise})" ]; then
  echo 'Error: kCluster_pairwise is not installed.' >&2
  exit 1
fi

OK="\e[32m[OK] \e[0m"
PROCESSING="\e[33m[RUNNING] \e[0m"

mkdir -p oma_data/maps
mkdir -p oma_seqs

# echo "Downloading gene2refseq.gz ..."
# wget -N ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

# echo "Downloading oma-ref-seq.txt.gz ..."
# wget -N https://omabrowser.org/All/oma-refseq.txt.gz

#######################################
#             DOWNLOAD                #
#######################################

echo -e "\e[33m\e[1mDownloading.. \e[0m"

FILE=./oma_data/eukaryotes.cdna.fa.gz
if [ -f "$FILE" ]; then
    echo -e "${OK} $FILE exist, skipping.."
else 
    echo "Downloading eukaryotes.cdna.fa.gz"
    wget -N https://omabrowser.org/All/eukaryotes.cdna.fa.gz -O ${FILE}
fi


FILE=./oma_data/oma-groups.txt.gz
if [ -f "$FILE" ]; then
    echo -e "${OK} $FILE exist, skipping.."
else 
    echo "Downloading oma-groups.txt.gz"
    wget -N https://omabrowser.org/All/oma-groups.txt.gz -O ${FILE}
fi

#######################################
#            GENERATING MAPS          #
#######################################

echo -e "\e[33m\e[1mGenerating maps .. \e[0m"
FILE=./oma_data/maps/oma_group_to_species.pickle
if [ -f "$FILE" ]; then
    echo -e "${OK} groups_to_species.pickle exists, skipping.."
else
    echo "Generating oma-group-id -> species_list"
    python scripts/generate_maps.py groups_to_species
fi

#######################################
#         Extracting Species          #
#######################################

echo -e "\e[33m\e[1mExtracting Species .. \e[0m"


# Exp_1 : HUMAN GORGO PANTR PANPA PONAB NOMLE
FILE=./oma_seqs/exp_1
if [ -d "$FILE" ]; then
    echo -e "${OK} Exp_1 sequences found, skipping the extraction.."
else
    echo "Extracting [HUMAN GORGO PANTR PANPA PONAB NOMLE]"
    python scripts/filter_oma_groups.py  HUMAN GORGO PANTR PANPA PONAB NOMLE
fi

#######################################
#              Indexing               #
#######################################

echo -e "\e[33m\e[1mIndexing .. \e[0m"

for dir in oma_seqs/*     # list directories in the form "/tmp/dirname/"
do

    exp_id=${dir%*/}
    exp_id=${exp_id##*/}
    
    if ls ${dir}/*map 1> /dev/null 2>&1; then
        echo -e "${OK} ${exp_id} already indexed, skipping.."
    else
        echo -e "${PROCESSING} Indexing ${exp_id##*/} .."
        kCluster index_kmers -f ${dir}/*.fa -n  ${dir}/*.fa.names -k 31
        mv idx* ${dir}
    fi

done

#######################################
#              Pairwise               #
#######################################

echo -e "\e[33m\e[1mCalculating Pairwise similarity .. \e[0m"

# virtualQs
# Set here all virtualQs you want to process
# ---------------------
QsList="30,31"
# ---------------------


for dir in oma_seqs/*     # list directories in the form "/tmp/dirname/"
do

    exp_id=${dir%*/}
    exp_id=${exp_id##*/}
    exp_no=$(echo "${exp_id//[!0-9]/}")
    idx_prefix=${dir}/idx_exp${exp_no}

    if ls ${dir}/*kCluster.tsv 1> /dev/null 2>&1; then
        echo -e "${OK} ${exp_id} pairwise matrix already exist, skipping.."
    else
        echo -e "${PROCESSING} Generating ${exp_id##*/} pairwise TSV .."
        ${kCluster_pairwise} --idx=${idx_prefix} --qs=${QsList}
    fi

done


#######################################
#              Pivoting               #
#######################################

echo -e "\e[33m\e[1mPivoting .. \e[0m"

for dir in oma_seqs/*     # list directories in the form "/tmp/dirname/"
do

    exp_id=${dir%*/}
    exp_id=${exp_id##*/}
    exp_no=$(echo "${exp_id//[!0-9]/}")
    idx_prefix=${dir}/idx_exp${exp_no}

    if ls ${dir}/*pivoted.tsv 1> /dev/null 2>&1; then
        echo -e "${OK} ${exp_id} pivoted pairwise matrix already exist, skipping.."
    else
        echo -e "${PROCESSING} Pivoting ${exp_id##*/} pairwise TSV .."
        ${kCluster_pairwise} pivote --idx=${dir}/idx_exp${exp_no} --qs=${QsList}
    fi

done

#######################################
#            Clustering               #
#######################################

for dir in oma_seqs/*     # list directories in the form "/tmp/dirname/"
do

    exp_id=${dir%*/}
    exp_id=${exp_id##*/}
    exp_no=$(echo "${exp_id//[!0-9]/}")
    idx_prefix=${dir}/idx_exp${exp_no}
    # echo -e "${PROCESSING} Clustering ${dir}/idx_exp${exp_no}_pivoted.tsv"
    echo -e "\e[33m\e[1mClustering ${dir}/idx_exp${exp_no}_pivoted.ts .. \e[0m"

    for THRESHOLD in {0..100..1};
    do
        THRESHOLD=$(printf "%02d" $THRESHOLD)
        FILE=${dir}/clusters_0.${THRESHOLD}%_idx_exp${exp_no}_pivoted.tsv
        if [ -f "$FILE" ]; then
            echo -e "${OK} ${FILE} exists, skipping.."
        else
            echo "Threshold ${THRESHOLD}%"
            kCluster cluster --qs ${QsList} --tsv ${dir}/idx_exp${exp_no}_pivoted.tsv --cutoff 0.${THRESHOLD}
            mv clusters_0.${THRESHOLD}%_idx_exp${exp_no}_pivoted.tsv ${dir}
        fi
    done

done

