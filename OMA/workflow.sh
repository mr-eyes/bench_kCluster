#!/usr/bin/env bash
set -euo pipefail


OK="\e[32m[OK] \e[0m"

mkdir -p oma_data/maps

# echo "Downloading gene2refseq.gz ..."
# wget -N ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

# echo "Downloading oma-ref-seq.txt.gz ..."
# wget -N https://omabrowser.org/All/oma-refseq.txt.gz

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


# echo "Extracting Homo Sapiens & Macaca Mulatta"
# zcat eukaryotes.cdna.fa.gz | seqkit grep -r -p HUMAN -p MACMU > human_mulatta.fa

echo -e "\e[33m\e[1mGenerating maps .. \e[0m"
# Choose from [groups_to_species, groups_to_refseq]

FILE=./oma_data/maps/oma_group_to_species.pickle
if [ -f "$FILE" ]; then
    echo -e "${OK} groups_to_species.pickle exists, skipping.."
else
    echo "Generating oma-group-id -> species_list"
    python scripts/generate_maps.py groups_to_species
fi
