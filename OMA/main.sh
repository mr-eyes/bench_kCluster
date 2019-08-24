#!/usr/bin/env bash
set -euo pipefail

echo "Downloading gene2refseq.gz ..."
wget -N ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz

echo "Downloading eukaryotes.cdna.fa.gz"
wget -N https://omabrowser.org/All/eukaryotes.cdna.fa.gz


FILE=./refseq_to_tax.pickle
if test -f "$FILE"; then
    echo "Pickle file already exist, skipping.."
else 
    echo "Generating refseq_id -> tax_id pickle file"
    python extract_tax_refseq.py
fi

