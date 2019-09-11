#!/usr/bin/env bash
set -euo pipefail

PROCESSING="\e[33m[RUNNING] \e[0m"
OK="\e[32m[OK] \e[0m"

if ! [ -x "$(command -v kCluster)" ]; then
  echo 'kCluster is not installed.' >&2
  echo -e "${PROCESSING} installing .." >&2
  pip install git+https://github.com/mr-eyes/kCluster#egg=kCluster

else
    echo -e "${OK} kCluster exist, skipping.."
fi


FILE=./kCluster2/build
if [ -d "$FILE" ]; then
    echo -e "${OK} kCluster2 (pairwise) exist, skipping.."
else 
    git clone https://github.com/mr-eyes/kCluster2.git
    cd kCluster2
    git submodule update --init --recursive
    mkdir build && cd build
    cmake ..
    make
fi

echo -e "---------"
echo -e "${OK} Everthing has been installed successfully"
echo ""
echo "\e[33m\e[1mTest by invoking: ./kCluster2/build/kCluster_pairwise"