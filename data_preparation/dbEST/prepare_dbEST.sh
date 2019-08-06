#!/usr/bin/env bash

# Stop on errors.
set -euo pipefail

########################
#      VARIABLES       #
#----------------------#

curr_pwd=$(pwd)
RUNLOG=${curr_pwd}/dbEST_data_preparation.log
dbEST_dir=${curr_pwd}/dbEST_reports
dbEST_fasta=${curr_pwd}/dbEST_fasta



echo "Run by `whoami` on `date`" > $RUNLOG


########################
#      DOWNLOAD        #
#----------------------#

# mkdir ${dbEST_dir}
# mkdir ${dbEST_fasta}

cd ${dbEST_dir}

for  i  in {1..82};
    do 
        echo "downloading https://ftp.ncbi.nih.gov/repository/dbEST/dbEST.reports.000000.${i}.gz"
        wget -c https://ftp.ncbi.nih.gov/repository/dbEST/dbEST.reports.000000.${i}.gz >>${RUNLOG} 2>&1;
done;



########################
#   Convert to FASTA   #
#----------------------#

cd ..
cd ${dbEST_fasta}

FILTER="homo sapien"

for  i  in {1..82};
    do 
        echo "processing dbEST.reports.000000.${i}.gz"
        python ../dbEST_report_to_fasta.py ${dbEST_dir}/dbEST.reports.000000.${i}.gz ${FILTER} >>${RUNLOG} 2>&1;
done;
