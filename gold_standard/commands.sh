#Indexing
FASTA="seq.fa"
NAMES="${FASTA}.names"

#Pairwise
IDX="${FASTA##*/}"
IDX="idx_${IDX%.*}"

KSIZE=25
MIN_Q=5
MAX_Q=${KSIZE}
STEP_Q=5



rm -rf *sqlite
kCluster index_kmers -n ${NAMES} -f ${FASTA} -k ${KSIZE}
kCluster pairwise --min-q ${MIN_Q} --max-q ${MAX_Q} --step-q ${STEP_Q} -i ${IDX} # Does not print anything
kCluster dump --db ${IDX}_kCluster.sqlite > sqlite_dump.tsv
kCluster cluster -m ${MIN_Q} -M ${MAX_Q} -s ${STEP_Q} --db ${IDX}_kCluster.sqlite --cutoff 0.0
kCluster cluster -m ${MIN_Q} -M ${MAX_Q} -s ${STEP_Q} --db ${IDX}_kCluster.sqlite --cutoff 0.9