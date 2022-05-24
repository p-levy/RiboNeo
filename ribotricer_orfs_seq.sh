#!/bin/bash

# Generate sequence for ORFs in ribotricer's index

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { echo "Usage: \$1=candidate orf file, \$2=genome fasta, \$3=output name"; exit 1; }
fi

set -uex

# Variables
INDEX=$1 #candidate_orfs_translating_final.tsv
FASTA=$2 #GRCh38.primary_assembly.genome.fa
OUTPUT=$3 #candidate_orfs_translating_final_aa_seq.txt

ribotricer orfs-seq --ribotricer_index $INDEX --fasta $FASTA --protein --saveto $OUTPUT

printf DONE
