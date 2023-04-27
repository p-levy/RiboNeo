#!/bin/bash

# Generate candidate ORF list filtered for translating ORFs only
# reduces the size of the candidate ORFs, easier to manipulate later in R for e.g.
# Run on server (too long on personal computer)

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { echo "Usage: \$1=translating orf file (repre.valid.pred.pvalue.parameters.txt), \
\$2=outdir, \
\$3=candidate orf genepred"; exit 1; }
fi

set -uex

# Isolated list of ORF_IDs from translating_ORFs file: translating_ORF_ID.txt
cat $1 | cut -f 1 | grep -v orfID > $2/translating_ORF_ID.txt

# Selected candidate ORFs only translating: candidate_orfs_translating.tsv
# Using join command (need to sort the id file and table file first)

join -t$'\t' <(sort -d $2/translating_ORF_ID.txt) <(sort -d -t$'\t' -k1,1 $3) > $2/translating.candidateORF.genepred.txt


# remove intermediate files (Keep for filtering of fasta with seqkit grep)
# rm -rf $2/translating_ORF_ID.txt

