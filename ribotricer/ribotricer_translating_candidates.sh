#!/bin/bash

# Generate candidate ORF list filtered for translating ORFs only
# reduces the size of the candidate ORFs, easier to manipulate later in R for e.g.

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { echo "Usage: \$1=translating orf file, \$2=candidate orfs file, \$3=outdir"; exit 1; }
fi

set -uex

# Isolated list of ORF_IDs from translating_ORFs file: translating_ORF_ID.txt
cat $1 | cut -f 1 | grep -v ORF_ID > $3/translating_ORF_ID.txt

# Selected candidate ORFs only translating: candidate_orfs_translating.tsv
# Using join command (need to sort the id file and table file first)
sort -d  $3/translating_ORF_ID.txt > $3/sorted_translating_ids.txt
sort -d -t$'\t' -k1,1 $2 > $3/sorted_candidate_orfs.tsv
join -t $'\t' $3/sorted_translating_ids.txt $3/sorted_candidate_orfs.tsv > $3/candidate_orfs_translating.tsv

# Add back first line (column header)
cat $2 | head -1 | cat - $3/candidate_orfs_translating.tsv > $3/candidate_orfs_translating_final.tsv
	# the - flag tells cat to read from stdin

# remove intermediate files
rm -rf $3/translating_ORF_ID.txt sorted_translating_ids.txt $3/sorted_candidate_orfs.tsv $3/candidate_orfs_translating.tsv
