#!/bin/bash

# Script for processing ribotricer output (after reformat step) to compute non-canonical specific read counts

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { printf "\n Usage: \n \$1=candidate orfs \n \$2=reformatted unique translating orfs \
\n \$3=prefix \n \$4=pos wig \n \$5=neg wig \n \$6=threads \n\n"; exit 1; }
fi

$(dirname $0)/create_bed.R $1 $2 $3

bedtools subtract -b $3_annotated_orfs.bed -a $3_non_annotated_orfs.bed -s > $3_novel.bed

$(dirname $0)/orf_count.R $3_novel.bed $4 $5 $2 $3 $6
