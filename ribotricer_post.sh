#!/bin/bash

# Script for processing ribotricer output (after reformat step) to compute non-canonical specific read counts

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { printf "\n Usage: \n \$1=reformatted unique translating orfs \
\n \$2=prefix \n \$3=pos wig \n \$4=neg wig \n\n"; exit 1; }
fi

$(dirname $0)/create_bed.R $1 $2

bedtools subtract -b $2_annotated_orfs.bed -a $2_non_annotated_orfs.bed -s > $2_novel.bed

$(dirname $0)/orf_count.R $2_novel.bed $3 $4 $1