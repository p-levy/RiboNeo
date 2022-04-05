#!/bin/bash

$(dirname $0)/create_bed.R $1 $2

bedtools subtract -b $2_annotated_orfs.bed -a $2_non_annotated_orfs.bed -s > $2_novel.bed

$(dirname $0)/orf_count.R $2_novel.bed $3 $4 $1