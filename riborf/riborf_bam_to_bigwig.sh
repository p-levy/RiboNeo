#!/bin/bash

# Converts bam file to a bigwig coverage file

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { echo "Usage: \$1=bam file (offset_corrected.sorted.bam)
\$2=outdir path"; exit 1; }
fi

set -uex

BASENAME=$(basename -s .bam $1)

echo creating bigWig coverage files

bamCoverage -b $1 --filterRNAstrand forward --binSize 1 -o ${2}/${BASENAME}.neg.bigWig \
& bamCoverage -b $1 --filterRNAstrand reverse --binSize 1 -o ${2}/${BASENAME}.pos.bigWig

echo converting to wig

bigwigtowig ${2}/${BASENAME}.neg.bigWig ${2}/${BASENAME}.neg.wig \
& bigwigtowig ${2}/${BASENAME}.pos.bigWig ${2}/${BASENAME}.pos.wig
