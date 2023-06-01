#!/bin/bash

######################################################################################################
# Author: Pierre Levy (levy.pierre@yahoo.fr)
# Date: 27-04-23
# Post-processing of RibORF data: combining with variant analysis (VCF) from whole-exome or whole
# genome sequencing to identify the mutated translating noncanonical ORFs and report their epitope
# aminoacid sequence. 
# Run on server (too long on personal computer)
######################################################################################################

if [ $# -eq 0  ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then { printf "Usage: \$1=VCF file to annotate, \n
\$2=output annotated VCF \n
\$3=Reference genome fasta \n
\$4=custom gff file (derived from traslating ORFs)\n"; exit 1; }
fi

set -uex

####### STEP 5: Run VEP ##############################################################################################################################

vep --input_file ${1} \
    --output_file ${2} \
    --vcf \
    --fasta ${3} \
    --custom ${4},,gff \
    --fork 40 \
    --format vcf
