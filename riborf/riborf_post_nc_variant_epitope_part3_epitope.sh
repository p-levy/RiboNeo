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
then { printf "Usage: \$1=VCF file (VEP annotated with RibORF-derived GFF), \n
\$2=outdir\n"; exit 1; }
fi

set -uex

####### STEP 6: Run riborf_vcfpy.py to 1) filter according to VAF, num callers etc. and ##############################################################
####### 2) extract relevant info from annotated VCF to prepare for epitope computation ###############################################################

python /RiboNeo/riborf/variant_epitope/riborf_vcfpy.py -o ${2} ${1}
# output = {outdir}/output_variants.tsv

####### STEP 7: Run riborf_nc_variants.R ##############################################################################################################

Rscript /RiboNeo/riborf/variant_epitope/riborf_nc_variants.R \
    ${2}/output_variants.tsv \
    ${2}/translating.candidateORF.genepred.txt \
    ${2}
# output = {outdir}/non_canonical_variants.tsv

####### STEP 8: Create dictionnaries ###################################################################################################################

python /RiboNeo/riborf/variant_epitope/riborf_fasta_to_aminoacid.py ${2}
python /RiboNeo/riborf/variant_epitope/riborf_fasta_to_dict.py ${2}

####### STEP 9: Create epitopes ########################################################################################################################

python /RiboNeo/riborf/variant_epitope/riborf_create_epitopes.py \
    ${2}/non_canonical_variants.tsv \
    -c ${2}/translating_candidateORF.dict \
    -d /RiboNeo/refs/Homo_sapiens.GRCh38.cdna.all.json \
    -a ${2}/translating_candidateORF_AA.dict \
    -o ${2}
# output = {outdir}/epitopes.tsv