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
then { printf "Usage: \$1=translating orf file (repre.valid.pred.pvalue.parameters.txt), \n
\$2=outdir \n
\$3=candidate orf genepred \n
\$4=candidate orf fasta \n
\$5=Mart file transcript id and biotype\n"; exit 1; }
fi

set -uex

####### STEP 1: Generate candidate ORF list filtered for translating ORFs only ######################################################################
# 
# reduces the size of the candidate ORFs, easier to manipulate later in R
#

# Isolated list of ORF_IDs from translating_ORFs file: translating_ORF_ID.txt
cat $1 | cut -f 1 | grep -v orfID > $2/translating_ORF_ID.txt

# Selected candidate ORFs only translating: candidate_orfs_translating.tsv
# Using join command (need to sort the id file and table file first)
LC_ALL=C join -t$'\t' <(LC_ALL=C sort $2/translating_ORF_ID.txt) <(LC_ALL=C sort -t$'\t' -k1,1 $3) > $2/translating.candidateORF.genepred.txt

####### STEP 2: Generate candidate ORF fasta filtered for translating ORFs only ######################################################################

seqkit grep -f $2/translating_ORF_ID.txt $4 > $2/translating_candidateORF.fa

####### STEP 3: Translating candidate ORF GTF ########################################################################################################

genePredToGtf file $2/translating.candidateORF.genepred.txt $2/translating.candidateORF.gtf

####### STEP 4: Run R script to prepare GTF for VEP annotation and finalize GTF ######################################################################
### CHECK PATH OF RSCRIPT

Rscript /RiboNeo/riborf/variant_epitope/riborf_gtf_to_vep.R $2/translating.candidateORF.gtf $5 $2/translating.candidateORF.vep.gtf
# Add gene entry and sort chromosomes alphabetically
gffread -EF --keep-genes --sort-alpha $2/translating.candidateORF.vep.gtf -o- > $2/translating.candidateORF.vep.gff3
grep -v "#" $2/translating.candidateORF.vep.gff3 |\
  sort -k1,1 -k4,4n -k5,5n -t$'\t' |\
  bgzip -c > $2/translating.candidateORF.vep.gff3.gz && \
  tabix -p gff $2/translating.candidateORF.vep.gff3.gz


