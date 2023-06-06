#!/usr/bin/env python3
"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2023-01-25
Purpose: Create epitope for each variant
"""

import argparse
import re
import os
import json
import pandas as pd
from Bio.Seq import translate
from Bio.Seq import reverse_complement


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Create epitope for each variant',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('variants',
                        metavar='variants',
                        help='path to non_canonical_variants.tsv (output of nc_variants.R')

    parser.add_argument('-c',
                        '--cds',
                        help='CDS dictionary',
                        metavar='\b',
                        type=str,
                        default='/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF.dict')

    parser.add_argument('-d',
                        '--cdna',
                        help='cDNA dictionary',
                        metavar='\b',
                        type=str,
                        default='/RiboNeo/refs/Homo_sapiens.GRCh38.cdna.all.json')

    parser.add_argument('-a',
                        '--aa',
                        help='AA dictionary',
                        metavar='\b',
                        type=str,
                        default='/Volumes/labagros/pierre/proj/riboseq/riborf/360RIO136_gencode_out/translating_candidateORF_AA.dict')

    parser.add_argument('-o',
                        '--out',
                        help='Output dir name',
                        metavar='\b',
                        type=str,
                        default='out')

    return parser.parse_args()


# --------------------------------------------------
def main():

    args = get_args()
    variants = args.variants
    cds = args.cds
    aa = args.aa
    cdna = args.cdna
    outdir = args.out

    # create output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Determine length of epitope (start+end)
    mer_start = 13
    mer_end = 12

    # Open CDS, AA and cDNA dictionaries
    cds_dict = json.load(open(cds))
    aa_dict = json.load(open(aa))
    cdna_dict = json.load(open(cdna))

    variant_table = pd.read_table(variants)
    # create subdictionnary for cdnas with frameshift
    ## list of transripts with fs: 
    fs_transcripts = []
    for i in variant_table.iterrows():
        if "frameshift" in i[1]["type"]:
            fs_transcripts.append(re.sub("\.\d+?","",i[1]["transcript"])) # remove version after ENST name
    ## subdictionnary for cdnas with frameshift
    cdna_fs_dict = {k: cdna_dict[k] for k in fs_transcripts}

    # Create columns to receive wt and mut epitopes
    variant_table = variant_table.assign(wt_mer = '', mut_mer = '', flags = '')

    for index, row in variant_table.iterrows(): # Iterrate by row to get compute epitopes (wt and mut)     
        ref = row["ref_alt"].split(">")[0]
        cdna_pos = int(re.sub("-\d+", "", str(row["cdna_pos"]))) #removes part after "-" for inframe indels
        aa_pos = int(re.sub("-\d+", "", str(row["protein_pos"]))) #removes part after "-" for inframe indels
        start_codon =int(row["startCodonPosition"])
        orfID = row["orfID"]
        fs = len(ref)-1 # -1 because the ref shows also the nucleotide that is still present in the alt allele
        start = aa_pos - mer_start
        flags = ""
        if start < 0:
            flags += ' Start_of_sequence_is_shorter_than_12aa_from_mutation'
            start = 0
        
        # Check what is the type of variant and compute wt_mer and mut_mer accordingly
        if 'missense' in row["type"]:
            var_AA = re.sub("./","",row["aa_change"])
            protein_seq = aa_dict[orfID]
            end = aa_pos + mer_end if aa_pos + mer_end < len(protein_seq) else None
            if end == None:
                flags += ' Stop_codon_is_less_than_12aa_from_mutation'
            wt_mer = protein_seq[start:end]
            mut_mer = protein_seq[start:aa_pos - 1] + var_AA + protein_seq[aa_pos:end]
        
        elif 'frameshift_variant' in row["type"] or 'stop_lost' in row["type"]:
            ins_seq = ""
            mut_cds = ""
            var_nt = ""
            cdna_seq = cdna_fs_dict[re.sub("\.\d+?","",row["transcript"])]
            if len(ref) == 1 and ('stop_lost' not in row["type"]): #true if a fs insertion
                ins_seq = row["ref_alt"].split(">")[1][1:] if row["strand"] == "+" else reverse_complement(row["ref_alt"].split(">")[1][1:]) # [1:] because first nt is the ref
                mut_cds = cdna_seq[start_codon-1:cdna_pos] + ins_seq + cdna_seq[cdna_pos:] # ins_seq[1:] start at index 1 to not include ref in inserted seq
            elif len(ref) == 1 and ('stop_lost' in row["type"]): #true if a stop_lost
                var_nt = row["ref_alt"].split(">")[1] if row["strand"] == "+" else reverse_complement(row["ref_alt"].split(">")[1])
                mut_cds = cdna_seq[start_codon-1:cdna_pos-1] + var_nt + cdna_seq[cdna_pos:]
            elif len(ref) > 1: #true if a fs deletion        
                mut_cds = cdna_seq[start_codon-1:cdna_pos-1] + cdna_seq[cdna_pos + fs - 1:]
            wt_mer = aa_dict[orfID][start:aa_pos+mer_start]
            mut_fasta = str(translate(mut_cds.replace(' ', ''), to_stop=True))
            mut_mer = mut_fasta[start:]
            if len(mut_mer) < len(wt_mer):
                flags += ' Fs_resulted_in_a_stop_codon_before_12_aa_from_mutation'
            
        elif 'inframe' in row["type"]:
            n_AA = int()
            protein_seq = aa_dict[orfID]
            end = aa_pos + mer_end if aa_pos + mer_end < len(protein_seq) else None
            if end == None:
                flags += ' Stop_codon_is_less_than_12aa_from_mutation'
            wt_mer = protein_seq[start:end]
            if len(ref) > 1: #true if deletion
                n_AA = int((len(row["ref_alt"].split(">")[0])-1)/3) #-1 to remove ref and /3 to get number of AAs
                mut_mer = protein_seq[start:aa_pos - 1] + protein_seq[aa_pos+n_AA-1:end]
            elif len(ref) == 1: #true if insertion
                ins_AA = translate(row["ref_alt"].split(">")[1][1:]) if row["strand"] == "+" else translate(reverse_complement(row["ref_alt"].split(">")[1][1:]))
                mut_mer = protein_seq[start:aa_pos - 1] + ins_AA + protein_seq[aa_pos:end]
                
        variant_table._set_value(index,'wt_mer',wt_mer) # add wt epitope seq to wt_mer column 
        variant_table._set_value(index,'mut_mer',mut_mer) # add mut epitope seq to mut_mer column
        variant_table._set_value(index,'flags',flags)

    # Create columns to check each variant locus with previous locus and see if they are within 39-nt distance
    variant_table['chrom'] = variant_table['locus'].transform(lambda x: x.split(':')[0])
    variant_table['pos'] = variant_table['locus'].transform(lambda x: int(x.split(':')[1]))
    variant_table['mut_diff'] = variant_table.groupby('chrom')['pos'].diff()
    variant_table['mut_diff'] = variant_table.groupby('locus')['mut_diff'].transform(max)
    variant_table.loc[(variant_table['mut_diff']<39) & (variant_table['mut_diff']!=0), 'flags'] = variant_table['flags'].transform(lambda x: x+' Passing_variant_within_39nt')
    # remove extra columns 
    variant_table = variant_table.drop(columns=['chrom', 'pos', 'mut_diff'])

    variant_table.to_csv(f'{outdir}/epitopes.tsv', sep = '\t', index=False)
# --------------------------------------------------
if __name__ == '__main__':
    main()
