#!/usr/bin/env python3

"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2023-01-11
Purpose: Create AminoAcid and Nucleotide (CDS) dictionnaries
"""

import argparse
from Bio import SeqIO
import json

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Create AminoAcid and Nucleotide (CDS) dictionnaries from Fastas',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('outdir',
                        metavar='outdir',
                        help='Path to output directory')

    return parser.parse_args()

# --------------------------------------------------

# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    outdir = args.outdir

    # AA dict
    aa_dict = dict()
    with open(f"{outdir}/translating_candidateORF_AA.fa") as input_file:
        my_list = list(SeqIO.parse(input_file, "fasta"))

        for record in my_list:
            new_record = record.id.replace(":", "|")
            aa_dict[new_record] = str(record.seq)

        json.dump(aa_dict, open(f"{outdir}/translating_candidateORF_AA.dict",'w'))

    # CDS dict
    cds_dict = dict()
    with open(f"{outdir}/translating_candidateORF.fa") as input_file:
        my_list = list(SeqIO.parse(input_file, "fasta"))

        for record in my_list:
            new_record = record.id.replace(":", "|")
            cds_dict[new_record] = str(record.seq)

        json.dump(cds_dict, open(f"{outdir}/translating_candidateORF.dict",'w'))

# --------------------------------------------------
if __name__ == '__main__':
    main()