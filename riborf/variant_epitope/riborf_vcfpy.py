#!/usr/bin/env python3

"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2023-01-11
Purpose: Variant selection from VEP-annotated VCF file using GTF file containing only translating ORFs from Riborf analysis
"""

import argparse
import os
import vcfpy
import numpy as np
import re
from collections import namedtuple

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Variant selection from VEP-annotated VCF file using GTF file \
        containing only translating ORFs from Riborf analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('vcf_path',
                        metavar='vcf_path',
                        help='Path to VEP-annotated VCF file')

    parser.add_argument('-o',
                        '--out',
                        help='Output dir name',
                        metavar='\b',
                        type=str,
                        default='out')

    return parser.parse_args()

# --------------------------------------------------

# -------------Variant selection filters----------------

t2n_ratio = 5
tumor_var_freq = 7
tumor_var_depth = 4
normal_var_freq = 100
normal_coverage = 10
tumor_coverage = 10
num_callers = 1
num_callers_indel = 1

# --------------------------------------------------

# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    vcf_path = args.vcf_path
    out = args.out

    # create output directory out and set as current wd
    if not os.path.exists(f'{out}'):
        os.mkdir(f'{out}') # create dir
    outdir = f'{out}' # create outdir variable
    # os.chdir(outdir) # set as new wd (uncomment if desired)

    vcf_in = vcfpy.Reader.from_path(vcf_path)

    for x in vcf_in.header.lines:
        if "CSQ" in str(x):
            csq = str(x)

    try:
        csq = re.sub(r'^.*?Format: .*?Format: ', '', csq)
    except UnboundLocalError:
        print("Are you sure this VCF file was annotated with VEP?")
        exit()

    csq = re.sub('\|', ' ', csq)
    csq = re.sub(' translating.candidateORF.vep.gff3.gz', ' GFF', csq)
    csq = re.sub("\'\}\)", '', csq)

    Record_INFO = namedtuple('INFO', csq)

    class Variant:
        def __init__(self):
            self.chrom = None
            self.start = None
            self.ref = None
            self.alt = None
            self.epitopes = None
            self.callers = None
            self.num_callers = None
            self.status = None
            self.gene = None
            self.feature = None
            self.consequence = None
            self.cdna_pos = None
            self.cds_pos = None
            self.protein_pos = None
            self.aa = None
            self.codon = None

        def __str__(self):
            return '{}:{}\t{}>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                self.chrom, self.start, self.ref, self.alt, self.gene, self.feature, self.consequence, self.status, self.cdna_pos, self.cds_pos, self.protein_pos, self.aa, self.codon)

    variants = list()
    for record in vcf_in:
            for info in record.INFO['CSQ']:
                record_INFO = Record_INFO(*info.split('|'))
                funcensGene = record_INFO.Consequence
                has_func_ens = 'missense' in funcensGene or 'frame' in funcensGene or 'stop_lost' in funcensGene
                if has_func_ens:
                    called = {x.sample: x.data for x in record.calls if x.called}
                    filtered = dict()
                    pass_snp = 0
                    pass_indel = 0
                    try:
                        if 'NORMAL.mutect' in called and 'TUMOR.mutect' in called and 'PASS' in record.FILTER:
                            normal_DP = int(called['NORMAL.mutect']['DP'])
                            token = called['NORMAL.mutect']['AD']
                            normal_AD = int(token[1]) if type(token) is list else int(token)
                            token = called['NORMAL.mutect']['AF']
                            value = token[0] if type(token) is list else token
                            normal_VAF = np.around(float(value) * 100,3) if normal_DP > 0.0 else 0.0
                            tumor_DP = int(called['TUMOR.mutect']['DP'])
                            token = called['TUMOR.mutect']['AD']
                            tumor_AD = int(token[1]) if type(token) is list else int(token)
                            token = called['TUMOR.mutect']['AF']
                            value = token[0] if type(token) is list else token
                            tumor_VAF = np.around(float(value) * 100, 3)
                            tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                            if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                    and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                    and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                                pass_snp += 1
                            filtered['mutect'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                            normal_AD,
                                                                            normal_VAF,
                                                                            tumor_DP,
                                                                            tumor_AD,
                                                                            tumor_VAF)
                        if 'NORMAL.somaticsniper' in called and 'TUMOR.somaticsniper' in called:
                            normal_DP = int(called['NORMAL.somaticsniper']['DP'])
                            normal_AD = sum(called['NORMAL.somaticsniper']['DP4'][2:])
                            normal_VAF = np.around((normal_AD / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                            tumor_DP = int(called['TUMOR.somaticsniper']['DP'])
                            tumor_AD = sum(called['TUMOR.somaticsniper']['DP4'][2:])
                            tumor_VAF = np.around((tumor_AD / float(tumor_DP)) * 100, 3)
                            tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                            is_somatic = int(called['TUMOR.somaticsniper']['SS']) == 2
                            if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                    and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                    and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio and is_somatic:
                                pass_snp += 1
                            if is_somatic:
                                filtered['somaticsniper'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                                       normal_AD,
                                                                                       normal_VAF,
                                                                                       tumor_DP,
                                                                                       tumor_AD,
                                                                                       tumor_VAF)
                        if ('NORMAL.varscan' in called and 'TUMOR.varscan' in called) \
                                or ('NORMAL.varscan_indel' in called and 'TUMOR.varscan_indel' in called) \
                                and 'PASS' in record.FILTER and 'SOMATIC' in record.INFO:
                            label_index = 'varscan' if 'NORMAL.varscan' in called else 'varscan_indel'
                            normal_DP = int(called['NORMAL.{}'.format(label_index)]['DP'])
                            normal_AD = sum(called['NORMAL.{}'.format(label_index)]['DP4'][2:])
                            token = called['NORMAL.{}'.format(label_index)]['FREQ']
                            value = token[0] if type(token) is list else token
                            normal_VAF = float(value.replace('%', ''))
                            tumor_DP = int(called['TUMOR.{}'.format(label_index)]['DP'])
                            tumor_AD = sum(called['TUMOR.{}'.format(label_index)]['DP4'][2:])
                            token = called['TUMOR.{}'.format(label_index)]['FREQ']
                            value = token[0] if type(token) is list else token
                            tumor_VAF = float(value.replace('%', ''))
                            tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                            if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                    and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                    and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                                if 'indel' in label_index:
                                    pass_indel += 1
                                else:
                                    pass_snp += 1
                            filtered[label_index] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                               normal_AD,
                                                                               normal_VAF,
                                                                               tumor_DP,
                                                                               tumor_AD,
                                                                               tumor_VAF)
                        if 'NORMAL.strelka' in called and 'TUMOR.strelka' in called and 'PASS' in record.FILTER:
                            ref_index = record.REF + 'U'
                            alt_index = str(record.ALT[0].serialize()) + 'U'
                            # normal_DP = int(called['NORMAL.strelka']['DP'])
                            token = called['NORMAL.strelka'][ref_index]
                            normal_AD1 = int(token[0]) if type(token) is list else int(token)
                            token = called['NORMAL.strelka'][alt_index]
                            normal_AD2 = int(token[0]) if type(token) is list else int(token)
                            normal_DP = normal_AD1 + normal_AD2
                            normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                            # tumor_DP = int(called['TUMOR.strelka']['DP'])
                            token = called['TUMOR.strelka'][ref_index]
                            tumor_AD1 = int(token[0]) if type(token) is list else int(token)
                            token = called['TUMOR.strelka'][alt_index]
                            tumor_AD2 = int(token[0]) if type(token) is list else int(token)
                            tumor_DP = tumor_AD1 + tumor_AD2
                            tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                            tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                            if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                    and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                                    and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                                pass_snp += 1
                            filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                             normal_AD2,
                                                                             normal_VAF,
                                                                             tumor_DP,
                                                                             tumor_AD2,
                                                                             tumor_VAF)
                        if 'NORMAL.strelka_indel' in called and 'TUMOR.strelka_indel' in called and 'PASS' in record.FILTER:
                            # normal_DP = int(called['NORMAL.strelka_indel']['DP'])
                            token = called['NORMAL.strelka_indel']['TAR']
                            normal_AD1 = int(token[0]) if type(token) is list else int(token)
                            token = called['NORMAL.strelka_indel']['TIR']
                            normal_AD2 = int(token[0]) if type(token) is list else int(token)
                            normal_DP = normal_AD1 + normal_AD2
                            normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                            # tumor_DP = int(called['TUMOR.strelka_indel']['DP'])
                            token = called['TUMOR.strelka_indel']['TAR']
                            tumor_AD1 = int(token[0]) if type(token) is list else int(token)
                            token = called['TUMOR.strelka_indel']['TIR']
                            tumor_AD2 = int(token[0]) if type(token) is list else int(token)
                            tumor_DP = tumor_AD1 + tumor_AD2
                            tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                            tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                            if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                    and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                                    and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                                pass_indel += 1
                            filtered['strelka_indel'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                                   normal_AD2,
                                                                                   normal_VAF,
                                                                                   tumor_DP,
                                                                                   tumor_AD2,
                                                                                   tumor_VAF)
                    except KeyError:
                        continue
                    
                    variant = Variant()
                    variant.chrom = record.CHROM
                    variant.start = record.POS
                    variant.ref = record.REF
                    variant.alt = record.ALT[0].serialize()
                    variant.callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
                    variant.num_callers = len(filtered)
                    variant.status = pass_snp >= num_callers or pass_indel >= num_callers_indel
                    variant.gene = record_INFO.Gene
                    variant.feature = record_INFO.Feature
                    variant.consequence = record_INFO.Consequence
                    variant.cdna_pos = record_INFO.cDNA_position
                    variant.cds_pos = record_INFO.CDS_position
                    variant.protein_pos = record_INFO.Protein_position
                    variant.aa = record_INFO.Amino_acids
                    variant.codon = record_INFO.Codons
                    variants.append(variant)

    output = open(f'{outdir}/output_variants.tsv', "w")
    for line in variants:
        output.write(f"{line}\n")

    output.close()
    vcf_in.close()

# --------------------------------------------------
if __name__ == '__main__':
    main()

