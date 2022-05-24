#!/usr/bin/env python3
"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2022-02-23
Purpose: Runs Ribotricer on bam file from Ribo-Seq data
"""

import argparse
import subprocess
import os
import logging
import sys


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Process ribo-seq bam files with ribotricer tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('sample',
                        metavar='sample',
                        help='sample name')

    parser.add_argument('bam',
                        metavar='bam',
                        help='star output bam file')

    parser.add_argument('-f',
                        '--fasta',
                        metavar='\b',
                        help='Path to ref genome fasta file',
                        default='/Volumes/labagros/pierre/refs/gencode/GRCh38.primary_assembly.genome.fa')

    parser.add_argument('-g',
                        '--gtf',
                        metavar='\b',
                        help='Path to GTF file',
                        default='/Volumes/labagros/pierre/refs/combined.gencode.mitranscriptome.unannotated/combined.gencode.v38.mitranscriptome.v2.unannotated.sorted.gtf')

    parser.add_argument('-c',
                        '--candidate',
                        metavar='\b',
                        help='Path to candidate orf file',
                        default='/Volumes/labagros/pierre/proj/riboseq/ribotricer/v4_candidate_orfs.tsv')

    parser.add_argument('-i',
                        '--index_prefix',
                        metavar='\b',
                        help='Ribotricer index prefix',
                        default='ribotricer')

    parser.add_argument('-o',
                        '--output_prefix',
                        metavar='\b',
                        help='Ribotricer output prefix',
                        default='ribotricer')

    parser.add_argument('-T',
                        '--threads',
                        help='Number of threads to use',
                        metavar='\b', # metavar \b to not show any metavar except the short and long flag
                        type=int,
                        default=4)

    return parser.parse_args()


# --------------------------------------------------

def exec_command(cmd):
    logger = logging.getLogger()
    logger.info(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, executable='/bin/bash')
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            logger.error(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            logger.error(line.rstrip())
        sys.exit()

# --------------------------------------------------
def main():
    """Main function"""

    # assign arguments variable names
    args = get_args()
    sample = args.sample
    bam = args.bam
    fasta = args.fasta
    gtf = args.gtf
    candidate = args.candidate
    index_prefix = args.index_prefix
    output_prefix = args.output_prefix
    threads = args.threads

    # create output directory out and set as current wd
    if not os.path.exists(f'{os.getcwd()}/{sample}_out'):
        os.mkdir(f'{os.getcwd()}/{sample}_out') # create dir
    outdir = f'{os.getcwd()}/{sample}_out' # create outdir variable
    print(outdir)
    os.chdir(outdir) # set as new wd


    # Create log file (called ribotricer.log)
    logging.basicConfig(format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                        logging.FileHandler("ribotricer.log"),
                        logging.StreamHandler() # these 2 handlers allow 1) to have the log file created and 2) to stream to the terminal
                        ])

    logger = logging.getLogger() # creates riboseq logger to add entries to the log


    # Ribotricer prepare-orfs, if no candidate orf file is provided
    if not os.path.exists(candidate):
        logger.info("Start `ribotricer prepare-orfs` step")
        cmd = f'ribotricer prepare-orfs \
        --gtf {gtf} \
        --fasta {fasta} \
        --prefix {index_prefix} \
        --min_orf_length 24 \
        --start_codons ATG,TTG,CTG,GTG \
        --longest'
        exec_command(cmd)
        candidate <- f'{index_prefix}_candidate_orfs.tsv'

    # Ribotricer: detecting translating ORFs
    logger.info("Starting Ribotricer `detect-orfs` step")
    cmd = f'ribotricer detect-orfs \
             --bam {bam} \
             --ribotricer_index {candidate} \
             --prefix {output_prefix}'
    exec_command(cmd)

    print("DONE")


    # THE END
    logger.info(f"****** Ribotricer analysis completed!  ******")


# --------------------------------------------------
if __name__ == '__main__':
    main()