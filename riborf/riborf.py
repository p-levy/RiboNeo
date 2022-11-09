#!/usr/bin/env python3
"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2022-02-23
Purpose: Riborf

requirements:
Bam file from the ribolace_process.py script (STAR-aligned)
Samtools
gtfToGenePred to make Transcript annotation file in genePred format
Genome assembly file in Fasta format
Linux high performance computing cluster
Perl program installation
R program installation in the PATH
RibORF package: https://github.com/zhejilab/RibORF/
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
        description='Process ribo-seq fastq files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('sample',
                        metavar='sample',
                        help='sample name')

    parser.add_argument('bam',
                        metavar='bam',
                        help='bam output from STAR alignment')

    parser.add_argument('-g',
                        '--genome',
                        metavar='\b',
                        help='Path to genome ref fasta file',
                        default='/ref/GRCh38.primary_assembly.genome.fa')

    parser.add_argument('-p',
                        '--genePred',
                        metavar='\b',
                        help='Path to genePred annotation file',
                        default='/ref/combined.gencode.v39.mitranscriptome.unique_id.v2.unannotated.sorted.genePred')

    parser.add_argument('-t',
                        '--transcriptGenePred',
                        metavar='\b',
                        help='Path to cDNA genePred annotation file',
                        default='/ref/gencode.v39.protein_coding.genePred')

    parser.add_argument('-d',
                        '--readlength',
                        metavar='\b',
                        help='read lengths to consider for readDist step',
                        default='28,29,30,31,32,33,34,35,36')


    parser.add_argument('-T',
                        '--threads',
                        help='Number of threads to use',
                        metavar='\b', # metavar \b to not show any metavar except the short and long flag
                        type=int,
                        default=30)

    parser.add_argument('-a',
                        '--annotate',
                        help='optionnal: runs ORF annotate step',
                        action='store_true')

    return parser.parse_args()


# --------------------------------------------------
# def exec_command(cmd):
#     """Execute commands in shell and logs it"""
#     logger = logging.getLogger()
#     logger.info(cmd)
#     subprocess.run(cmd, shell=True, executable='/bin/bash')

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
    genome = args.genome
    genePred = args.genePred
    transcript = args.transcriptGenePred
    readlength = args.readlength
    threads = args.threads
    orf_annotate = args.annotate

    # create output directory out and set as current wd
    if not os.path.exists(f'{os.getcwd()}/{sample}_out'):
        os.mkdir(f'{os.getcwd()}/{sample}_out') # create dir
    outdir = f'{os.getcwd()}/{sample}_out' # create outdir variable
    print(outdir)
    os.chdir(outdir) # set as new wd


    # Create log file (called riboseq.log)
    logging.basicConfig(format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                        logging.FileHandler("riboseq.log"),
                        logging.StreamHandler() # these 2 handlers allow 1) to have the log file created and 2) to stream to the terminal
                        ])

    logger = logging.getLogger() # creates riboseq logger to add entries to the log

    # ******************************************************************************************

    # Log start time
    logger.info("(Re-)Start RibORF pipeline")

    if orf_annotate:
    # ORFannotate.pl
        logger.info("****** Step 1 = Generate Candidate ORFs ******")
        cmd = f'perl /RibORF.1.0/ORFannotate.pl -g {genome} -t {genePred} -o {outdir}'
        exec_command(cmd)

    # Convert STAR Bam file to Sam file, input for readDist script
    logger.info("****** Step 2 = Convert STAR Bam file to Sam file, input for readDist script ******")
    cmd = f"samtools view -h -o {outdir}/{sample}.dedup.sam {bam}"
    exec_command(cmd)

    # ReadDist.pl on Star output
    logger.info("****** Step 3 = readDist.pl using protein coding genePred ******")
    cmd = f"perl /RibORF.1.0/readDist.pl -f {outdir}/{sample}.dedup.sam -g {transcript} -o {outdir} -d {readlength}"
    exec_command(cmd)


    #####################################################################
    ## STOP HERE TO MAKE MANUALLY THE offset.correction.parameters.txt ##
    #####################################################################

    # # offsetCorrect.pl: correct read locations based on offset distances between 5’ ends and ribosomal A-sites
    # logger.info("****** Step 4 = offsetCorrect.pl ******")
    # cmd = f"perl /RibORF.1.0/offsetCorrect.pl -r {outdir}/{sample}.dedup.sam -p {outdir}/offset.corretion.parameters.txt -o {outdir}/{sample}.offset_corrected.sam"
    # exec_command(cmd)

    # # readDist.pl on offset-corrected reads
    # logger.info("****** Step 5 = readDist.pl on offset-correctd reads ******")
    # cmd = f"perl /RibORF.1.0/readDist.pl -f {outdir}/{sample}.offset_corrected.sam -g {transcript} -o {outdir} -d 1"
    # exec_command(cmd)

    # # ribORF.pl to identify translated ORFs
    # logger.info("****** Step 6 = ribORF.pl to identify translated ORFs ******")
    # cmd = f"perl /RibORF.1.0/ribORF.pl -f {outdir}/{sample}.offset_corrected.sam -c /riboseq/riborf/annotate/3 -o {outdir}"
    # exec_command(cmd)

    # logger.info("****** ribORF.pl done ******")


    # logger.info("****** THE END ******")


# --------------------------------------------------
if __name__ == '__main__':
    main()
