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

RUN pipeline from the RibORF directory
"""

import argparse
import os
import shutil
import logging
import sys
from pipeline.common import *

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='RibORF analysis of translating ORFs',
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
                        help='Path to genePred annotation file')

    parser.add_argument('-t',
                        '--transcriptGenePred',
                        metavar='\b',
                        help='Path to cDNA genePred annotation file',
                        default='/ref/gencode.v39.protein_coding.genePred')
    
    parser.add_argument('-o',
                        '--output',
                        metavar='\b',
                        help='Path to output dir',
                        default='~/Riborf_out')

    parser.add_argument('-d',
                        '--readlength',
                        metavar='\b',
                        help='read lengths to consider for readDist step',
                        default='32,33,34,35,36')
    
    parser.add_argument('-f',
                        '--offset',
                        metavar='\b',
                        help='optionnal: provide to create custom offset.correction.parameters.txt',
                        default='13,13,13,13,13')
    
    parser.add_argument('-c',
                        '--candidates',
                        metavar='\b',
                        help='path to candidateORF.genepred.txt')
    
    parser.add_argument('-r',
                        '--gtfToGenePred',
                        help='optionnal: runs gtfToGenePred on genome annotation file',
                        action='store_true')
    
    parser.add_argument('-w',
                        '--gtf',
                        help='optionnal: gtf annotation file if gtfToGenePred step required',
                        action='store_true')

    parser.add_argument('-a',
                        '--annotate',
                        help='optionnal: runs ORF annotate step',
                        action='store_true')

    return parser.parse_args()

# --------------------------------------------------
def main():
    """Main function"""

    # assign arguments variable names
    args = get_args()
    sample = args.sample
    bam = args.bam
    genome = args.genome
    gtf = args.gtf
    genePred = args.genePred
    transcript = args.transcriptGenePred
    readlength = args.readlength
    orf_annotate = args.annotate
    candidates = args.candidates
    output = args.output
    offset = args.offset
    gtfToGenePred = args.gtfToGenePred


    # create output directory
    outdir = f'{output}/{sample}_out' # create outdir variable
    if not os.path.exists(outdir):
        os.mkdir(outdir) # create dir

    # Create log file (called riborf.log)
    logging.basicConfig(format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                        logging.FileHandler("riborf.log"),
                        logging.StreamHandler() # these 2 handlers allow 1) to have the log file created and 2) to stream to the terminal
                        ])

    logger = logging.getLogger() # creates riboseq logger to add entries to the log

    # ******************************************************************************************

    # Log start time
    logger.info("Start RibORF pipeline")

    if gtfToGenePred:
        genePred = re.sub("\.gtf.*", "", gtf) + ".genePred"
        cmd = f'gtfToGenePred {gtf} {genePred}'

    if orf_annotate:
    # ORFannotate.pl
        logger.info("****** Generate Candidate ORFs ******")
        cmd = f'perl ORFannotate.pl -g {genome} -t {genePred} -o {outdir}'
        exec_command(cmd)

    # Convert STAR Bam file to Sam file, input for readDist script
    logger.info("****** Step 1 = Convert STAR Bam file to Sam file, input for readDist script ******")
    cmd = f"samtools view -h -o {outdir}/{sample}.dedup.sam {bam}"
    exec_command(cmd)

    # ReadDist.pl on Star output
    logger.info("****** Step 2 = readDist.pl using protein coding genePred ******")
    cmd = f"perl readDist.pl -f {outdir}/{sample}.dedup.sam -g {transcript} -o {outdir} -d {readlength}"
    exec_command(cmd)

    # Create offset.corretion.parameters.txt if different from default or copy it from pipeline dir if default
    if offset != "13,13,13,13,13":
        offsetCorrectionFile = open(f"{outdir}/offset.corretion.parameters.txt", "w")
        readlengths=readlength.split(",")
        offsets=offset.split(",")
        offsetCorrection=""
        for rl,off in zip(readlengths,offsets):
            offsetCorrection = offsetCorrection + rl + "\t" + off + "\n"
        offsetCorrectionFile.write(offsetCorrection)
        offsetCorrectionFile.close()
    else:
        shutil.copyfile("/RiboNeo/offset.corretion.parameters.txt", f"{outdir}/offset.corretion.parameters.txt")

    # offsetCorrect.pl: correct read locations based on offset distances between 5’ ends and ribosomal A-sites
    logger.info("****** Step 3 = offsetCorrect.pl ******")
    cmd = f"perl offsetCorrect.pl -r {outdir}/{sample}.dedup.sam -p {outdir}/offset.corretion.parameters.txt -o {outdir}/{sample}.offset_corrected.sam"
    exec_command(cmd)

    # readDist.pl on offset-corrected reads
    logger.info("****** Step 4 = readDist.pl on offset-correctd reads ******")
    cmd = f"perl readDist.pl -f {outdir}/{sample}.offset_corrected.sam -g {transcript} -o {outdir} -d 1"
    exec_command(cmd)

    # ribORF.pl to identify translated ORFs
    logger.info("****** Step 5 = ribORF.pl to identify translated ORFs ******")
    if orf_annotate:
        cmd = f"perl ribORF.pl -f {outdir}/{sample}.offset_corrected.sam -c {outdir}/candidateORF.genepred.txt -o {outdir}"
    else:
        cmd = f"perl ribORF.pl -f {outdir}/{sample}.offset_corrected.sam -c {candidates} -o {outdir}"
    exec_command(cmd)

    logger.info("****** ribORF.pl done ******")

    # Create offset-corrected bam file from sam file
    logger.info("****** Step 6 = Create offset-corrected bam file from sam file and index ******")
    cmd = f"samtools sort {outdir}/{sample}.offset_corrected.sam | samtools view -Sb -o {outdir}/{sample}.offset_corrected.bam \
    && samtools index {outdir}/{sample}.offset_corrected.bam"
    exec_command(cmd)

    logger.info("****** THE END ******")


# --------------------------------------------------
if __name__ == '__main__':
    main()
