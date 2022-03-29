#!/usr/bin/env python3
"""
Author : Pierre Levy <levy.pierre@yahoo.fr>
Date   : 2022-02-23
Purpose: Process ribo-seq fastq files from Ribolace seq until STAR alignment
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

    parser.add_argument('fastq',
                        metavar='fastq',
                        help='raw fastq')

    parser.add_argument('-r',
                        '--rrna',
                        metavar='\b',
                        help='Path to rRNA ref fasta',
                        default='/ref/rrna_trna/Rfam_rRNA.fa')

    parser.add_argument('-t',
                        '--trna',
                        metavar='\b',
                        help='Path to tRNA ref fasta',
                        default='/ref/rrna_trna/Rfam_tRNA.fa')


    parser.add_argument('-i',
                        '--rrnai',
                        metavar='\b',
                        help='Path to rRNA bowtie2 index prefix',
                        default='/ref/rrna_trna/rRNA')

    parser.add_argument('-j',
                        '--trnai',
                        metavar='\b',
                        help='Path to tRNA bowtie2 index prefix',
                        default='/ref/rrna_trna/tRNA')


    parser.add_argument('-s',
                        '--stari',
                        metavar='\b',
                        help='Path to STAR index',
                        default='/ref/STARIndex')


    parser.add_argument('-g',
                        '--gtf',
                        metavar='\b',
                        help='Path to GTF file',
                        default='/ref/combined.gencode.v38.mitranscriptome.v2.unannotated.sorted.gtf')


    parser.add_argument('-T',
                        '--threads',
                        help='Number of threads to use',
                        metavar='\b', # metavar \b to not show any metavar except the short and long flag
                        type=int,
                        default=30)

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
    fastq = args.fastq
    threads = args.threads
    rrna = args.rrna
    trna = args.trna
    rrna_i = args.rrnai
    trna_i = args.trnai
    starindex = args.stari
    gtf = args.gtf

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


    # Log start time
    logger.info("Start Ribo-Seq Process")

    # qc for raw fastqs
    logger.info("****** Step 1 = Fastqc ******")
    if not os.path.exists(f'{outdir}/qc'):
        os.mkdir(f'{outdir}/qc') # create dir
    cmd = f'fastqc --nogroup --threads {threads} {fastq} -o {outdir}/qc' # define shell command
    exec_command(cmd) # execute shell command (see above for exec_command function def)

    # cutadapt (1)
    logger.info("****** Step 2 = Cutadapt (1) ******")
    cmd = f'cutadapt --minimum-length 20 -a TCTCCTTGCATAATCACCAACC --discard-untrimmed --cores={threads} \
    -o trim.fastq {fastq} >> riboseq.log'
    exec_command(cmd)

    # umi_tools (extract UMIs)
    logger.info("****** Step 3 = Umi_tools (extract UMIs) ******")
    cmd = f"umi_tools extract -I {outdir}/trim.fastq --bc-pattern='(?P<umi_1>.{{4}}).+(?P<umi_2>.{{4}})$' \
    --extract-method=regex -S {outdir}/extract.fq"
    exec_command(cmd)

    # cutadapt (2)
    logger.info("****** Step 4 = Cutadapt (2) ******")
    cmd = f"cutadapt --cut 1 -o {outdir}/trim2.fq {outdir}/extract.fq --cores={threads} >> riboseq.log"
    exec_command(cmd)

    # Bowtie2 exclusion of rRNA reads  ******* CHANGE PATH TO INDEX!!! MAKE VARIABLE FLAG *********
    logger.info("****** Step 5 = Bowtie2 exclusion of rRNA reads ******")
    cmd = f"bowtie2 -p {threads} -N 1 --no-1mm-upfront \
                 -q {outdir}/trim2.fq \
                 --un {outdir}/no_rRNA.fq \
                 -x {rrna_i} \
                 > /dev/null 2> >(tee -a {outdir}/riboseq.log >&2)" #sends stdout to null, stderr to log file and redirects it to stderr to show on terminal screen
    exec_command(cmd)

    # Bowtie2 exclusion of tRNA reads
    logger.info("****** Step 6 = Bowtie2 exclusion of tRNA reads ******")
    cmd = f"bowtie2 -p {threads} -N 1 --no-1mm-upfront \
            -q {outdir}/no_rRNA.fq \
            --un {outdir}/no_rRNA_no_tRNA.fq \
            -x {trna_i} \
            > /dev/null 2> >(tee -a {outdir}/riboseq.log >&2)" #sends stdout to null, stderr to log file and redirects it to stderr to show on terminal screen
    exec_command(cmd)

    # STAR mapping to genome
    logger.info("****** Step 7 = STAR mapping to genome ******")
    cmd = f"STAR --genomeDir {starindex} \
     --readFilesIn {outdir}/no_rRNA_no_tRNA.fq \
     --runThreadN {threads} \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --alignEndsType EndToEnd \
     -outFilterMultimapNmax 1 \
     --sjdbGTFfile {gtf} \
     --alignIntronMin 20 \
     --alignIntronMax 100000 \
     --outFilterMismatchNmax 1 \
     --outFilterType BySJout \
     --outFilterMismatchNoverLmax 0.04 \
     --twopassMode Basic"
    exec_command(cmd)
    
    # Samtools index
    logger.info("****** Step 8 = Samtools index of STAR output bam ******")
    cmd = f"samtools index {outdir}/Aligned.sortedByCoord.out.bam {outdir}/Aligned.sortedByCoord.out.bam.bai"
    exec_command(cmd)

    # umi_tools (remove duplicates)
    logger.info("****** Step 9 = Umi_tools (remove duplicates) ******")
    cmd = f"umi_tools dedup -I {outdir}/Aligned.sortedByCoord.out.bam -S {outdir}/{sample}.dedup.bam \
    --method=directional --log2stderr 2> >(tee -a {outdir}/riboseq.log >&2)"
    exec_command(cmd)

    # Samtools index 2
    logger.info("****** Step 10 = Samtools index final bam ******")
    cmd = f"samtools index {outdir}/{sample}.dedup.bam {outdir}/{sample}.dedup.bam.bai"
    exec_command(cmd)


    # THE END
    logger.info(f"****** Pipeline completed! ---> {sample}.dedup.bam file generated correctly ******")


# --------------------------------------------------
if __name__ == '__main__':
    main()
