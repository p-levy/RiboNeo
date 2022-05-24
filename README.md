# Ribo-seq Pipeline

TO BE CONTINUED...

## riboseq_to_bam.py

This pipeline processes fastq files from **RiboLACE seq** (Immagina Biotech) into BAM files aligned to the human genome. 
- Step 1: `fastqc` (qc of raw fastq files)
- Step 2: `cutadapt` (read trimming 1)
- Step 3: `umi_tools` (extract UMIs)
- Step 4: `cutadapt` (read trimming 2)
- Step 5: `bowtie2` (exclusion of rRNA reads)
- Step 6: `bowtie2` (exclusion of tRNA reads)
- Step 7: `star` (alignment to human genome (hg38 by default))
- Step 8: `samtools index` (required for next step umi_tools dedup)
- Step 9: `umi_tools dedup` (remove duplicate reads using UMIs)
- Step 10: `samtools index` (index final bam file)

Main outputs:
sample.dedup.bam
sample.dedup.bam.bai

## ribotricer
in progress...

## ribotricer_post.sh

This script processes `ribotricer` output table of *translating ORFs* and retrieves counts of **non-annotated ORF-specific intervals**.

**Usage**:
`ribotricer_post.sh /path/to/orf_file.txt prefix pos_wig neg_wig`

