# Ribo-seq Pipeline


## Pre-processing and mapping to genome

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
- Step 11: `ribowaltz` (creates different qc plots of riboseq data)

Main outputs:
<br>
`sample.dedup.bam`, `sample.dedup.bam.bai` and `ribowaltz` pdf files.

Run: `./preprocess.py [-h] [-i] [-j] [-l] [-u] [-r] [-s] [-g] [-a] [-T] sample fastq`

```
positional arguments:
  sample             sample name
  fastq              path to fastq

optional arguments:
  -h, --help         show this help message and exit
  -i, --rrnai    Path to rRNA bowtie2 index prefix (default: refs/rRNA)
  -j, --trnai    Path to tRNA bowtie2 index prefix (default: refs/tRNA)
  -l, --lowlen   lower length filter RiboWaltz (default: 28)
  -u, --uplen    upper length filter RiboWaltz (default: 36)
  -r, --rds      Path to RiboWalzt annotation RDS (default: None)
  -s, --stari    Path to STAR index (default: None)
  -g, --gtf      Path to GTF file (default: None)
  -a, --annotate     optionnal: runs riboWalz create_annotation(gtf) (default:
                     False)
  -T, --threads  Number of threads to use (default: 4)
```

## RibORF

Run: `./riborf.py [-h] [-g] [-p] [-t] [-d] [-T] [-a] sample bam`

```
positional arguments:
  sample                sample name
  bam                   bam output from STAR alignment

optional arguments:
  -h, --help            show this help message and exit
  -g, --genome      Path to genome ref fasta file (default:
                        /ref/GRCh38.primary_assembly.genome.fa)
  -p, --genePred    Path to genePred annotation file (default: /ref/combin
                        ed.gencode.v39.mitranscriptome.unique_id.v2.unannotate
                        d.sorted.genePred)
  -t, --transcriptGenePred 
                        Path to cDNA genePred annotation file (default:
                        /ref/gencode.v39.protein_coding.genePred)
  -d, --readlength  read lengths to consider for readDist step (default:
                        28,29,30,31,32,33,34,35,36)
  -T, --threads     Number of threads to use (default: 30)
  -a, --annotate        optionnal: runs ORF annotate step (default: False)
```



