# riboseq_pipeline

This pipeline processes fastq files from **ribosome profiling** runs into BAM files aligned to the human genome. 
- Step 1: `fastp` (read trimming and qc)
- Step 2: `star` alignment to human genome (hg38 by default)
