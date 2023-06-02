#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)==0 || args[1]=="-h" || args[1]=="--help") {
  cat("
How to use:
arg1: translating.candidateORF.gtf
arg2: MART file: trancript_id_gene_id_trancript_type_gene_type.txt
arg3: output GTF path

      ")
  quit()
}

gtf_path <- args[1]
mart_path <- args[2]
output_gtf <- args[3]

gtf <- rtracklayer::import(gtf_path)

gtf2 <- as.data.frame(gtf)

# Split orfID column into multiple fields (to get transcript id)
gtf2 <- gtf2 %>% mutate(gene_ID = gene_id) %>%  separate(., gene_ID,
                                      c("transcript", 
                                        "chrom", 
                                        "strand2", 
                                        "RankNumber", 
                                        "transcriptLength", 
                                        "startCodonPosition", 
                                        "stopCodonPosition", 
                                        "candidateORFType", 
                                        "startCodonType"), sep = '[|:]')

gtf2$transcript <- gsub("\\..*", "", gtf2$transcript)

#Import Mart trancript type and gene id table
transcript_type <- fread(mart_path) %>% dplyr::rename("transcript" = "Transcript stable ID")

# Merge with GTF
gtf2 <- gtf2 %>% left_join(transcript_type, by = "transcript")
names(gtf2)[24] <- "biotype"
#names(gtf2)[25] <- "gene_biotype"

# change gene_id, transcript_id and exon_id columns
gtf2$gene_id <- gtf2$`Gene stable ID`
gtf2$exon_id <- gsub("\\..*TG\\.", "\\.", gtf2$exon_id)

gtf2$source <- "riborf"

gtf3 <- gtf2 %>% select(1:13, 24)

rtracklayer::export(gtf3, output_gtf)