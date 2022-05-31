#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(trailingOnly=TRUE)
candidate <- args[1]
longest <- args[2]
output_prefix <- args[3]

candidate_orfs <- fread(candidate)
longest_orfs <- fread(candidate)

annotated_orfs <- candidate_orfs[ candidate_orfs$ORF_type == "annotated",]
non_annotated_orfs <- longest_orfs[ longest_orfs$ORF_type != "annotated",]

annotated_orfs <- annotated_orfs %>% mutate(coordinate = strsplit(as.character(coordinate), ',')) %>% unnest(coordinate)
non_annotated_orfs <- non_annotated_orfs %>% mutate(coordinate = strsplit(as.character(coordinate), ',')) %>% unnest(coordinate)

annotated_orfs <- annotated_orfs %>% separate(., coordinate, c("start", "end"), sep = '-')
non_annotated_orfs <- non_annotated_orfs %>% separate(., coordinate, c("start", "end"), sep = '-')

annotated_orfs$Score <- rep(0, nrow(annotated_orfs))
non_annotated_orfs$Score <- rep(0, nrow(non_annotated_orfs))

annotated_orfs <- annotated_orfs[,c(4,8,9,1,11, 6)]
non_annotated_orfs <- non_annotated_orfs[,c(4,8,9,1,11, 6)]

write.table(annotated_orfs, paste(output_prefix, "annotated_orfs.bed", sep = "_"), sep = "\t", row.names=F, col.names=F, quote=F)
write.table(non_annotated_orfs, paste(output_prefix, "non_annotated_orfs.bed", sep = "_"), sep = "\t", row.names=F, col.names=F, quote=F)
