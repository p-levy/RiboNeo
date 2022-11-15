#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(parallel))

# Usage
if (length(args)==0 | args[1]=="-h" | args[1]=="--help") {
  cat("
How to use:
arg1: path to non_canonical_specific.bed
arg2: path to pos wig file
arg3: path to neg wig file
arg4: path to riborf repre.valid.pred.pvalue.parameters.txt
arg5: output prefix
arg6: number of cores to use for parallel process

      ")
  quit()
}

args = commandArgs(trailingOnly=TRUE)
novel_bed <- args[1]
pos_wig <- args[2]
neg_wig <- args[3]
riborf <- args[4]
output_prefix <- args[5]
threads <- strtoi(args[6])
max_threads <- detectCores()

if (threads >= max_threads){
	stop(sprintf("Number of threads (%s) requested larger than available threads (%s). Stopping execution", threads, max_threads))
}

all_orfs <- fread(riborf)
non_annotated_orfs <- all_orfs %>% filter(str_detect(orfID, "canonical", negate = T))
novel <- fread(novel_bed)
colnames(novel) <- c("chrom", "start", "end", "ORF", "score", "strand")
novel$coordinate <- paste(novel$start, novel$end, sep="-")

wig_neg <- import.wig(neg_wig)
wig_pos <- import.wig(pos_wig)

positive_orfs <- novel[ novel$strand != "-",]
negative_orfs <- novel[ novel$strand != "+",]

cl <- makeCluster(threads)
invisible(clusterEvalQ(cl, c(library(GenomicRanges))))
clusterExport(cl, c("score", "wig_neg", "wig_pos"), envir=environment())

positive_orfs$counts <- parallel::parApply(cl, positive_orfs, 1, function(x) sum(score(wig_pos[seqnames(wig_pos) == x["chrom"] & x["start"] <= start(wig_pos) & end(wig_pos) <= x["end"]])))
negative_orfs$counts <- parallel::parApply(cl, negative_orfs, 1, function(x) sum(score(wig_neg[seqnames(wig_neg) == x["chrom"] & x["start"] <= start(wig_neg) & end(wig_neg) <= x["end"]])))

positive_orfs <- positive_orfs %>% group_by(ORF) %>% mutate(reads = sum(counts)) %>% select(-c(counts, start, end, score)) %>% mutate(novel_coordinate = paste(coordinate, collapse=", ")) %>% select(-c(coordinate)) %>% dplyr::slice(1) %>% ungroup()
negative_orfs <- negative_orfs %>% group_by(ORF) %>% mutate(reads = sum(counts)) %>% select(-c(counts, start, end, score)) %>% mutate(novel_coordinate = paste(coordinate, collapse=", ")) %>% select(-c(coordinate)) %>% dplyr::slice(1) %>% ungroup()

all_new_orfs <- rbind(positive_orfs, negative_orfs)

final_orf_data <- merge(non_annotated_orfs, all_new_orfs, by.x = c("orfID", "chrom", "strand"), by.y = c("ORF", "chrom", "strand"))

write.table(final_orf_data, paste(output_prefix, "translating_ORFs.txt", sep = "_"), sep = "\t", row.names=F, col.names=T, quote=F)
