#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)==0 | args[1]=="-h" | args[1]=="--help") {
  cat("
How to use:
arg1: path to ribotricer 'translating_ORFs.tsv'
arg2: path to 'candidate_orfs.tsv'
arg3: genome fasta
arg4: outdir

      ")
  quit()
}


# Argument variables
pathToTranslating <- args[1]
pathToCandidates <- args[2]
genome <- args[3]
outdir <- args[4]
# 
# # Filter candidate ORF table for translating ORFs only (to reduce size)
# system(paste0("./ribotricer_translating_candidates.sh", " ",
#               pathToTranslating, " ", 
#               pathToCandidates, " ", 
#               outdir))
# 
# # Get AA-sequence of filtered translating candidate ORFs
# system(paste0("./ribotricer_orfs_seq.sh", " ",
#               outdir, "candidate_orfs_translating_final.tsv", " ",
#               genome, " ",
#               outdir, "candidate_orfs_translating_final_aa_seq.tsv"))

# load data tables
ribotricer <- fread(args[1])
candidates <- fread(paste0(outdir, "candidate_orfs_translating_final.tsv"))
aa_seq <- fread(paste0(outdir, "candidate_orfs_translating_final_aa_seq.tsv"))

# Merge translating orfs with coordinate column of candidate orf file
ribotricer <- ribotricer %>% left_join(candidates %>% select(ORF_ID, coordinate), by = "ORF_ID")

# add position columns
ribotricer <- ribotricer %>% rowwise() %>% 
              mutate(start = str_split(ORF_ID, "_")[[1]][2], end = str_split(ORF_ID, "_")[[1]][3]) # We have to use dplyr::rowwise() here, otherwise str_split only works with the first row occurence

# add aa_length column and rename length to nt_length
ribotricer <- ribotricer %>% mutate(aa_length = round(length/3)) %>% rename(nt_length = length)

# reorder
ribotricer <- ribotricer %>%  select(ORF_ID, gene_name, ORF_type, chrom, start, end, nt_length, aa_length, read_density, phase_score, valid_codons_ratio, transcript_type, everything())

# Keep only unique ORFs
ribotricer <- ribotricer %>% distinct(chrom, start, end, .keep_all = TRUE)

# join translating ORF DT and aa_seq
ribotricer_aa <- ribotricer %>% left_join(aa_seq, by = "ORF_ID")

# Filter longest unique ORFs for each gene
df_uniq_longest_orf <- data.frame(matrix(nrow=0,ncol=24))
colnames(df_uniq_longest_orf) <- colnames(ribotricer_aa)

for (gene in unique(ribotricer_aa$gene_name)) {
  index <- c()
  df_gene <- ribotricer_aa %>% filter(gene_name == gene) %>% arrange(aa_length)
  if (length(df_gene$sequence) > 1) {
    for (i in 1:(length(df_gene$sequence)-1)) {
      for (j in 1:(length(df_gene$sequence)-i)) {
        if (str_detect(df_gene$sequence[i+j], df_gene$sequence[i])) {
          index <- c(index, i)
        }
      }
    }
    unique(index)
    if (is.numeric(index)) {
      longest_df_gene <- df_gene[-unique(index),]
    } else {
      longest_df_gene <- df_gene
    }
  }
  if (length(df_gene$sequence) == 1) {
    longest_df_gene <- df_gene
  }
  df_uniq_longest_orf <- rbind(df_uniq_longest_orf, longest_df_gene)
}
print("Crazy loop DONE!")

# Keep only "essential" columns
df_uniq_longest_orf <- df_uniq_longest_orf %>% select(1,2,3,4,7,19,21,22,23)

# Export table of unique longest translating ORFs
write.table(df_uniq_longest_orf, file = paste0(outdir, "uniq_longest_orf.tsv"), quote = F, row.names = F, sep = "\t")


#################################


