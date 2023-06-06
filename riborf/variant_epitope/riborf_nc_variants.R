#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
  cat("
How to use:
arg1: path to output_variants.tsv (output of riborf_vcfpy.py)
arg2: translating.candidateORF.genepred.txt
arg3: output dir

      ")
  quit()
}

outvar <- args[1]
candidates <- args[2]
output_dir <- args[3]

# Load results from variant annotation using translating orfs (riborf analysis) 
variants <- fread(outvar)
colnames(variants) <- c("locus", "ref_alt", "gene", "orfID", "type", "pass", "cdna_pos", "cds_pos", "protein_pos", "aa_change", "codon_change")
variants$orfID <- gsub("[&:]", "|", variants$orfID) #replace & and : signs by |

# Load riborf results about translating orfs with coordinates
tr_orf_candidate <- fread(candidates)
colnames(tr_orf_candidate)[1:7] <- c("orfID", "chr", "strand", "transcript_start", "transcript_end", "orf_start", "orf_end")
tr_orf_candidate$orfID <- gsub("[:]", "|", tr_orf_candidate$orfID)

# Merge to have orf coordinates in variant table
variants <- variants %>% left_join(tr_orf_candidate[,c(1,6:7)], by = "orfID")


# Split orfID column into multiple fields
variants <- variants %>% mutate(orf_ID = orfID) %>%  separate(., orf_ID,
                                                              c("transcript", 
                                                                "chrom", 
                                                                "strand", 
                                                                "RankNumber", 
                                                                "transcriptLength", 
                                                                "startCodonPosition", 
                                                                "stopCodonPosition", 
                                                                "candidateORFType", 
                                                                "startCodonType"), sep = '[|:]')

# Keep only variants without any canonical expression
nc_variants <- data.frame()

for (gene_x in unique(variants$gene)) {
  gene_df <- variants %>% filter(gene == gene_x)
  if (!("canonical" %in% gene_df$candidateORFType)) {
    nc_variants <- rbind(nc_variants, gene_df)
  }
}

# Filter for variants passing our filters
nc_pass <- nc_variants %>% filter(pass == "TRUE")

# Export table 
write.table(nc_pass, file=paste0(output_dir, "/non_canonical_variants.tsv"), quote = F, row.names = F, sep = "\t")