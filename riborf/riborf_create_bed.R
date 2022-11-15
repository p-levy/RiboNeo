#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)==0 | args[1]=="-h" | args[1]=="--help") {
  cat("
How to use:
arg1: path to riborf repre.valid.pred.pvalue.parameters.txt
arg2: path to translating.candidateORF.genepred.txt
arg3: output prefix
arg4: number of cores to use for parallel loop

      ")
  quit()
}
 
input_df <- args[1]
candidates <- args[2]
output_prefix <- args[3]
n.cores <- args[4]

# Load input files as data tables
total_orfs <- fread(input_df)
candidates_tr <- fread(candidates) %>% select(1, c(8:10)) %>% rename("orfID" = "V1", "exon_start"="V9", "exon_end"="V10", "n_exon"="V8")
# Merge with candidate ORF file to get coordinates of each ORF
total_orfs <- total_orfs %>% inner_join(candidates_tr, by = "orfID")


# Split orfID column into multiple fields
total_orfs <- total_orfs %>% mutate(orf_ID = orfID) %>%  separate(., orfID,
                                      c("transcript", 
                                        "chrom", 
                                        "strand", 
                                        "RankNumber", 
                                        "transcriptLength", 
                                        "startCodonPosition", 
                                        "stopCodonPosition", 
                                        "candidateORFType", 
                                        "startCodonType"), sep = '[|:]')

# Separate df for canonical and non-canonical
canonical_orfs <- total_orfs[total_orfs$candidateORFType == "canonical",]
non_canonical_orfs <- total_orfs[total_orfs$candidateORFType != "canonical",]

#create the cluster (for parallel creation of the bed files)
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


# Function: Loop to create bed file for canonical ORFs
bedLoop <- function(x) {
  bed <- data.table()
  bed <- foreach(i = x$orf_ID, .combine = "rbind") %dopar% {
    library(dplyr)
    library(tidyr)
    df <- x %>% filter(orf_ID == i)
    c <- mapply(paste,unlist(strsplit(as.character(df$exon_start), ",")),unlist(strsplit(as.character(df$exon_end), ",")), sep="-")
    if (df$n_exon > 1){
      df <- rbind(df, df[rep(1, n_exon-1), ])
    }
    df <- df %>% mutate(coordinate = c)
    df$Score <- rep(0, nrow(df))
    df <- df %>% separate(., coordinate, c("start", "end"), sep = '-')
    df <- df %>% select(2, 27, 28, 26, 29, 3)
  }
  write.table(bed, paste0(output_prefix, "_", deparse(substitute(x)), "_canonical_orfs.bed"), sep = "\t", row.names=F, col.names=F, quote=F)
}

# Create bed file for canonical ORFs
bedLoop(canonical_orfs)

# Create bed file for non-canonical ORFs
bedLoop(non_canonical_orfs)

stopCluster(my.cluster)
