#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)==0 | args[1]=="-h" | args[1]=="--help") {
  cat("
How to use:
arg1: path to gtf (protein coding gencode)
args2: output dir

      ")
  quit()
}

# Variables
gtf <- args[1]
output_dir <- args[2]

# Libraries
suppressPackageStartupMessages(library(riboWaltz))

annotation_dt = create_annotation(gtfpath = gtf)

saveRDS(annotation_dt, file = paste0(output_dir, "annotation_ribowaltz.rds"))