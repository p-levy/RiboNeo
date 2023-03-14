#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# Usage
if (length(args)==0 || args[1]=="-h" || args[1]=="--help") {
  cat("
How to use:
arg1: patient name/code (should start with a letter, not a number)
arg2: path to directory where bam files are
arg3: path to ribowaltz annotation.rds
arg4: lower length filter
arg5: upper length filter
arg6: output dir

      ")
  quit()
}

# Variables
patient <- args[1]
bam_dir <- args[2]
annotation <- args[3]
low_filter <- as.numeric(args[4])
up_filter <- as.numeric(args[5])
output_dir <- args[6]


# Libraries
suppressPackageStartupMessages(library(riboWaltz))

annotation_dt <- readRDS(annotation)

# Load reads
## We first define the "name_samples" character string vector as follow:
name_of_bams <- c(paste0("sample_", patient))
names(name_of_bams) <- c("Aligned.toTranscriptome.out")

## Acquiring input files
unfiltered_reads_list <- bamtolist(
  bamfolder = bam_dir,
  name_samples = name_of_bams,
  annotation = annotation_dt)

reads_list_filt = length_filter(
  data = unfiltered_reads_list,
  length_filter_mode = "custom",
  length_range = low_filter:up_filter
)

# P-sites calculation
# Unfiltered
psite_offset = psite(unfiltered_reads_list)
reads_psite_list = psite_info(unfiltered_reads_list, psite_offset)
# length-filtered list
psite_offset_filt = psite(reads_list_filt)
reads_psite_list_filt = psite_info(reads_list_filt, psite_offset_filt)

# Plots
# Reads length distribution
# Unfiltered
length_dist_1 = rlength_distr(unfiltered_reads_list, sample = paste0("sample_", patient))
pdf(file = paste0(output_dir, "/length_dist.pdf"))
length_dist_1$plot
dev.off()
# length-filtered
length_dist_2 = rlength_distr(reads_list_filt, sample = paste0("sample_", patient))
pdf(file = paste0(output_dir, "/length_dist_filt.pdf"))
length_dist_2$plot
dev.off()

# Ends heatmaps
# Unfiltered
ends_heatmap = rends_heat(
  unfiltered_reads_list, 
  annotation_dt, 
  sample = paste0("sample_", patient), 
  cl = 95, # specifying a confidence level for
#'   restricting the plot to a sub-range of read lengths. The new range is
#'   associated to the most abundant populations of reads accounting for the cl%
#'   of the sample. Default is 95.
  utr5l = 25, 
  cdsl = 40, 
  utr3l = 25
)
pdf(file = paste0(output_dir, "/ends_heatmap.pdf"))
ends_heatmap$plot 
dev.off()

# length-filtered
ends_heatmap_filt = rends_heat(
  reads_list_filt, 
  annotation_dt, 
  sample = paste0("sample_", patient), 
  cl = 100, # all read lengths plotted
  utr5l = 25, 
  cdsl = 40, 
  utr3l = 25
)
pdf(file = paste0(output_dir, "/ends_heatmap_filt.pdf"))
ends_heatmap_filt$plot
dev.off()

# P-sites per region
# length-filtered
psite_region_1 = region_psite(reads_psite_list_filt, annotation_dt, sample = paste0("sample_", patient))
pdf(file = paste0(output_dir, "/psite_region_filt.pdf"))
psite_region_1$plot
dev.off()

# Trinucleotide periodicity
# unfiltered
frames_stratified = frame_psite_length(
  reads_psite_list, 
  sample = paste0("sample_", patient),
  region = "all", 
  cl = 95
)
pdf(file = paste0(output_dir, "/frames_stratified.pdf"))
frames_stratified$plot
dev.off()
# length-filtered
frames_stratified_1 = frame_psite_length(
  reads_psite_list_filt, 
  sample = paste0("sample_", patient),
  region = "all", 
  cl = 100
)
pdf(file = paste0(output_dir, "/frames_stratified_filt.pdf"))
frames_stratified_1$plot
dev.off()

frames_1 = frame_psite(reads_psite_list_filt, sample = paste0("sample_", patient), region = "all")
pdf(file = paste0(output_dir, "/frame_psite_filt.pdf"))
frames_1[["plot"]]
dev.off()

# Metaplots
# length-filtered
metaprofile = metaprofile_psite(
  reads_psite_list_filt, 
  annotation_dt, 
  sample = paste0("sample_", patient),
  utr5l = 20,
  cdsl = 40,
  utr3l = 20,
  plot_title = paste0(patient, " transcript")
)
pdf(file = paste0(output_dir, "/metaprofile_filt.pdf"))
metaprofile[[paste0("plot_", patient)]]
dev.off()

