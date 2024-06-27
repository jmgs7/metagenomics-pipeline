#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(argparse)
library(data.table)
library(parallel)

r1_files <- list.files("rawReads", pattern = "R1_001.fastq.gz", full.names = TRUE)

sample_ids <- gsub("_001.fastq.gz", "", r1_files)
sample_ids <- unlist(lapply(sample_ids, function(x) {
    basename(unlist(strsplit(x, "_"))[1])
}))

# Parse command line arguments
parser <- ArgumentParser()

# Required arguments
parser$add_argument("-f", "--filter",action='store_true', help = "filter")

# Parse the arguments
args <- parser$parse_args()

samplesheet <- rbindlist(mclapply(sample_ids, function(x) {
    fastq_files <- list.files("rawReads", pattern = x, full.names = TRUE)
    fastq_files <- unlist(lapply(fastq_files, function(f) {
        if (basename(unlist(strsplit(f, "_"))[1]) == x){
            return(f)
        }        
    }))
    r1_fastq_files <- fastq_files[grep("R1", fastq_files)]
    r2_fastq_files <- fastq_files[grep("R2", fastq_files)]
    return(data.frame("sampleID" = x, "forwardReads" = r1_fastq_files, "reverseReads" = r2_fastq_files))
}, mc.cores = 10))

if (args$filter) {
    # select random samples
    sample_sheet <- samplesheet %>% head(5)
    sample_ids <- sample_sheet$sampleID
}

dir.create("usedRawReads", showWarnings = FALSE)

samplesheet <- rbindlist(mclapply(sample_ids, function(x) {
    sample_sheet_x <- samplesheet[sampleID == x]
    r1_fastq_files <- sample_sheet_x$forwardReads
    r2_fastq_files <- sample_sheet_x$reverseReads
    file.copy(r1_fastq_files, file.path("usedRawReads", glue("{x}_R1.fastq.gz")), overwrite = TRUE)
    file.copy(r2_fastq_files, file.path("usedRawReads", glue("{x}_R2.fastq.gz")), overwrite = TRUE)
    r1_fastq_files <- file.path("usedRawReads", glue("{x}_R1.fastq.gz"))
    r2_fastq_files <- file.path("usedRawReads", glue("{x}_R2.fastq.gz"))
    return(data.frame("sampleID" = x, "forwardReads" = r1_fastq_files, "reverseReads" = r2_fastq_files))
}, mc.cores = 10))

write.table(file="sample_sheet.tsv", samplesheet, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
