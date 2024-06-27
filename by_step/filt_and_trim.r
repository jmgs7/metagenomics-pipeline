#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(parallel))

# read usedRawReads/sample_ids.json
fastq_files <- fromJSON("usedRawReads/sample_ids.json")

# read threshold.json
trunc <- fromJSON("rawQC/threshold.json")

filt_and_trim <- function(sampleid, fastq_files, trunc){
    fastq_files <- fastq_files[[sampleid]]
    fwd_file <- fastq_files$R1
    rev_file <- fastq_files$R2

    fwd_file_out <- paste0("trimmedReads/", sampleid, "_R1.fastq.gz")
    rev_file_out <- paste0("trimmedReads/", sampleid, "_R2.fastq.gz")

    out <- filterAndTrim(
        fwd_file, fwd_file_out, rev_file, rev_file_out,
        truncLen = c(trunc$R1, trunc$R2),
        maxEE=c(2,2),rm.phix=TRUE,
        minLen=175, multithread = 10
    )
    out <- cbind(out, ID = row.names(out))
    return(list(R1=fwd_file_out, R2=rev_file_out))
}

dir.create("trimmedReads", showWarnings = FALSE)

trim_fq_files <- mclapply(names(fastq_files), filt_and_trim, fastq_files = fastq_files, trunc = trunc, mc.cores = 10)
names(trim_fq_files) <- names(fastq_files)
write_json(trim_fq_files, "trimmedReads/sample_ids.json", pretty = TRUE)