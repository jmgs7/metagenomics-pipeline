#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(argparse))

# create command line parser
parser <- ArgumentParser()

# add read argument
parser$add_argument(
    "-r", "--read",
    help = "R1 or R2"
)

# parse the command line arguments
args <- parser$parse_args()

#read usedRawReads/sample_ids.json
fastq_files <- fromJSON("trimmedReads/sample_ids.json")


err_learning <- function(read, fastq_files){
    fastq_files <- unlist(lapply(fastq_files, function(x){
        x[[read]]
    }))
    errF <- learnErrors(fastq_files, multithread=10, nbases = 1e8)
}


err_learning_reads <- err_learning(args$read, fastq_files)
dir.create("dada2", showWarnings = FALSE)
saveRDS(err_learning_reads, glue("dada2/err_{args$read}.rds"))
