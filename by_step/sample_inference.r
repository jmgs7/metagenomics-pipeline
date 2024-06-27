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


sample_inferring <- function(read, fastq_files){
    fastq_files <- unlist(lapply(fastq_files, function(x){
        x[[read]]
    }))
    errF <- readRDS(glue("dada2/err_{read}.rds"))
    dadaFs <- dada(fastq_files, err=errF, multithread=10)
}


dada_object <- sample_inferring(args$read, fastq_files)
dir.create("dada2", showWarnings = FALSE)
saveRDS(dada_object, glue("dada2/inference_{args$read}.rds"))