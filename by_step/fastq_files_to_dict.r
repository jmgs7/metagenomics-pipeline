#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(jsonlite))


# create command line parser
parser <- ArgumentParser()

# add argument fastq containing folder
parser$add_argument(
    "-f", "--fastq",
    help = "fastq containing folder"
)

# parse the command line arguments
args <- parser$parse_args()


# get list of fastq files
fq_files <- list.files(args$fastq, pattern = "*.fastq.gz", full.names = TRUE)

assign_names <- function(r1_file){
    basename(unlist(strsplit(r1_file, "_R1"))[1])
}

get_sample_names <- function(fq_files){
    r1_files <- fq_files[grep("_R1.fastq.gz", fq_files)]
    r1_names <- unlist(lapply(r1_files, assign_names))
    r2_files <- fq_files[grep("_R2.fastq.gz", fq_files)]
    fq_files <- lapply(seq_along(r1_names), function(x){
        list(R1=r1_files[x], R2=r2_files[x])
    })
    names(fq_files) <- r1_names
    return(fq_files)
}

fq_files <- get_sample_names(fq_files)
## write fq_files object in a json file
write_json(fq_files, glue("{args$fastq}/sample_ids.json"), pretty = TRUE)