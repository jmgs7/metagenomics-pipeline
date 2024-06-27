#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dada2))

fastq_files <- fromJSON("trimmedReads/sample_ids.json")


fwd_reads <- unlist(lapply(names(fastq_files), function(x){
    fastq_files[[x]]$R1
}))

rev_reads <- unlist(lapply(names(fastq_files), function(x){
    fastq_files[[x]]$R2
}))

dada_fwd <- readRDS("dada2/inference_R1.rds")
dada_rev <- readRDS("dada2/inference_R2.rds")

mergers <- mergePairs(dada_fwd, fwd_reads, dada_rev, rev_reads)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)

saveRDS(seqtab.nochim, "dada2/seqtab.nochim.rds")
