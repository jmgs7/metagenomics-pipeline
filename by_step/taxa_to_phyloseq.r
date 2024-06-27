#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))

# create command line parser
parser <- ArgumentParser()

# add pheno table argument
parser$add_argument(
    "-p", "--pheno",
    help = "pheno table"
)

# parse the command line arguments
args <- parser$parse_args()

samdf <- fread(args$pheno, header = TRUE, stringsAsFactors = FALSE) %>% as.data.frame()

seqtab.nochim <- readRDS("dada2/seqtab.nochim.rds")
taxa <- readRDS("dada2/taxa.rds")

rownames(samdf) <- samdf$sample

samdf <- samdf[rownames(seqtab.nochim), ]

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                sample_data(samdf),
                tax_table(taxa))
dir.create("phyloseq", showWarnings = FALSE)
saveRDS(ps, "phyloseq/raw.rds")