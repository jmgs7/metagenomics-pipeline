#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dada2))

# read seqtab.nochim.rds
seqtab.nochim <- readRDS("dada2/seqtab.nochim.rds")

taxa <- assignTaxonomy(seqtab.nochim, "references/silva_nr99_v138.1_train_set.fa.gz", multithread=10)
taxa <- addSpecies(taxa, "references/silva_species_assignment_v138.1.fa.gz")

saveRDS(taxa, "dada2/taxa.rds")