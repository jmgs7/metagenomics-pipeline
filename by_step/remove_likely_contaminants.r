#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(decontam))


ps <- readRDS("phyloseq/raw.rds")

# 'is.neg' is a logical vector in the sample_data indicating which samples are mock communities
sample_data(ps)$is.neg <- sample_data(ps)$status == "Control_positivo"

# Applying the prevalence method
contamdf <- isContaminant(ps, method="prevalence", neg="is.neg")

# Filter out contaminants
ps.clean <- prune_taxa(!contamdf$contaminant, ps)
saveRDS(ps.clean, "phyloseq/clean.rds")
