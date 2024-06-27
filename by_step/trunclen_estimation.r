#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))


# read rawQC/fastqc_results.json
fastqc_results <- fromJSON("rawQC/fastqc_results.json")


# concat all tables
fastqc_stats <- rbindlist(lapply(fastqc_results, function(x){
    fread(x)
})) %>% mutate(Base = factor(Base, unique(Base)))

# set threshold
qmin <- 25

# get median values for each cycle by read
median_values <- fastqc_stats %>% group_by(Base, read) %>% summarise(Median = median(Median), n=n()) %>% ungroup() %>% 
    mutate(BaseNum = as.numeric(gsub("-.*", "", Base)))

breaks_show <- unique(median_values$Base)
breaks_show <- split(breaks_show, ceiling(seq_along(breaks_show)/8))
breaks_show <- unlist(lapply(1:length(breaks_show), function(x){
    shunk <- breaks_show[[x]]
    first <- shunk[1]
}))

aaa <- ggplot(median_values, aes(x = Base, y = Median, color = read)) + geom_point() + scale_x_discrete(breaks = breaks_show) + theme_minimal() + geom_hline(yintercept = qmin, color = "red")
ggsave("rawQC/threshold.pdf", aaa, width = 10, height = 5, units = "in")

# First Basenum of read R1 where Median is less than 25
idx_1 <- which(median_values$read == "R1" & median_values$Median < qmin)[1]
idx_2 <- which(median_values$read == "R2" & median_values$Median < qmin)[1]

# get the base number
base1 <- ifelse(is.na(idx_1), max(median_values$BaseNum), median_values$BaseNum[idx_1])
base2 <- ifelse(is.na(idx_2), max(median_values$BaseNum), median_values$BaseNum[idx_2])

# list the thresholds
threshold <- list("R1"=base1, "R2"=base2)

# save the threshold in json
write_json(threshold, "rawQC/threshold.json")
