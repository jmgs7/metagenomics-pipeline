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

# add argument for outdir
parser$add_argument(
    "-o", "--outdir",
    help = "outdir"
)

# parse the command line arguments
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE)


# read json file
fq_files <- fromJSON(glue("{args$fastq}/sample_ids.json"))

read_fastqc <- function(file) {

    # list the folder inside zip file
    zip_list <- unzip(file, list = TRUE)

    # path to fastqc_data.txt
    fastqc_data <- zip_list$Name[grepl("fastqc_data.txt", zip_list$Name)]

    # read the fastqc_data.txt inside the zip without decompress
    fastqc_data <- readLines(unz(file, fastqc_data))

    # we need only the module Per base sequence quality
    # to do this we keep all lines from ">>Per base sequence quality"
    # and end when we found a line with ">>END_MODULE"

    index_0 <- which(grepl(">>Per base sequence quality", fastqc_data))
    index_1 <- which(grepl(">>END_MODULE", fastqc_data))
    index_1 <- min(index_1[index_1 > index_0])

    fastqc_data_quality <- fastqc_data[(index_0 + 1):(index_1 - 1)]

    fastqc_data_quality_df <- rbindlist(lapply(1:length(fastqc_data_quality), function(i) {
        line <- fastqc_data_quality[i]
        line <- unlist(strsplit(line, "\t"))
        # generate dataframe
        data.frame(t(line))
    })) %>% as.data.frame()

    colnames(fastqc_data_quality_df) <- fastqc_data_quality_df[1,]
    fastqc_data_quality_df <- fastqc_data_quality_df[-1,]
    fastqc_data_quality_df <- fastqc_data_quality_df %>%
        mutate(Median = as.numeric(Median), filename = basename(file)) %>%
        rename(Base = "#Base") %>% mutate(Base = factor(Base, levels = unique(Base)),
        read = ifelse(grepl("R1", filename), "R1", "R2")) %>%
        select(Base, Median, filename, read)
}


# run fastqc
save_res <- mclapply(names(fq_files), function(x){
    r1_file <- fq_files[[x]]$R1
    r2_file <- fq_files[[x]]$R2
    system(glue("fastqc {r1_file} {r2_file} -o {args$outdir} -t 10"))
    r1_file <- glue("{args$outdir}/{x}_R1_fastqc.zip")
    r2_file <- glue("{args$outdir}/{x}_R2_fastqc.zip")
    fastqc_data <- rbindlist(lapply(c(r1_file, r2_file), read_fastqc))
    breaks_show <- unique(fastqc_data$Base)
    breaks_show <- split(breaks_show, ceiling(seq_along(breaks_show)/8))
    breaks_show <- unlist(lapply(1:length(breaks_show), function(x){
        shunk <- breaks_show[[x]]
        first <- shunk[1]
    }))

    aaa <- ggplot(fastqc_data, aes(x = Base, y = Median, color = read)) + geom_point() + scale_x_discrete(breaks = breaks_show) + theme_minimal()
    ggsave(glue("{args$outdir}/{x}.pdf"), aaa, width = 10, height = 5, units = "in")
    write.table(fastqc_data, glue("{args$outdir}/{x}_stats.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    return(glue("{args$outdir}/{x}_stats.tsv"))
}, mc.cores = 10)

names(save_res) <- names(fq_files)

## write save_res object in a json file
write_json(save_res, glue("{args$outdir}/fastqc_results.json"), pretty = TRUE)
