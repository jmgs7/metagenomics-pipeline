#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(glue))
suppressMessages(library(dada2))
suppressMessages(library(DECIPHER))
suppressMessages(library(decontam))
suppressMessages(library(phyloseq))
suppressMessages(library(DESeq2))
suppressMessages(library(vegan))

r1_files <- list.files("rawReads", pattern = "R1_001.fastq.gz", full.names = TRUE)

sample_ids <- gsub("_001.fastq.gz", "", r1_files)
sample_ids <- unlist(lapply(sample_ids, function(x) {
    basename(unlist(strsplit(x, "_"))[1])
}))

create_sample_sheet <- function(sample_id) {
    fastq_files <- list.files("rawReads", pattern = sample_id, full.names = TRUE)
    fastq_files <- unlist(lapply(fastq_files, function(f) {
        if (basename(unlist(strsplit(f, "_"))[1]) == sample_id) {
            return(f)
        }
    }))
    r1_fastq_files <- fastq_files[grep("R1", fastq_files)]
    r2_fastq_files <- fastq_files[grep("R2", fastq_files)]
    return(data.frame("sampleID" = sample_id, "forwardReads" = r1_fastq_files, "reverseReads" = r2_fastq_files))
}

samplesheet <- rbindlist(lapply(sample_ids, create_sample_sheet))

samplesheet

dir.create("usedRawReads", showWarnings = FALSE)

copy_files <- function(sample_id) {
    sample_sheet_x <- samplesheet[sampleID == sample_id]
    r1_fastq_files <- sample_sheet_x$forwardReads
    r2_fastq_files <- sample_sheet_x$reverseReads
    file.copy(r1_fastq_files, file.path("usedRawReads", glue("{sample_id}_R1.fastq.gz")), overwrite = TRUE)
    file.copy(r2_fastq_files, file.path("usedRawReads", glue("{sample_id}_R2.fastq.gz")), overwrite = TRUE)
    r1_fastq_files <- file.path("usedRawReads", glue("{sample_id}_R1.fastq.gz"))
    r2_fastq_files <- file.path("usedRawReads", glue("{sample_id}_R2.fastq.gz"))
    return(data.frame("sampleID" = sample_id, "forwardReads" = r1_fastq_files, "reverseReads" = r2_fastq_files))
}

samplesheet <- rbindlist(mclapply(sample_ids, copy_files, mc.cores = 10))

write.table(file = "sample_sheet.tsv", samplesheet, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

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

    colnames(fastqc_data_quality_df) <- fastqc_data_quality_df[1, ]
    fastqc_data_quality_df <- fastqc_data_quality_df[-1, ]
    fastqc_data_quality_df <- fastqc_data_quality_df %>%
        mutate(Median = as.numeric(Median), filename = basename(file)) %>%
        dplyr::rename(Base = "#Base") %>%
        mutate(
            Base = factor(Base, levels = unique(Base)),
            read = ifelse(grepl("R1", filename), "R1", "R2")
        ) %>%
        select(Base, Median, filename, read)
}

launch_fastqc <- function(sample_name, r1_file, r2_file) {
    r1_fastqc_file <- glue("{outdir}/{sample_name}_R1_fastqc.zip")
    r2_fastqc_file <- glue("{outdir}/{sample_name}_R2_fastqc.zip")
    if (!file.exists(r1_fastqc_file)) {
        system(glue("fastqc {r1_file} {r2_file} -o {outdir} -t 10 -q"))
    }
    fastqc_data <- rbindlist(lapply(c(r1_fastqc_file, r2_fastqc_file), read_fastqc))
    breaks_show <- unique(fastqc_data$Base)
    breaks_show <- split(breaks_show, ceiling(seq_along(breaks_show) / 8))
    breaks_show <- unlist(lapply(1:length(breaks_show), function(x) {
        shunk <- breaks_show[[x]]
        first <- shunk[1]
    }))

    aaa <- ggplot(fastqc_data, aes(x = Base, y = Median, color = read)) +
        geom_point() +
        scale_x_discrete(breaks = breaks_show) +
        theme_classic()
    fastqc_data$sample_name <- sample_name
    ggsave(glue("{outdir}/{sample_name}.pdf"), aaa, width = 10, height = 5, units = "in")
    write.table(fastqc_data, glue("{outdir}/{sample_name}_stats.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    return(fastqc_data)
}

outdir <- "fastqc"
dir.create(outdir, showWarnings = F)
outdir <- "fastqc/raw"
dir.create(outdir, showWarnings = F)

fastqc_stats <- rbindlist(mclapply(seq(1, nrow(samplesheet)), function(i) {
    sample_name <- samplesheet$sampleID[i]
    r1_file <- samplesheet$forwardReads[i]
    r2_file <- samplesheet$reverseReads[i]
    launch_fastqc(sample_name, r1_file, r2_file)
}, mc.cores = 10))

head(fastqc_stats)

# set threshold
qmin <- 25

# get median values for each cycle by read
median_values <- fastqc_stats %>%
    group_by(Base, read) %>%
    summarise(Median = median(Median), n = n()) %>%
    ungroup() %>%
    mutate(BaseNumMin = as.numeric(gsub("-.*", "", Base))) %>%
    mutate(BaseNumMax = as.numeric(gsub(".*-", "", Base)))

# the following code is useful to visualize the data
breaks_show <- unique(median_values$Base)
breaks_show <- split(breaks_show, ceiling(seq_along(breaks_show) / 8))
breaks_show <- unlist(lapply(1:length(breaks_show), function(x) {
    shunk <- breaks_show[[x]]
    first <- shunk[1]
}))

aaa <- ggplot(median_values, aes(x = Base, y = Median, color = read)) +
    geom_point() +
    scale_x_discrete(breaks = breaks_show) +
    theme_classic() +
    geom_hline(yintercept = qmin, color = "red")
ggsave("fastqc/raw/threshold.pdf", aaa, width = 10, height = 5, units = "in")


# First BaseNumMin of read R1 where Median is less than 25
idx_1 <- which(median_values$read == "R1" & median_values$Median < qmin)[1]
idx_2 <- which(median_values$read == "R2" & median_values$Median < qmin)[1]

# get the base number
# If all bases passed the phred score threshold, trim at least the last 10 bp just in case.
minimumTrimming <- max(median_values$BaseNumMax) - 10
base1 <- ifelse(is.na(idx_1), minimumTrimming, min(minimumTrimming, median_values$BaseNumMin[idx_1]))
base2 <- ifelse(is.na(idx_2), minimumTrimming, min(minimumTrimming, median_values$BaseNumMin[idx_2]))

threshold <- list("R1" = base1, "R2" = base2)


dir.create("trimmedReads", showWarnings = FALSE)

fwd_reads <- samplesheet$forwardReads
rev_reads <- samplesheet$reverseReads

fwd_post_reads <- file.path("trimmedReads", paste0(samplesheet$sampleID, "_R1.fastq.gz"))
rev_post_reads <- file.path("trimmedReads", paste0(samplesheet$sampleID, "_R2.fastq.gz"))

names(fwd_post_reads) <- samplesheet$sampleID
names(rev_post_reads) <- samplesheet$sampleID

if (!file.exists(fwd_post_reads[1])) {
    out <- filterAndTrim(fwd_reads, fwd_post_reads, rev_reads, rev_post_reads,
        truncLen = c(threshold$R1, threshold$R2),
        maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
        compress = TRUE, multithread = 10
    )
    head(out)
}

samplesheet$forwardReads <- fwd_post_reads
samplesheet$reverseReads <- rev_post_reads

outdir <- "fastqc/trimmed"
dir.create(outdir, showWarnings = F)

fastqc_stats <- rbindlist(mclapply(seq(1, nrow(samplesheet)), function(i) {
    sample_name <- samplesheet$sampleID[i]
    r1_file <- samplesheet$forwardReads[i]
    r2_file <- samplesheet$reverseReads[i]
    launch_fastqc(sample_name, r1_file, r2_file)
}, mc.cores = 10))

head(fastqc_stats)

derep_fn <- function(read) {
    fastq_files <- samplesheet %>%
        select(!!as.symbol(read)) %>%
        pull(!!as.symbol(read))
    derep <- derepFastq(fastq_files)
}

dir.create("rdata", showWarnings = F)

if (!file.exists(file.path("rdata", "derep.rdata"))) {
    derep <- mclapply(c("forwardReads", "reverseReads"), derep_fn, mc.cores = 2)
    names(derep) <- c("forwardReads", "reverseReads")
    save(derep, file = file.path("rdata", "derep.rdata"))
}

load(file.path("rdata", "derep.rdata"))

err_learning <- function(read) {
    errF <- learnErrors(derep[[read]], multithread = 10, nbases = 1e8)
}

if (!file.exists(file.path("rdata", "err.rdata"))) {
    error_vars <- lapply(c("forwardReads", "reverseReads"), err_learning)
    names(error_vars) <- c("forwardReads", "reverseReads")
    save(error_vars, file = file.path("rdata", "err.rdata"))
}

load(file.path("rdata", "err.rdata"))

sample_inference <- function(read) {
    dadaFs <- dada(derep[[read]], err = error_vars[[read]], multithread = 10)
}

if (!file.exists(file.path("rdata", "dada.rdata"))) {
    dada_reads <- lapply(c("forwardReads", "reverseReads"), sample_inference)
    names(dada_reads) <- c("forwardReads", "reverseReads")
    save(dada_reads, file = file.path("rdata", "dada.rdata"))
}

load(file.path("rdata", "dada.rdata"))

if (!file.exists(file.path("rdata", "merged.rdata"))) {
    mergers <- mergePairs(dada_reads$forwardReads, derep$forwardReads, dada_reads$reverseReads, derep$reverseReads, verbose = TRUE)
    head(mergers[[1]])
    save(mergers, file = file.path("rdata", "merged.rdata"))
}

load(file.path("rdata", "merged.rdata"))

if (!file.exists(file.path("rdata", "seqtab.rdata"))) {
    seqtab <- makeSequenceTable(mergers)
    save(seqtab, file = file.path("rdata", "seqtab.rdata"))
}

load(file.path("rdata", "seqtab.rdata"))


if (!file.exists(file.path("rdata", "seqtab.nochim.rdata"))) {
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 10, verbose = TRUE)
    save(seqtab.nochim, file = file.path("rdata", "seqtab.nochim.rdata"))
}

load(file.path("rdata", "seqtab.nochim.rdata"))

if (!file.exists(file.path("rdata", "taxa-info.rdata"))) {
    download.file(url = "http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile = file.path("rdata", "taxa_info.rdata"))
    load(file.path("rdata", "taxa_info.rdata"))
    dna <- DNAStringSet(getSequences(seqtab.nochim))
    tax_info <- IdTaxa(test = dna, trainingSet = trainingSet, strand = "both", processors = NULL)
    save(tax_info, file = file.path("rdata", "taxa-info.rdata"))
}

load(file.path("rdata", "taxa-info.rdata"))

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep = "\t", quote = F, col.names = NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern = ">", replacement = "", x = asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote = F, col.names = NA)

colnames(asv_tab)
vector_for_decontam <- c(rep(FALSE, 32), rep(TRUE, 1)) # TRUE is the possitive control (mock)

contam_df <- isContaminant(t(asv_tab), neg = vector_for_decontam)

table(contam_df$contaminant) # identified 7 as contaminants

contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[-dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
    sep = "\t", quote = F, col.names = NA
)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
    sep = "\t", quote = F, col.names = NA
)


#### phyloseq

count_tab <- read.table("ASVs_counts-no-contam.tsv", header = T, row.names = 1, check.names = F, sep = "\t")

# remove control samples
count_tab <- count_tab[, -grep("Control", colnames(count_tab))]
colnames(count_tab) <- gsub("_R1.fastq.gz","", colnames(count_tab))

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv",
    header = T,
    row.names = 1, check.names = F, sep = "\t"
))

sample_info_tab <- rbindlist(lapply(colnames(count_tab), function(col.name){
  shunks <- unlist(strsplit(col.name,"-"))
  last_shunk <- shunks[length(shunks)]
  f_letter <- last_shunk
  # f_letter <- unlist(strsplit(last_shunk,""))[1]
  data.frame(sample_name = col.name, group = f_letter)
})) %>% as.data.frame()

rownames(sample_info_tab) <- sample_info_tab$sample_name

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~group) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <- phyloseq(vst_count_phy, sample_data(sample_info_tab))


# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="group") + 
  geom_point(size=1) + 
  # labs(col="group") + 
  # geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_bw() +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# alpha diversity
rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(count_tab)))))


# Richness and diversity estimates
# first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_data(sample_info_tab))

# and now we can call the plot_richness() function on our phyloseq object
plot_richness(ASV_physeq, color="group", measures=c("Chao1", "Shannon")) +
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


plot_richness(ASV_physeq, x="group", color="group", measures=c("Chao1", "Shannon")) + 
  theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# using phyloseq to make a count table that has summed all ASVs
# that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="phylum")) 

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,"phylum"]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# we also have to account for sequences that weren't assigned any
# taxonomy even at the phylum level 
# these came into R as 'NAs' in the taxonomy table, but their counts are
# still in the count table
# so we can get that value for each sample by subtracting the column sums
# of this new table (that has everything that had a phylum assigned to it)
# from the column sums of the starting count table (that has all
# representative sequences)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in
# broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

# making count table broken down by class (contains classes beyond the
# Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class")) 

# making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)

# making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])

# changing the row names like above so that they correspond to the taxonomy,
# rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 

# making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# there are also possibly some some sequences that were resolved to the level
# of Proteobacteria, but not any further, and therefore would be missing from
# our class table
# we can find the sum of them by subtracting the proteo class count table
# from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

# and to check we didn't miss any other sequences, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 42 rows, which might be a little busy for a
# summary figure
# many of these taxa make up a very small percentage, so we're going to
# filter some out
# this is a completely arbitrary decision solely to ease visualization and
# intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
# now we have 12, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

## don't worry if the numbers or taxonomy vary a little, this might happen due to different versions being used 
## from when this was initially put together

# first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()

# take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
# manipulating tables like this is something you may need to do frequently in R

# change X from the start
filt_major_taxa_proportions_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g %>% mutate(Sample = gsub("\\.","-", gsub("X","", Sample)))

# now we want a table with "color" and "characteristics" of each sample to
# merge into our plotting table so we can use that more easily in our plotting
# function
# here we're making a new table by pulling what we want from the sample
# information table
sample_info_for_merge <- data.frame("Sample"=row.names(sample_info_tab), "Group"=sample_info_tab$group, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
# (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

# and now we're ready to make some summary figures with our wonderfully
# constructed table

## a good color scheme can be hard to find, i included the viridis package
## here because it's color-blind friendly and sometimes it's been really
## helpful for me, though this is not demonstrated in all of the following :/ 

# one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")


ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(Group), shape=factor(Group)), size=2, width=0.15, height=0) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")
