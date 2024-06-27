#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(Biostrings))


ps <- readRDS("phyloseq/clean.rds")

# Remove mock samples
ps <- prune_samples(sample_data(ps)$status != "Control_positivo", ps) # Remove mock sample
ps <- prune_samples(sample_data(ps)$status != "Control_negativo", ps) # Remove negative control sample

# Renaming ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Compute alpha diversity

# Calculate the relevant alpha-diversity measures
erDF = estimate_richness(ps, split=TRUE, measures=c("Observed", "Shannon", "Simpson"))
# Measures may have been renamed in `erDF`. Replace it with the name from erDF
measures = colnames(erDF)
# Define "measure" variables and s.e. labels, for melting.
ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
# Remove any S.E. from `measures`
measures = measures[!measures %in% ses]

erDF <- erDF %>% rownames_to_column("sample") %>% 
    pivot_longer(cols = measures, names_to = "measure", values_to = "value") %>%
    # remove X at the beginning of sample name if it exists
    mutate(sample = gsub("^X", "", sample)) %>%
    # change dots in sample names to underscores
    mutate(sample = gsub("\\.", "-", sample)) %>%
    left_join(sample_data(ps), by = "sample")

alpha_plot <- ggplot(erDF, aes(x = status, y = value, fill = status, color = status)) +
    geom_boxplot(alpha = 0.2) +
    geom_point() +
    facet_wrap(~ measure, scales = "free") +
    theme_minimal()

dir.create("figures", showWarnings = FALSE)

ggsave("figures/alpha_diversity.pdf", alpha_plot, width = 10, height = 5, units = "in", bg="white")


# Compute beta diversity

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# Calcultae Bray-Curtis distances and convert to matrix
abrel_bray <- phyloseq::distance(ps.prop, method = "bray") %>% as.matrix()

# Store Bray Curtis distances between metagenomes from the same groups
sub_dist <- list()
groups_all <- sample_data(ps.prop)$status

for (group in unique(groups_all)) { 
    row_group <- which(groups_all == group)
    sample_group <- sample_names(ps.prop)[row_group]
    sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))

beta_plot <- ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
    geom_jitter() + 
    geom_boxplot(alpha=0.6) +  
    theme(legend.position="none") +
    ylab("Bray-Curtis diversity") +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))

ggsave("figures/beta_diversity.pdf", beta_plot, width = 10, height = 5, units = "in", bg="white")

# PcoA plot
ord = ordinate(ps.prop, method="PCoA", distance = "bray")
pcoa_plot <- plot_ordination(ps.prop, ord, color = "status") + 
  geom_point(size=4) + 
  stat_ellipse(aes(group=status))

ggsave("figures/pcoa.pdf", pcoa_plot, width = 10, height = 5, units = "in", bg="white")
