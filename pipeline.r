#!/usr/bin/env Rscript

# This script is heavily based on https://github.com/ycl6/16S-Demo

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
suppressMessages(library(phangorn))
suppressMessages(library(ShortRead))
source("./lib/qc_analysis.r")
source("./lib/dada2.r")

folders.list <- list.dirs(recursive = FALSE)
fastq.folders.list <- folders.list[grepl("rawReads", folders.list)]

# cutadapt path
cutadapt <- "/home/jose/anaconda3/envs/metagenomics/bin/cutadapt" # cutadapt version 4.9
FWD <- "GTGYCAGCMGCCGCGGTAA" # 515F
REV <- "GGACTACNVGGGTWTCTAAT" # 806R
FWD.RC <- dada2::rc(FWD)
REV.RC <- dada2::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)


dbpath <- "db/"
if (!file.exists(dbpath)) dir.create(dbpath)

download_file <- function(url, destfile){
  if (!file.exists(destfile)){
    download.file(url, destfile = destfile)
  }
}

download.file("http://evomics.org/wp-content/uploads/2016/01/taxa_summary.R.gz", "taxa_summary.R.gz")
system("gunzip -f taxa_summary.R.gz")
source("taxa_summary.R", local = TRUE)

url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
download_file(url, destfile = paste0(dbpath, "/silva_nr99_v138.1_train_set.fa.gz"))

url <- "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz"
download_file(url, destfile = paste0(dbpath, "/silva_species_assignment_v138.1.fa.gz"))

url <- "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
download_file(url, destfile = paste0(dbpath, "/16SMicrobial.tar.gz"))
ref1 <- paste0(dbpath, "silva_nr99_v138.1_train_set.fa.gz")
ref2 <- paste0(dbpath, "silva_species_assignment_v138.1.fa.gz")
ref3 <- paste0(dbpath, "16SMicrobial.tar.gz")

create_fa_from_blast <- function(ref3 ){
  outfile <- paste0(dbpath, "16SMicrobial.fa.gz")
  if (!file.exists(outfile)) {
    system(glue("tar zxf {ref3} -C {dbpath}"))
    system(glue("/home/jose/anaconda3/envs/metagenomics/bin/blastdbcmd -db {dbpath}16S_ribosomal_RNA -entry all -out {dbpath}16SMicrobial.fa"))
    system(glue("gzip {dbpath}16SMicrobial.fa"))
  }
  return(outfile)
}
ref3 <- create_fa_from_blast(ref3)

lapply(fastq.folders.list, function(fastq.folder) {
    print(paste0("Processing: ", fastq.folder))
  
    # create folders
    folder.id <- unlist(strsplit(basename(fastq.folder), "_"))[1]
    trimmed.folder <- paste0(folder.id, "_trimmed")
    filt.folder <- paste0(folder.id, "_filtered")
    outfiles.folder <- paste0(folder.id, "_outfiles")
    images.folder <- paste0(folder.id, "_images")

    if (!dir.exists(trimmed.folder)) dir.create(trimmed.folder)
    if (!dir.exists(filt.folder)) dir.create(filt.folder)
    if (!dir.exists(outfiles.folder)) dir.create(outfiles.folder)
    if (!dir.exists(images.folder)) dir.create(images.folder)

    # list fastq files
    fns <- sort(list.files(fastq.folder, full.names = TRUE))
    fns <- fns[grep("MDMesa", fns)]
    fns <- fns[grep(".fastq.gz", fns)]
    fnFs <- fns[grep("_R1_", fns)]
    fnRs <- fns[grep("_R2_", fns)]
    sample.names <- basename(fnFs)
    sample.names <- unlist(lapply(sample.names, function(x) {
        y <- unlist(strsplit(x, "_MDMesa"))[1]
        return(y)
    }))

    # output of cutadapt
    fnFs.cut <- file.path(trimmed.folder, paste0(sample.names, "_R1.fastq.gz"))
    fnRs.cut <- file.path(trimmed.folder, paste0(sample.names, "_R2.fastq.gz"))

    log.cut <- gsub("_R1.fastq.gz", ".log", fnFs.cut)

    # launch cutadapt
    keep_res <- mclapply(seq_along(fnFs), function(i) {
        if (!file.exists(fnFs.cut[i])) {
            system2(cutadapt,
                stdout = log.cut[i], stderr = log.cut[i], # log file
                args = c(
                    R1.flags, R2.flags,
                    "-n 2", # -n 2 required to remove FWD and REV from reads
                    "--match-read-wildcards", # enable IUPAC nucleotide codes (wildcard characters)
                    "--length 300", # Truncate reads to 300 bp
                    "-m 150", # discard reads shorter than LEN (avoid length zero sequences)
                    "--overlap 10", # min overlap between read and adapter for an adapter to be found
                    "-j 0", # auto-detection of CPU cores, only available on Python 3
                    "--discard-untrimmed",
                    "-o", fnFs.cut[i], "-p", fnRs.cut[i], # trimmed files
                    fnFs[i], fnRs[i]
                ) # input files
            )
        }
    }, mc.cores = detectCores())

    # list trimmed files
    fns <- sort(list.files(trimmed.folder, full.names = TRUE))
    fnFs <- fns[grep("_R1.fastq.gz", fns)]
    fnRs <- fns[grep("_R2.fastq.gz", fns)]

    # create samplesheet
    sample_sheet <- data.frame(
        sampleID = sample.names,
        forwardReads = fnFs,
        reverseReads = fnRs
    )

    # qc and trunc estimation
    raw_fqc <- fastq_quality_control(sample_sheet, images.folder, raw = TRUE)
    trunc_parameters <- lapply(raw_fqc, trunc_estimation, min_qual = 30)
    cat(paste0("\nTrunc parameters: ", trunc_parameters, "\n"))

    # set filter files
    filtFs <- file.path(filt.folder, basename(fnFs))
    filtRs <- file.path(filt.folder, basename(fnRs))

    # launch filterAndTrim
    out <- filter_and_trim(fnFs, filtFs, fnRs, filtRs, trunc_parameters, outfiles.folder)

    # save out results
    out <- as.data.frame(out)
    rownames(out) <- sample.names
    write.table(file = "quality_trimming.tsv", out, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    # launch err_function
    errs <- err_function(filtFs, filtRs, outfiles.folder)
    errF <- errs[[1]]
    errR <- errs[[2]]

    # launch dada_function
    dadas <- dada_function(filtFs, filtRs, errF, errR, outfiles.folder)
    dadaFs <- dadas[[1]]
    dadaRs <- dadas[[2]]

    # launch mergepairs
    mergers <- merge_function(dadaFs, filtFs, dadaRs, filtRs, outfiles.folder)
    seqtab <- makeSequenceTable(mergers)
    
    # launch remove bimera
    seqtab.nochim = removebimera_function(seqtab, outfiles.folder)
    rownames(seqtab.nochim) = sample.names
    
    # launch decontam.
    # You need to adjust the number of FALSES and TRUES and their order according to you sample distribution.
    vector_for_decontam <-  grepl("control-negativo", rownames(seqtab.nochim), ignore.case = TRUE) # TRUE is the negative control.
    contam_df <- isContaminant(seqtab.nochim, neg = vector_for_decontam)
    contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
    seqtab.nochim.nocontam <- seqtab.nochim[,!colnames(seqtab.nochim) %in% contam_asvs]
    
    # assign tax
    taxatab <- assign_taxonomy(seqtab.nochim.nocontam, ref1, ref2, ref3, outfiles.folder)
    spec_silva <- taxatab[[2]]
    spec_ncbi <- taxatab[[3]]
    seqs <- taxatab[[4]]
    taxtab <- taxatab[[1]]

    ##############################################################################
    # lo que hace es un merge de addSpecies de Silva y de NCBI
    SVformat = paste("%0",nchar(as.character(ncol(seqtab.nochim.nocontam))),"d", sep = "")
    svid = paste0("ASV_", sprintf(SVformat, seq(ncol(seqtab.nochim.nocontam))))

    s_silva = as.data.frame(spec_silva, stringsAsFactors = FALSE)
    rownames(s_silva) = svid

    s_ncbi = as.data.frame(spec_ncbi, stringsAsFactors = FALSE)
    rownames(s_ncbi) = svid
    s_ncbi$Genus = gsub("\\[|\\]", "", s_ncbi$Genus)

    s_merged = cbind(s_ncbi, s_silva)
    colnames(s_merged) = c("nGenus","nSpecies","sGenus","sSpecies")
    s_merged1 = s_merged[!is.na(s_merged$nSpecies),]
    colnames(s_merged1)[1:2] = c("Genus","Species")
    s_merged2 = s_merged[is.na(s_merged$nSpecies) & !is.na(s_merged$sSpecies),]
    colnames(s_merged2)[3:4] = c("Genus","Species")
    s_merged3 = s_merged[is.na(s_merged$nSpecies) & is.na(s_merged$sSpecies),]
    colnames(s_merged3)[3:4] = c("Genus","Species")

    s_final = rbind(s_merged1[,c("Genus","Species")], s_merged2[,c("Genus","Species")],
                    s_merged3[,c("Genus","Species")])
    s_final = s_final[order(row.names(s_final)),]
    s_final = as.matrix(s_final)
    if("Genus" %in% colnames(taxtab$tax)) {
      gcol = which(colnames(taxtab$tax) == "Genus")
    } else {
      gcol = ncol(taxtab$tax)
    }

    matchGenera <- function(gen.tax, gen.binom, split.glyph = "/") {
      if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
      if((gen.tax == gen.binom) ||
         grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) ||
         grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }

    gen.match = mapply(matchGenera, taxtab$tax[,gcol], s_final[,1])
    taxtab$tax = cbind(taxtab$tax, s_final[,2])
    colnames(taxtab$tax)[ncol(taxtab$tax)] = "Species"
    print(paste(sum(!is.na(s_final[,2])), "out of",
                nrow(s_final), "were assigned to the species level."))
    taxtab$tax[!gen.match,"Species"] = NA
    print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])),
                "had genera consistent with the input table."))
    # #' 
    # #' # Multiple sequence alignment
    # #' 
    # #' Prepare a `data.frame` `df` from `seqtab.nochim` and `taxtab`.
    # #' 
    # ## ----prepare-df---------------------------------------------------------------
    df = data.frame(sequence = seqs, abundance = colSums(seqtab.nochim.nocontam), stringsAsFactors = FALSE)
    df$id = svid
    df = merge(df, as.data.frame(taxtab), by = "row.names")
    rownames(df) = df$id
    df = df[order(df$id),2:ncol(df)]
    
    # #' 
    # #' Performs alignment of multiple unaligned sequences.
    # #' 
    # ## ----alignseqs, cache = TRUE--------------------------------------------------
    
    align_seqs <- function(df){
      rds_align <- file.path(outfiles.folder, "align.rds")
      if (!file.exists(rds_align)){
        alignment = AlignSeqs(Biostrings::DNAStringSet(setNames(df$sequence, df$id)), anchor = NA, processors = detectCores())
        saveRDS(alignment, file = rds_align)
      }
      alignment <- readRDS(rds_align)
      return(alignment)
    }
    
    alignment <- align_seqs(df)
    
    # #' 
    # #' Export alignment.
    # #' 
    # ## ----export-alignment, cache = TRUE-------------------------------------------
    
    export_alignment <- function(alignment, outfiles.folder){
      rds_phang <- file.path(outfiles.folder,"phan.rds")
      if (!file.exists(rds_phang)){
        phang.align = phyDat(as(alignment, "matrix"), type = "DNA")
        saveRDS(phang.align, file = rds_phang)
      }
      phang.align <- readRDS(rds_phang)
      if (!file.exists(glue("{outfiles.folder}/alignment.aln"))){
        write.phyDat(phang.align, file = glue("{outfiles.folder}/alignment.fasta"), format = "fasta")
        write.phyDat(phang.align, file = glue("{outfiles.folder}/alignment.aln"), format = "phylip")
      }
      return(phang.align)
    }
    
    phang.align = export_alignment(alignment, outfiles.folder)
    
    
    # #' 
    # #' # Construct phylogenetic tree
    # #' 
    # #' Set the full-path to the RAxML and RAxML-NG.
    # #' 
    # ## ----set-raxml-path, eval = FALSE---------------------------------------------
    raxml = "/home/jose/anaconda3/envs/metagenomics/bin/raxmlHPC-PTHREADS-SSE3"		# e.g. /usr/local/bin/raxmlHPC-PTHREADS-SSE3
    raxmlng = "/home/jose/anaconda3/envs/metagenomics/bin/raxml-ng"	# e.g. /usr/local/bin/raxml-ng
    
    construct_tree <- function(raxml, raxmlng){
      raxml_tree <- glue("{outfiles.folder}/RAxML_binaryModelParameters.raxml_tree_GTRCAT")
      if (!file.exists(raxml_tree)){
        system2(raxml, args = c("-T $THREADS", "-f E", "-p 1234", "-x 5678", "-m GTRCAT", "-N 1",
                                glue("-s {outfiles.folder}/alignment.aln"), "-n raxml_tree_GTRCAT"))
        system(glue("mv RAxML* {outfiles.folder}"))
      }
      
      raxmlng_tree <- glue("{outfiles.folder}/GTRCAT.raxml.bestModel")
      if (!file.exists(raxmlng_tree)){
        system2(raxmlng, args = c("--evaluate", "--force", "--seed 1234", "--log progress", "--threads $THREADS",
                                  glue("--msa {outfiles.folder}/alignment.fasta"),  "--model GTR+G", "--brlen scaled",
                                  glue("--tree {outfiles.folder}/RAxML_fastTree.raxml_tree_GTRCAT"), "--prefix GTRCAT"))
        system(glue("mv GTRCAT* {outfiles.folder}"))
      }
      raxml_tree = read_tree(glue("{outfiles.folder}/GTRCAT.raxml.bestTree"))
      system(glue("cp {outfiles.folder}/GTRCAT.raxml.bestTree ./phylogenetic_tree.tree"))
      return(raxml_tree)
    }
    
    # #' 
    # #' Import tree using the `read_tree` function from `phyloseq`.
    # #' 
    # ## ----import-tree--------------------------------------------------------------
    raxml_tree = construct_tree(raxml, raxmlng)
    
    # Obtenemos la información de cada muestra a través de su nombre de archivo, que hemos modificado para que sean altamente informativos.
    filenames <- list.files(fastq.folder, pattern = "R1_001.fastq.gz", full.names = FALSE)
    filenames <- gsub("_R1_001.fastq.gz", "", filenames)
    
    sample_tab <- rbindlist(lapply(filenames, function(col.name){
      shunks <- unlist(strsplit(col.name,"_"))
      sample_name <- shunks[1]
      batch <- shunks[4]
      
      shunks <- unlist(strsplit(sample_name,"-"))
      sample <- shunks[1]
      sample_number <- shunks[2]
      time <- shunks[3]
      group <- shunks[length(shunks)]
      time_group = paste0(time, "_", group)
      data.frame(sample_name = sample_name, sample = sample, sample_number = sample_number, time = time, group = group, time_group = time_group, batch = batch)
    })) %>% as.data.frame()
    
    rownames(sample_tab) <- sample_tab$sample_name
    write.table(sample_tab, "sample_sheet.tsv", sep = "\t", quote = F, col.names = NA)

    # #' 
    # #' # Handoff to `phyloseq`
    # #' 
    # #' Prepare `new_seqtab` and `tax` data.frames containing abundance and taxonomy information respectively.
    # #' 
    # ## ----create-tax---------------------------------------------------------------
    new_seqtab = seqtab.nochim.nocontam
    colnames(new_seqtab) = df[match(colnames(new_seqtab), df$sequence),]$id
    write.table(t(new_seqtab), "ASVs_counts.tsv",
                sep = "\t", quote = F, col.names = NA
    )
    
    # Create taxonomy table.
    new_taxtab = taxtab
    rownames(new_taxtab$tax) = df[match(rownames(new_taxtab$tax), df$sequence),]$id
    tax = as.data.frame(new_taxtab$tax)
    tax$Family = as.character(tax$Family)
    tax$Genus = as.character(tax$Genus)
    write.table(tax, "ASVs_taxonomy.tsv",
                sep = "\t", quote = F, col.names = NA
    )
    write.table(new_taxtab, "ASVs_taxonomy_scores.tsv",
                sep = "\t", quote = F, col.names = NA
    )
      
    # giving our seq headers more manageable names (ASV_1, ASV_2...)
    asv_seqs <- colnames(seqtab.nochim.nocontam)
    asv_headers <- paste0(">", df[match(colnames(seqtab.nochim.nocontam), df$sequence),]$id)
    
    # making and writing out a fasta of our final ASV seqs:
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, "ASVs_sequences.fa")

    ps = phyloseq(tax_table(as.matrix(tax)),
                  sample_data(sample_tab),
                  otu_table(new_seqtab, taxa_are_rows = TRUE),
                  phy_tree(raxml_tree))
    save(ps, file = file.path(outfiles.folder, "phyloseq_object.rdata"))
})
