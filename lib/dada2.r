suppressMessages(library(dada2))

filter_and_trim <- function(fnFs, filtFs, fnRs, filtRs, trunc_parameters, outfiles.folder) {
    rds_file <- file.path(outfiles.folder, "filt_and_trim.rds")
    if (!file.exists(rds_file)) {
        # Perform filtering and trimming
        out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
            # Need to keep paramters consistent between runs of the same study
            truncLen = c(trunc_parameters$forwardReads, trunc_parameters$reverseReads),
            maxN = 0, truncQ = 2, maxEE = c(2, 2),
            rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE
        )
        saveRDS(out, file = rds_file)
    }
    out <- readRDS(rds_file)
    return(out)
}

err_function <- function(filtFs, filtRs, outfiles.folder){
  rds_err_r1 <- file.path(outfiles.folder, "err_r1.rds")
  rds_err_r2 <- file.path(outfiles.folder, "err_r2.rds")
  if (!file.exists(rds_err_r1) | !file.exists(rds_err_r2)) {
    errF = learnErrors(filtFs, multithread = TRUE)
    errR = learnErrors(filtRs, multithread = TRUE)
    saveRDS(errF, file = rds_err_r1)
    saveRDS(errR, file = rds_err_r2)
  }
  errF = readRDS(rds_err_r1)
  errR = readRDS(rds_err_r2)
  return(list(errF, errR))
}

dada_function <- function(filtFs, filtRs, errF, errR, outfiles.folder){
  rds_dadaFs <- file.path(outfiles.folder, "dadaFs.rds")
  rds_dadaRs <- file.path(outfiles.folder, "dadaRs.rds")
  if (!file.exists(rds_dadaFs) | !file.exists(rds_dadaRs)) {
    dadaFs = dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
    dadaRs = dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)
    saveRDS(dadaFs, file = rds_dadaFs)
    saveRDS(dadaRs, file = rds_dadaRs)
  }
  dadaFs = readRDS(rds_dadaFs)
  dadaRs = readRDS(rds_dadaRs)
  return(list(dadaFs, dadaRs))
}

merge_function <- function(dadaFs, filtFs, dadaRs, filtRs, outfiles.folder){
  rds_merge <- file.path(outfiles.folder, "mergers.rds")
  if (!file.exists(rds_merge)) {
    mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
    saveRDS(mergers, file = rds_merge)
  }
  mergers = readRDS(rds_merge)
  return(mergers)
}

removebimera_function <- function(seqtab, outfiles.folder){
  rds_nochim <- file.path(outfiles.folder, "seqtab.nochim.rds")
  if (!file.exists(rds_nochim)) {
    seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE,
                                       verbose = TRUE)
    saveRDS(seqtab.nochim, file = rds_nochim)
  }
  seqtab.nochim = readRDS(rds_nochim)
  return(seqtab.nochim)
}

assign_taxonomy <- function(seqtab, ref1, ref2, outfiles.folder){
  seqs = getSequences(seqtab)
  rds_taxa <- file.path(outfiles.folder, "taxa.rds")
  if (!file.exists(rds_taxa)) {
    set.seed(12345)
    taxatab = assignTaxonomy(seqs, refFasta = ref1, minBoot = 80, tryRC = TRUE, verbose = TRUE,
                            outputBootstraps = FALSE, multithread = TRUE)
    taxatab = addSpecies(taxatab, refFasta = ref2, allowMultiple = FALSE, tryRC = TRUE)
    saveRDS(taxatab, file = rds_taxa)
  }
  taxatab = readRDS(rds_taxa)
  return(list(taxatab, seqs))
}
