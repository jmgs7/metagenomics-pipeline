#' Get Quantile from a Distribution
#'
#' This function retrieves a quantile value from a cumulative distribution.
#'
#' @param xx A vector of values.
#' @param yy A vector of counts corresponding to the values.
#' @param q The quantile to retrieve.
#' @return The value corresponding to the specified quantile.
get_quant <- function(xx, yy, q) {
  xx[which(cumsum(yy) / sum(yy) >= q)][[1]]
}

#' Summarize Quality Control Metrics
#'
#' This function summarizes quality control metrics from sequencing data.
#'
#' @param f The file path to the FASTQ file.
#' @param n The number of reads to consider for the summary.
#' @return A list containing data frames for plotting, statistics, and annotations.
sum_qc <- function(f, n) {
  # Perform quality assessment on the file
  srqa <- qa(f, n = n)
  # Extract quality scores per cycle
  df <- srqa[["perCycle"]]$quality
  # Calculate mean quality scores per cycle
  means <- rowsum(df$Score * df$Count, df$Cycle) / rowsum(df$Count, df$Cycle)
  # Get total read counts
  rc <- sum(srqa[["readCounts"]]$read)
  # Calculate quantiles for quality scores per cycle
  q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify = TRUE)
  q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify = TRUE)
  q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify = TRUE)
  # Calculate cumulative counts per cycle
  cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify = TRUE)
  # Create a data frame for plotting
  plotdf <- df %>% mutate(file = basename(f))
  # Create a data frame for statistics
  statsdf <- data.frame(
    Cycle = as.integer(rownames(means)), Mean = means,
    Q25 = as.vector(q25s), Q50 = as.vector(q50s), Q75 = as.vector(q75s), Cum = 10 * as.vector(cums) / min(rc, n), file = basename(f)
  )
  # Create a data frame for annotations
  anndf <- data.frame(minScore = min(df$Score), label = basename(f), rc = rc, file = basename(f))
  return(list(plotdf = plotdf, statsdf = statsdf, anndf = anndf))
}

#' Analyze Read Quality
#'
#' This function analyzes the quality of reads from a sample sheet and generates summary statistics and plots.
#'
#' @param col.name The column name in the sample sheet to analyze (either "forwardReads" or "reverseReads").
#' @param samplesheet The sample sheet data frame.
#' @param result_folder The folder to save the results.
#' @param raw A boolean indicating whether the reads are raw or filtered.
#' @return A data frame with cycle, count, and median quality scores.
analyse_read <- function(col.name, samplesheet, result_folder, raw) {
  message(glue("Analyse {col.name}"))
  # Extract file paths from the sample sheet
  fls <- samplesheet[[col.name]]
  raw <- ifelse(raw, "raw", "filt")
  dest_file <- glue("{result_folder}/{raw}_{col.name}.tsv")
  n <- 5e+04
  if (!file.exists(dest_file)) {
    # Perform quality control on the files using parallel processing
    qc_res <- mclapply(fls, sum_qc, n = n, mc.cores = min(length(fls), 10))
    plotdf <- rbindlist(lapply(qc_res, function(x) x$plotdf))
    statsdf <- rbindlist(lapply(qc_res, function(x) x$statsdf))
    anndf <- rbindlist(lapply(qc_res, function(x) x$anndf))
    plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
    plotdf.summary$label <- col.name
    # Calculate summary statistics
    means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, plotdf.summary$Cycle) / rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify = TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify = TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify = TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), simplify = TRUE)
    statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), Q75 = as.vector(q75s), Cum = 10 * as.vector(cums) / sum(pmin(anndf$rc, n)))
    # Create a plot
    p <- ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) +
      geom_tile(aes(fill = Count)) +
      scale_fill_gradient(low = "#F5F5F5", high = "black") +
      geom_line(data = statdf.summary, aes(y = Mean), color = "#66C2A5") +
      geom_line(data = statdf.summary, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") +
      geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", size = 0.25) +
      geom_line(data = statdf.summary, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") +
      ylab("Quality Score") +
      xlab("Cycle") +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      guides(fill = FALSE) +
      facet_wrap(~label)
    data <- p$data
    # Create a data frame with cycle, count, and median quality scores
    df <- data.frame(Cycle = character(), Count = character(), Median = character(), stringsAsFactors = FALSE)
    cycles <- sort(unique(data$Cycle))
    for (cycle in cycles) {
      subdata <- data[data[, "Cycle"] == cycle, ]
      score <- list()
      # Convert to list to calculate median
      for (j in 1:nrow(subdata)) {
        score <- unlist(c(score, rep(subdata$Score[j], subdata$Count[j])))
      }
      temp <- data.frame(Cycle = cycle, Count = sum(subdata$Count), Median = median(score), stringsAsFactors = FALSE)
      df <- rbind(df, temp)
    }
    # Save the plot and data frame
    file_name <- gsub("\\.tsv$", ".pdf", dest_file)
    ggsave(filename = file_name, plot = p)
    write.table(df, file = dest_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # Read and return the data frame
  df <- read.delim(dest_file, stringsAsFactors = FALSE)
  return(df)
}

#' Perform Quality Control on FASTQ Files
#'
#' This function performs quality control on the forward and reverse reads in a sample sheet.
#'
#' @param samplesheet The sample sheet data frame.
#' @param result_folder The folder to save the results.
#' @param raw A boolean indicating whether the reads are raw or filtered.
#' @return A list of data frames with quality control results for forward and reverse reads.
fastq_quality_control <- function(samplesheet, result_folder, raw) {
  # Analyze forward and reverse reads
  res_qc <- lapply(c("forwardReads", "reverseReads"), analyse_read, samplesheet = samplesheet, result_folder = result_folder, raw = raw)
  names(res_qc) <- c("forwardReads", "reverseReads")
  return(res_qc)
}

#' Estimate Truncation Cycle
#'
#' This function estimates the truncation cycle for quality filtering based on a specified minimum quality score.
#'
#' @param fqc_table A data frame containing quality control metrics.
#' @param min_qual The minimum quality score to consider for truncation. Default is 25.
#' @return The estimated truncation cycle.
trunc_estimation <- function(fqc_table, min_qual = 25) {
  # Filter the quality control table for cycles with median quality below or equal to the minimum quality
  fqc_table_filt <- fqc_table %>%
    filter(Median <= min_qual)
  # If no rows are found after filtering, return the maximum cycle minus 10
  if (nrow(fqc_table_filt) == 0) {
    return(max(fqc_table$Cycle) - 10)
  } else {
    # Otherwise, return the minimum cycle from the filtered, but it will always trim at least the last 10 positions.
    return(min(max(fqc_table$Cycle) - 10, min(fqc_table_filt$Cycle)))
  }
}