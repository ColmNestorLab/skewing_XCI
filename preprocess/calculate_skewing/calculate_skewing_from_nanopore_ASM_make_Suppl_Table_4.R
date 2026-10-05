#!/usr/bin/env Rscript

library(NanoMethViz)
library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)
library(GenomicRanges)

source("/home/bjogy93/Desktop/XXHealth/ALentini_ggplot_functions.R")

#renv::install("bioc::NanoMethViz", prompt = F)
#renv::install("tidyverse", prompt = F)
#renv::install("doParallel", prompt = F)
#renv::install("foreach", prompt = F)


#### Defining functions ####
apply_cluster_reads_parallel <- function(mbr, bed, min_pts, num_cores = 4) {
  # Initialize a parallel backend
  #cl <- makeCluster(detectCores())
  #try with 4 cores first
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Create a progress bar
  pb <- progress_estimated(nrow(bed))
  
  # Define a function to process each row
  process_row <- function(row) {
    # Get the current row
    current_row <- bed[row, ]
    
    # Apply the cluster_reads function to the current row
    row_cluster <- tryCatch({
      NanoMethViz:::cluster_reads(mbr, current_row$chr, current_row$start, current_row$end, min_pts = min_pts)},
      error = function(err){
        # Handle the error (e.g., print a message)
        message(paste("Error in row", row, ":", err$message))
        # Skip to the next row
        return(NA)
      })
    
    if (all(is.na(row_cluster))) {
      return(NULL)
    }
    
    # Add the CGI_id and chr, start, end
    row_cluster <- row_cluster %>% mutate(CGI_id = paste0(current_row$chr, ":", current_row$start, "-", current_row$end), chr = current_row$chr, start = current_row$start, end = current_row$end)
    
    # Calculate the average methylation by cluster_id and add it as a new column
    row_cluster <- row_cluster %>% group_by(cluster_id) %>% mutate(avg_cluster_methylation = mean(mean))
    
    # If there are exactly 2 clusters, assign the cluster with the lowest average methylation to Xa and the other to Xi
    if (nlevels(row_cluster$cluster_id) == 2) { # Watch out for NA clusters...
      low_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% min()
      high_mC <- row_cluster %>% filter(cluster_id %in% c("1", "2")) %>% pull(avg_cluster_methylation) %>% max()
      row_cluster <- row_cluster %>% mutate(assigned_X = case_when(avg_cluster_methylation == low_mC ~ "Xa", avg_cluster_methylation == high_mC ~ "Xi", TRUE ~ "NA"))
    } else {
      row_cluster$assigned_X <- NA
    }
    
    return(row_cluster)
  }
  
  # Apply the function to each row in parallel
  res <- foreach(row = 1:nrow(bed), .combine = bind_rows, .packages = c("tidyverse", "NanoMethViz")) %dopar% {
    pb$tick()  # Update the progress bar
    process_row(row)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Combine the results
  res <- bind_rows(res)
  
  return(res)
}

calculate_skew_by_block <- function(clustered_reads, haplotyped_reads){
  #remove uninformative reads that don’t clusters or reads from CGIs that don’t have exactly 2 clusters
  clustered_reads <- clustered_reads %>% filter(assigned_X %in% c("Xa","Xi"))
  #remove reads that appear multiple times because they span several CGIs
  clustered_reads <- clustered_reads %>% distinct(read_name, .keep_all = TRUE)
  #merge methylation cluster information with haplotype and phase set information
  df2 <- left_join(clustered_reads,haplotyped_reads, relationship = "many-to-many")
  #remove the reads that couldn’t be haplotyped
  df2 <- df2 %>% filter(!is.na(HP))
  
  #count by haplotype blocks
  counts_by_block <- df2 %>% group_by(PS, assigned_X, HP) %>% summarise(counts = n())
  
  skew_by_block <- counts_by_block %>%
    unite(combi, assigned_X, HP) %>%
    mutate(combi = recode(combi, "Xa_1" = "H1_Xa", "Xa_2" = "H2_Xa", "Xi_1" = "H1_Xi", "Xi_2" = "H2_Xi")) %>%
    pivot_wider(id_cols = PS, names_from = combi, values_from = counts, values_fill = 0) %>%
    mutate(H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi),
           totalCount = H1_Xa + H1_Xi + H2_Xa + H2_Xi)
  
  return(skew_by_block)
}

compute_folded_skew_with_CI <- function(df) {
  
  k <- df$H1_Xa + df$H2_Xi
  n <- df$H1_Xa + df$H2_Xa + df$H1_Xi + df$H2_Xi
  
  logLik_fun <- function(P) {
    
    if(P <= 0 || P >= 0.5) return(-Inf)
    
    log_term1 <- dbinom(k, n, P, log=TRUE)
    log_term2 <- dbinom(k, n, 1-P, log=TRUE)
    
    sum(log(exp(log_term1) + exp(log_term2)))
  }
  
  # MLE
  opt <- optimize(logLik_fun, interval=c(1e-8, 0.5-1e-8), maximum=TRUE)
  P_hat <- opt$maximum
  logLik_max <- opt$objective
  
  # Likelihood cutoff
  cutoff <- logLik_max - 0.5 * qchisq(0.95, df=1)
  
  # CI bounds
  lower <- uniroot(function(P) logLik_fun(P) - cutoff,
                   interval=c(1e-8, P_hat))$root
  
  upper <- uniroot(function(P) logLik_fun(P) - cutoff,
                   interval=c(P_hat, 0.5-1e-8))$root
  
  # Scaling
  skew_scaled      <- 100 - 100 * P_hat
  skew_scaled_low  <- 100 - 100 * upper
  skew_scaled_high <- 100 - 100 * lower
  
  # Return as data frame
  data.frame(
    P_hat = P_hat,
    CI_lower = lower,
    CI_upper = upper,
    skew_scaled = skew_scaled,
    skew_scaled_CI_lower = skew_scaled_low,
    skew_scaled_CI_upper = skew_scaled_high
  )
}

# Store this function. Its the same as below but does not output 95% CI when calculating the skewing from the phase sets.
compute_folded_skew <- function(df) {
  
  k <- df$H1_Xa + df$H2_Xi
  n <- df$H1_Xa + df$H2_Xa + df$H1_Xi + df$H2_Xi
  
  logLik_fun <- function(P) {
    
    if(P <= 0 || P >= 0.5) return(-Inf)
    
    log_term1 <- dbinom(k, n, P, log=TRUE)
    log_term2 <- dbinom(k, n, 1-P, log=TRUE)
    
    sum(log(exp(log_term1) + exp(log_term2)))
  }
  
  # MLE
  opt <- optimize(logLik_fun, interval=c(1e-8, 0.5-1e-8), maximum=TRUE)
  P_hat <- opt$maximum
  logLik_max <- opt$objective
  
  # Likelihood cutoff
  cutoff <- logLik_max - 0.5 * qchisq(0.95, df=1)
  
  # CI bounds
  #lower <- uniroot(function(P) logLik_fun(P) - cutoff,
  #                 interval=c(1e-8, P_hat))$root
  
  #upper <- uniroot(function(P) logLik_fun(P) - cutoff,
  #                 interval=c(P_hat, 0.5-1e-8))$root
  
  # Scaling
  skew_scaled      <- 100 - 100 * P_hat
  #skew_scaled_low  <- 100 - 100 * upper
  #skew_scaled_high <- 100 - 100 * lower
  
  # Return as data frame
  data.frame(
    P_hat = P_hat,
    #CI_lower = lower,
    #CI_upper = upper,
    skew_scaled = skew_scaled#,
    #skew_scaled_CI_lower = skew_scaled_low,
    #skew_scaled_CI_upper = skew_scaled_high
  )
}

samplez <- 
  c("TRiXi","xxh034","xxh006","xxh014","xxh020","xxh027","xxh176")

fin_df <- data.frame()
fin_skew_df <- data.frame()

for (i in samplez){
  
  # Access the argument(s)
  lib <- i
  bam <- paste0("data/revision_new_experiment/bams/",i,"/",i,"_haplotag_chrX.bam")
  BED <- read_tsv("anno/CGIs_hg38_chrX.bed", col_names = c("chr", "start", "end"))
  BED <- BED[BED$chr == "chrX",c("chr",     "start",     "end")]
  ncpus <- 6
  
  haplotyped_reads <- read_tsv(paste0("data/revision_new_experiment/bams/",i,"/",i,".haplotyped_reads.tsv"), col_names = c("read_name", "HP", "PS"))
  
  haplotyped_reads$HP <- gsub(haplotyped_reads$HP, pattern = "HP:i:", replacement = "")
  haplotyped_reads$PS <- gsub(haplotyped_reads$PS, pattern = "PS:i:", replacement = "")
  
  #create the ModBamResult object
  mbr <- ModBamResult(
    methy = ModBamFiles(
      samples = lib,
      paths = bam
    ),
    samples = data.frame(
      sample = lib,
      group = 1
    )
  )
  
  clustered_reads <- apply_cluster_reads_parallel(mbr, BED, min_pts = 3, num_cores = ncpus)
  
  write_tsv(clustered_reads, paste0(lib, "_excl_unclustered_reads_CGIX_clustered_reads.tsv.gz"))
  
  haplotyped_reads <- haplotyped_reads[haplotyped_reads$read_name %in% clustered_reads$read_name,]
  
  skew <- calculate_skew_by_block(clustered_reads, haplotyped_reads)
  
  skew$H1_Xa_skew_scaled <- 100 - 100 * skew$H1_Xa_skew
  
  write_tsv(skew, paste0(lib,"_excl_unclustered_reads_CGIX_skew.tsv.gz"))
  
  # compute folded skew
  P_folded <- compute_folded_skew(skew)
  
  temp_df <- data.frame(P_folded, i)
  
  fin_df <- rbind(fin_df, temp_df)
  
  skew_df <- skew
  
  skew_df$sample <- i
  
  fin_skew_df <- rbind(fin_skew_df, skew_df)
  
}

# export ASM methylation df 
write.table(fin_df, "Supplementary_tables/ASM_skewing_table.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(fin_skew_df, "Supplementary_tables/Supplementary_Table4.tsv", quote = F, sep = "\t", row.names = F, col.names = T)


