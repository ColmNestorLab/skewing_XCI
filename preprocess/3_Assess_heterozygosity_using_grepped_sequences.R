# Load necessary library
library(dplyr)
library(ggplot2)
library(data.table)
library(biomaRt)
library(tidyverse)
library(karyoploteR)

# Before this, the two scripts '1_bash_preprocess_repeats.sh' and '2_Get_flanking_sequences_of_repeats_for_grepping.R' needs to be run.

# good to have.
source("plot_parameters_TRiX.R")

# repeat data
trip_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/trip_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat"))

quad_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quad_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat"))

quin_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/quin_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat"))

auto_data <- fread("/home/bjogy93/Desktop/TREX_proj/TREX_bash/repeat_output/autosomal_repeats_sequences_to_grep.tsv", 
                   col.names = c("contig", "repeat_start", "repeat_start_20_before_100_after", "repeat_type", "sequence_directly_after_repeat", "closest_gene"))


# Connect to the Ensembl database for human genome (hg38)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Retrieve gene information for the X chromosome
genes <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
               filters = 'chromosome_name', values = 'X', mart = ensembl)


# Function to find the closest gene for a given position
find_closest_gene <- function(pos) {
  # Calculate the distance to each gene (start or end position)
  genes <- genes %>%
    mutate(distance = pmin(abs(start_position - pos), abs(end_position - pos))) %>%
    arrange(distance)
  
  # If the position is within a gene, return that gene; otherwise, return the closest one
  within_gene <- genes %>% filter(start_position <= pos & end_position >= pos)
  
  if (nrow(within_gene) > 0) {
    return(within_gene$hgnc_symbol[1])
  } else {
    return(genes$hgnc_symbol[1])
  }
}

# Apply the function to each position in the repeat data.
trip_data <- trip_data %>%
  rowwise() %>%
  mutate(closest_gene = find_closest_gene(repeat_start))

# add gene name. The repeats are not necessarily in the actual gene, here I just add the closest gene.
trip_data[trip_data$repeat_start %in% 37737518,]$closest_gene <- "XK"
trip_data[trip_data$repeat_start %in% 38562575,]$closest_gene <- "TSPAN7"
trip_data[trip_data$repeat_start %in% 48631673,]$closest_gene <- "WDR13"
trip_data[trip_data$repeat_start %in% 69082485,]$closest_gene <- "PJA1"
trip_data[trip_data$repeat_start %in% 69083812,]$closest_gene <- "PJA1_2"
trip_data[trip_data$repeat_start %in% 123767927,]$closest_gene <- "THOC"

quad_data <- quad_data %>%
  rowwise() %>%
  mutate(closest_gene = find_closest_gene(repeat_start))

# add gene name. The repeats are not necessarily in the actual gene, here I just add the closest gene.
quad_data[quad_data$repeat_start %in% 833522,]$closest_gene <- "SHOX"
quad_data[quad_data$repeat_start %in% 1109474,]$closest_gene <- "CRLF2"
quad_data[quad_data$repeat_start %in% 68903924,]$closest_gene <- "EFNB1"
quad_data[quad_data$repeat_start %in% 127008981,]$closest_gene <- "PRR32"
quad_data[quad_data$repeat_start %in% 1658563,]$closest_gene <- "ASMT"
quad_data[quad_data$repeat_start %in% 4551340,]$closest_gene <- "FAM239B"
quad_data[quad_data$repeat_start %in% 39613271,]$closest_gene <- "MIR3937"
quad_data[quad_data$repeat_start %in% 48244717,]$closest_gene <- "SSX1"
quad_data[quad_data$repeat_start %in% 88417303,]$closest_gene <- "CPXCR1"
quad_data[quad_data$repeat_start %in% 137545369,]$closest_gene <- "LINC02931"
quad_data[quad_data$repeat_start %in% 137554666,]$closest_gene <- "LINC02931_2"
quad_data[quad_data$repeat_start %in% 137568340,]$closest_gene <- "ZIC3"

# Added closest gene manually to the autosomal data.


### TRIPLET ###
# Specify the folder containing the text files
folder_path_trip <- "/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/trip"

# List all text files in the folder
file_list_trip <- list.files(path = folder_path_trip, pattern = "\\.txt$", full.names = TRUE)

# Get file information, including size
file_info <- file.info(file_list_trip)

# Filter files to exclude those with size 0
file_list_trip <- file_list_trip[file_info$size > 0]

# Function to read a file and add the file name as a column
read_and_label_trip <- function(file) {
  # Read the text file with fill = TRUE to handle rows with different numbers of elements
  data <- fread(file, header = F, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  data$file_namez <- basename(file)
  data <- data %>% tidyr::separate(file_namez, into = c("name", "repeat_grepped", NA), sep = "_", remove = T)
  data_smash <- merge(data, trip_data, by = c("repeat_grepped"), by.y = "repeat_start_20_before_100_after")
  
  data_smash <- dplyr::select(data_smash, contig, repeat_start, closest_gene, repeat_type, repeat_grepped, sequence_directly_after_repeat, V1, name)
  
  # Add a column with the file name (without path)
  colnames(data_smash)[colnames(data_smash) == "V1"] <- "grep_result"
  
  return(data_smash)
}

# Apply the function to all files and combine the results into one data frame
df_trip <- bind_rows(lapply(file_list_trip, read_and_label_trip))

# Shorten the flanking sequence (to allow for more hits).
df_trip$bases_after_repeat_5 <- substr(df_trip$sequence_directly_after_repeat, 1, 5)

df_trip$match <- mapply(function(long, short) grepl(short, long), df_trip$grep_result, df_trip$bases_after_repeat_5)


# Function to find the longest unbroken repeat sequence followed by flanking sequence
count_longest_repeats <- function(df) {
  df$num_repeats <- sapply(1:nrow(df), function(i) {
    sequence <- df$grep_result[i]
    flanking <- df$bases_after_repeat_5[i]
    repeat_seq <- df$repeat_type[i]
    
    # Regular expression to match the longest unbroken repeat directly before the flanking sequence
    pattern <- paste0("(", repeat_seq, ")+", flanking)
    
    # Find all matches in the sequence
    matches <- gregexpr(pattern, sequence, perl = TRUE)[[1]]
    
    if (matches[1] != -1) {
      # Extract matching sequences
      matched_sequences <- regmatches(sequence, gregexpr(pattern, sequence, perl = TRUE))[[1]]
      
      # Find the longest unbroken repeat count for each match
      max_repeats <- max(sapply(matched_sequences, function(match) {
        repeat_only <- sub(flanking, "", match) # Remove the flanking part
        nchar(repeat_only) / nchar(repeat_seq) # Count repeats
      }))
      return(max_repeats)
    } else {
      return(0) # No match found
    }
  })
  return(df)
}


# Apply the function
result_df_trip <- count_longest_repeats(df_trip)

# Print the result
print(result_df_trip)

trip_repeats_distrubition <- ggplot(result_df_trip[result_df_trip$num_repeats > 6,], aes(x=num_repeats)) + geom_histogram(binwidth = 1) + facet_wrap(~closest_gene, scales = "free")










### quadLET ###
# Specify the folder containing the text files
folder_path_quad <- "/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/quad"

# List all text files in the folder
file_list_quad <- list.files(path = folder_path_quad, pattern = "\\.txt$", full.names = TRUE)

# Get file information, including size
file_info <- file.info(file_list_quad)

# Filter files to exclude those with size 0
file_list_quad <- file_list_quad[file_info$size > 0]

# Function to read a file and add the file name as a column
read_and_label_quad <- function(file) {
  # Read the text file with fill = TRUE to handle rows with different numbers of elements
  data <- fread(file, header = F, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  data$file_namez <- basename(file)
  data <- data %>% tidyr::separate(file_namez, into = c("name", "repeat_grepped", NA), sep = "_", remove = T)
  data_smash <- merge(data, quad_data, by = c("repeat_grepped"), by.y = "repeat_start_20_before_100_after")
  
  data_smash <- dplyr::select(data_smash, contig, repeat_start, closest_gene, repeat_type, repeat_grepped, sequence_directly_after_repeat, V1, name)
  
  # Add a column with the file name (without path)
  colnames(data_smash)[colnames(data_smash) == "V1"] <- "grep_result"
  
  return(data_smash)
}

# Apply the function to all files and combine the results into one data frame
df_quad <- bind_rows(lapply(file_list_quad, read_and_label_quad))

# Shorten the flanking sequence (to allow for more hits).
df_quad$bases_after_repeat_5 <- substr(df_quad$sequence_directly_after_repeat, 1, 8)

df_quad$match <- mapply(function(long, short) grepl(short, long), df_quad$grep_result, df_quad$bases_after_repeat_5)


# Function to find the longest unbroken repeat sequence followed by flanking sequence
count_longest_repeats <- function(df) {
  df$num_repeats <- sapply(1:nrow(df), function(i) {
    sequence <- df$grep_result[i]
    flanking <- df$bases_after_repeat_5[i]
    repeat_seq <- df$repeat_type[i]
    
    # Regular expression to match the longest unbroken repeat directly before the flanking sequence
    pattern <- paste0("(", repeat_seq, ")+", flanking)
    
    # Find all matches in the sequence
    matches <- gregexpr(pattern, sequence, perl = TRUE)[[1]]
    
    if (matches[1] != -1) {
      # Extract matching sequences
      matched_sequences <- regmatches(sequence, gregexpr(pattern, sequence, perl = TRUE))[[1]]
      
      # Find the longest unbroken repeat count for each match
      max_repeats <- max(sapply(matched_sequences, function(match) {
        repeat_only <- sub(flanking, "", match) # Remove the flanking part
        nchar(repeat_only) / nchar(repeat_seq) # Count repeats
      }))
      return(max_repeats)
    } else {
      return(0) # No match found
    }
  })
  return(df)
}

# Apply the function
result_df_quad <- count_longest_repeats(df_quad)

# Print the result
print(result_df_quad)

quad_repeats_distrubition <- ggplot(result_df_quad[result_df_quad$num_repeats > 6,], aes(x=num_repeats)) + geom_histogram(binwidth = 1) + facet_wrap(~closest_gene, scales = "free")







### autoLET ###
# Specify the folder containing the text files
folder_path_auto <- "/home/bjogy93/Desktop/TREX_proj/TREX_bash/grep_results/autosomal/"

# List all text files in the folder
file_list_auto <- list.files(path = folder_path_auto, pattern = "\\.txt$", full.names = TRUE)

# Get file information, including size
file_info <- file.info(file_list_auto)

# Filter files to exclude those with size 0
file_list_auto <- file_list_auto[file_info$size > 0]

# Function to read a file and add the file name as a column
read_and_label_auto <- function(file) {
  # Read the text file with fill = TRUE to handle rows with different numbers of elements
  data <- fread(file, header = F, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  data$file_namez <- basename(file)
  data <- data %>% tidyr::separate(file_namez, into = c("name", "repeat_grepped", NA), sep = "_", remove = T)
  data_smash <- merge(data, auto_data, by = c("repeat_grepped"), by.y = "repeat_start_20_before_100_after")
  
  data_smash <- dplyr::select(data_smash, contig, repeat_start, closest_gene, repeat_type, repeat_grepped, sequence_directly_after_repeat, V1, name)
  
  # Add a column with the file name (without path)
  colnames(data_smash)[colnames(data_smash) == "V1"] <- "grep_result"
  
  return(data_smash)
}

# Apply the function to all files and combine the results into one data frame
df_auto <- bind_rows(lapply(file_list_auto, read_and_label_auto))

# Shorten the flanking sequence (to allow for more hits).
df_auto$bases_after_repeat_5 <- substr(df_auto$sequence_directly_after_repeat, 1, 5)

df_auto$match <- mapply(function(long, short) grepl(short, long), df_auto$grep_result, df_auto$bases_after_repeat_5)


# Function to find the longest unbroken repeat sequence followed by flanking sequence
count_longest_repeats <- function(df) {
  df$num_repeats <- sapply(1:nrow(df), function(i) {
    sequence <- df$grep_result[i]
    flanking <- df$bases_after_repeat_5[i]
    repeat_seq <- df$repeat_type[i]
    
    # Regular expression to match the longest unbroken repeat directly before the flanking sequence
    pattern <- paste0("(", repeat_seq, ")+", flanking)
    
    # Find all matches in the sequence
    matches <- gregexpr(pattern, sequence, perl = TRUE)[[1]]
    
    if (matches[1] != -1) {
      # Extract matching sequences
      matched_sequences <- regmatches(sequence, gregexpr(pattern, sequence, perl = TRUE))[[1]]
      
      # Find the longest unbroken repeat count for each match
      max_repeats <- max(sapply(matched_sequences, function(match) {
        repeat_only <- sub(flanking, "", match) # Remove the flanking part
        nchar(repeat_only) / nchar(repeat_seq) # Count repeats
      }))
      return(max_repeats)
    } else {
      return(0) # No match found
    }
  })
  return(df)
}


# Apply the function
result_df_auto <- count_longest_repeats(df_auto)

# Print the result
print(result_df_auto)

auto_repeats_distrubition <- ggplot(result_df_auto[result_df_auto$num_repeats > 6,], aes(x=num_repeats)) + geom_histogram(binwidth = 1) + facet_wrap(~closest_gene, scales = "free")

# outputs repeat length distributions for all repeats
#pdf("/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/01_plots/Fig_1B_repeat_length_distribution_with_autosomal.pdf")
#trip_repeats_distrubition
#quad_repeats_distrubition
#auto_repeats_distrubition
#dev.off()

# export tables
#write.table(result_df_trip, "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_trip.tsv", quote = F, row.names = F, sep = "\t")
#write.table(result_df_quad, "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_quad.tsv", quote = F, row.names = F, sep = "\t")
#write.table(result_df_auto, "/home/bjogy93/Desktop/TREX_proj/TREX_R_remotes/data/result_df_auto.tsv", quote = F, row.names = F, sep = "\t")