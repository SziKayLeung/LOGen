#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot the distribution of Iso-Seq raw read lengths; CCS level, polished level etc 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## read_rarefaction_files
## rarefaction_distribution
##   
## ---------- Notes -----------------
##
## Pre-requisite: generate files using E.Tseng's cupcake script: get_seq_stats.py

## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))


## ------------------- tab_raw_isoseq_length

# Aim: tabulate the Iso-Seq raw read data 
# Input:
  # lengths_input_dir = path of directory with files generated from E.Tseng scripts for lengths 
  # suffix_name = str: common string for defining the input file 
  # type = str: <CCS> <Clustered> 
    # defines the type of Iso-Seq raw read data (for naming purposes)
# Output:
  # 1 dataframe: <sample> <length> <type>

tab_raw_isoseq_length <- function(lengths_input_dir,suffix_name, type){
  
    # read in list of files ending with the input suffix name
    filenames <- as.list(list.files(path = lengths_input_dir, pattern = paste0(suffix_name,"$"), full.names = TRUE))
    files <- lapply(filenames, read.table)
    names(files) <- list.files(path = lengths_input_dir, pattern = paste0(suffix_name,"$"))
    
    # differentiate files prior to merging by creating column with the name of the file
    for(i in 1:length(files)){files[[i]]$File <- names(files)[[i]]}
    
    # bind all files
    length_df <- bind_rows(files) %>% mutate(Sample = word(File, c(1), sep = fixed(".")), Type = type)
   
    # create column for genotype depending on the sample column
    length_df$Genotype <- sapply(length_df$Sample, function(x) classify_genotype(x, case_samples, control_samples))
    return(length_df)

}


## ------------------- plot_raw_isoseq_length

# Aim: plot distribution of the lengths of the Iso-Seq raw read data 
# Input:
  # length_df = output from tab_raw_isoseq_length
  # type = str: <CCS> <Clustered>;  same type value as tab_raw_isoseq_length
# Output: density plot of the read length between cases and control 

plot_raw_isoseq_length <- function(length_df, type){
  
  p = dplyr::filter(length_df, !grepl('Merged', Sample)) %>% 
    filter(Type == type) %>% 
    ggplot(., aes(V2, fill = Genotype)) +
    geom_density(alpha = 0.2) +
    labs(x = paste0(type, " Read Length (kb)"), y = "Density") +
    theme_bw() + mytheme +
    theme(legend.title = element_blank(), legend.position = c(0.9, 0.9)) +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    scale_fill_manual(values = c(label_colour(label_name("case")),label_colour(label_name("control"))))
  
  return(p)
  
}