#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot number of ERCC molecules detected in dataset and correlation with known amounts added 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## number_ERCC_detected
##   
## ---------- Notes -----------------
## 
## Pre-requisite: 
## 1. Ran alignment of dataset to ERCC reference, followed by SQANTI
## 2. Dataframe of information for each ERCC molecule: <ERCC.ID> <conc> <amoles> <amount_of_ERCC_added>
    # <amount_of_ERCC_added> = known amount of ERCC from calculations
## 3. function: draw_density.R from aesthetics_basics_plots


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- number_ERCC_detected

# Aim: plot ERCC-related isoforms in dataset
# Input:
  # class.files = df: classification file of associated ERCC molecules generated from SQANTI (following ERCC minimap2 alignment)
  # ercc_calculation = df: <ERCC.ID> <conc> <amoles> <amount_of_ERCC>; original concentration of each ERCC 
# Output: 
  # p1: bar-plot of number of isoforms per ERCC molecule 
  # p2: scatter-plot of amount of reads vs amount of each ERCC molecule
  # p3: correlation of amount of reads vs amount of each ERCC molecule
  
number_ERCC_detected <- function(class.files, ercc_calculation){
  
  cat("Total unique ERCCs:", length(unique(class.files$chrom)), paste0("(", round(length(unique(class.files$chrom))/92 * 100,2), "%)"))
  
  # redundant 
  redundant <- class.files[,c("chrom")] %>% table() %>% reshape2::melt(.)
  colnames(redundant) <- c("ERCC","num_isoforms")
  
  p1 <- redundant %>% 
    group_by(num_isoforms) %>% tally() %>% ggplot(., aes(x = as.factor(num_isoforms), y = n)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "ERCC", y = "Number of Isoforms")
  
  # isoform vs concentration
  isoform_conc <- merge(redundant, ercc_calculation, by.x = "ERCC", by.y = "ERCC.ID", all = TRUE) %>%
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    replace_na(list(num_isoforms = 0)) %>% 
    mutate(num_isoforms = as.factor(num_isoforms))
  
  p2 <- ggplot(isoform_conc, aes(x = num_isoforms, y = log2_amount_of_ERCC, colour = num_isoforms)) + 
    geom_jitter(width = 0.2) + 
    theme(legend.position = "none") + 
    mytheme + labs(x = "Number of Isoforms", y = "Amount of ERCC (Log2)") + theme(legend.position = "none")
  
  # correlation 
  ERCC_corr <- merge(class.files,ercc_calculation, by.x = "chrom", by.y = "ERCC.ID") %>% 
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    mutate(log2_FL_reads = log2(FL)) 
  
  p3 <- density_plot(ERCC_corr,"log2_amount_of_ERCC","log2_FL_reads", "Amount of ERCC (Log2)", "Number of FL Reads (Log2)","") 
  
  return(list(p1,p2,p3))
}