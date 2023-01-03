#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: identify the sensitity of long-read data 
## (i.e. how many samples detect unique isoforms, what proportion of isoforms are detected across all samples)
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## no_of_isoforms_sample
##   
## ---------- Notes -----------------
## 
## Pre-requisite: requires the SQANTI classification file to have full-length read counts ("FL" column) for each sample
##


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- no_of_isoforms_sample

# Aim: plot the number of samples with detected expression of isoforms
# Prerequisite: class.files from SQANTI includes abundance (counts) per sample
# Input:
  # classification file generated from SQANTI
# Output:
  # p1 = the proportion of isoforms detected across samples
  # p2 = the number of samples against median read count

no_of_isoforms_sample <- function(class){
  
  # median_FL: median full-length read counts for each isoform (row)
  # num_samples: the number of samples where 0 read counts for isoform
  # across each row (i.e isoform, count the number of occurences where reads are != 0)
  dat <- class %>% dplyr::select(starts_with("FL.")) %>% 
    mutate(median_FL = apply(.,1, function(x) median(x)), 
           num_samples = apply(.,1, function(x) length(x[which(x != "0")]))) 
  
  table(dat$num_samples)
  table(dat$num_samples)/sum(table(dat$num_samples))
  
  p1 <- ggplot(dat, aes(x = as.factor(num_samples))) + geom_bar(aes(y = (..count..)/sum(..count..))) + 
    scale_y_continuous(labels = perc_lab)  + mytheme + labs(x = "Number of Samples", y = "Isoforms (%)")
  
  p2 <- ggplot(dat, aes(x = as.factor(num_samples), y = log(median_FL))) + geom_boxplot() + 
    mytheme + labs(x = "Number of Samples", y = "FL Read Count(Log10)")
  
  return(list(p1,p2))
}