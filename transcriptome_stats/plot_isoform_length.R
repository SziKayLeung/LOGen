#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot distribution of isoform length
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## multi_isoform_length
##   
## ---------- Notes -----------------
## 
## Only of multi-exonic isoforms (mono-exonic isoforms are removed)


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- multi_isoform_length

# Aim: plot the length of multi-exonic isoforms 
# Input: 
  # class.files = df: classification file generated from SQANTI

multi_isoform_length <- function(class.files){
  # remove isoforms that are mono-exonic
  class.files <- class.files %>% filter(subcategory != "mono-exon")
  
  p <- ggplot(class.files, aes(x = length)) + geom_histogram(bins = 15, fill="gray", col="black") + 
    labs(x = "Transcript Length (kb)", y = "Number of Isoforms (K)") + mytheme +
    scale_x_continuous(labels = ks) + 
    scale_y_continuous(labels = ks) 
  
  # how many isoforms have length between 2 - 4kb
  two <- class.files[which(class.files$length >= 2000 & class.files$length <= 4000),] %>% nrow()
  print(paste0("Number of isoforms 2-4kb:", two, "(",round(two/nrow(class.files),2) *100,"%)"))
  
  return(p)
}