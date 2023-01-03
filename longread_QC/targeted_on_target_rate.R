#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Determine on-target rate of long-read targeted experiments 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## on_target_plot
##   
## ---------- Notes -----------------
##
## Pre-requisite: 
## 1. generated and read in reults from E.Tseng scripts for calculating on-probe target rate


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))


## ------------------- on_target_plot

# Aim: plot the on-target rate based on the number of bases which overlap with probes
# Input:
  # Probes_files = list: read in files of fasta.sam.probe_hit.txt
  # targetedpheno = df: phenotype of samples <Sample> <Phenotype> <Batch> <Barcode> 
# Note: function only relevant for targeted dataset with multiple samples that are barocoded and batched
# Output: box-plot of the on-target rate per batch and colour coded by phenotype

on_target_plot <- function(Probes_files, targetedpheno){
  
  # merge on the output across samples
  Probes <- merge(ldply(Probes_files, function(x) nrow(x)),
                  ldply(Probes_files, function(x) length(which(x$num_base_overlap != "0"))),by = ".id") %>%
    `colnames<-`(c("file", "Total_mapped_reads", "reads_probe_hit")) %>% 
    # generate columns for downstream plotting
    mutate(perc = reads_probe_hit/Total_mapped_reads * 100) %>%
    mutate(sample = word(.$file, c(1), sep = fixed("."))) %>% 
    full_join(., targetedpheno, by = c("sample" = "Sample"))
  
  # plot
  p1<- ggplot(Probes, aes(x = as.factor(Batch), y = perc, fill = as.factor(Phenotype))) + geom_boxplot() +
    geom_point(aes(colour = as.factor(Phenotype)), position = position_jitterdodge(), size = 3) + 
    mytheme + labs(y = "On-Target Rate (%)", x = "Batch") + 
    scale_fill_manual(values = c(alpha(label_colour("TG"),0.4),alpha(label_colour("WT"),0.4)), name = "Genotype") + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), guide="none") 
  
  for(i in 1:3){cat("Mean on target rate in Batch",i,":", 
                    mean(Probes[Probes$Batch == i,"perc"]),"\n")}
  return(p1)
}