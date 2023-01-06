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
## read_target_probes
##   
## ---------- Notes -----------------
##
## Pre-requisite: 
## 1. generated and read in results from E.Tseng scripts for calculating on-probe target rate


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))


## ------------------- read_target_probes

# Aim: read in the files generated from E.Tseng scripts for downstream plotting 
# Input:
  # dir = str: path of directory containing files 
  # filename = str: pattern to grep for files 
# Output:
  # list of read in files 

read_target_probes <- function(dir, filename){
  
  # list files in directory with pattern and read
  Probes_input <- list.files(path = dir, pattern = filename, full.names = T)
  print(Probes_input)
  Probes_files <- lapply(Probes_input, function(x) read.table(x, header=T, as.is=T, sep="\t"))
  names(Probes_files) <- list.files(path = dir, pattern = filename)
  
  return(Probes_files)
  
}


## ------------------- on_target_plot

# Aim: plot the on-target rate based on the number of bases which overlap with probes
# Input:
  # Probes_files = list: read in files of fasta.sam.probe_hit.txt
  # targetedpheno = df: phenotype of samples <Sample/BarcodedSample> <Phenotype> <Batch> <Barcode> 
  # type = str: <notbatched> <batched>, used to determine whether to include samples in batches
# Note: 
  # function only relevant for targeted dataset with multiple samples that are barcoded 
  # ensure that Probes_files have a "BatchX.BCX" in the file list name if using the "Batch" type
# Output: box-plot of the on-target rate per batch and colour coded by phenotype

on_target_plot <- function(Probes_files, targetedpheno, type){
  
  # merge on the output across samples
  Probes <- merge(ldply(Probes_files, function(x) nrow(x)),
                  ldply(Probes_files, function(x) length(which(x$num_base_overlap != "0"))),by = ".id") 
  
  if(type == "notbatched"){
    
    Probes <- Probes %>%
      `colnames<-`(c("file", "Total_mapped_reads", "reads_probe_hit")) %>% 
      mutate(perc = reads_probe_hit/Total_mapped_reads * 100) %>%
      mutate(sample = word(.$file, c(1), sep = fixed("."))) %>% 
      full_join(., targetedpheno, by = c("sample" = "Sample"))

  }else if (type == "batched"){
    
    Probes <- Probes %>%
      `colnames<-`(c("file", "Total_mapped_reads", "reads_probe_hit")) %>% 
      mutate(perc = reads_probe_hit/Total_mapped_reads * 100,
             Batch = as.factor(word(.$file, c(1), sep = fixed("."))),
             Sample = as.factor(word(.$file, c(2), sep = fixed("."))),
             BarcodedSample = as.factor(paste0(Batch,Sample))) %>%
      left_join(., targetedpheno %>% select(-Batch,-Sample,-sample), by = "BarcodedSample")
    
  }else{
    print("error: 3rd argument should be <notbatched> <batched>")
  }
   
  # plot
  p1 <- ggplot(Probes, aes(x = as.factor(Batch), y = perc, fill = as.factor(Phenotype))) + geom_boxplot() +
    geom_point(aes(colour = as.factor(Phenotype)), position = position_jitterdodge(), size = 3) + 
    mytheme + labs(y = "On-Target Rate (%)", x = "Batch") + 
    scale_fill_manual(values = c(alpha(label_colour("TG"),0.4),alpha(label_colour("WT"),0.4)), name = "Genotype") + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), guide=FALSE) 
  
  for(i in 1:3){cat("Mean on target rate in Batch",i,":", 
                    mean(Probes[Probes$Batch == i,"perc"]),"\n")}
  return(p1)
}