#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot rarefaction curves for Iso-Seq data 
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
## Pre-requisite: generate files using E.Tseng's cupcake script: make_file_for_subsampling_from_collapsed.py

## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- read_rarefaction_files

# Aim: Input rarefaction files generated from E.Tseng rarefaction scripts 
# Input:
  # input_rarefaction_dir = directory containing files generated from E.Tseng scripts (i.e <sample>_rarefaction.by_refgene.min_fl_2.txt)
  # name = prefix of files 
# Output:
  # List: <genes> <isoforms> <isoforms_cate> = lengths at the different levels

read_rarefaction_files <- function(input_rarefaction_dir, name){
  
  input_files <- list(
    genes = paste0(input_rarefaction_dir,name,'.rarefaction.by_refgene.min_fl_2.txt'),
    isoforms = paste0(input_rarefaction_dir,name,'.rarefaction.by_refisoform.min_fl_2.txt'),
    isoforms_cate = paste0(input_rarefaction_dir,name,'.rarefaction.by_refisoform.min_fl_2.by_category.txt')
  )
  
  input_files <- lapply(input_files, function(x) read.table(x, sep=' ',header=T,skip=1))
  
  # datawrangle
  input_files$genes <- input_files$genes %>% mutate(type = "Genes")
  input_files$isoforms <- input_files$isoforms %>% mutate(type = "Isoforms")
  input_files$isoforms_cate$category <- factor(input_files$isoforms_cate$category, 
                                                       levels = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", 
                                                                  "novel_not_in_catalog", "fusion", "antisense", "genic", "intergenic"),
                                                       labels = c("FSM", "ISM", "NIC", "NNC", "Fusion", "Antisense", "Genic", "Intergenic"))
  
  for(j in 1:nrow(input_files$isoforms_cate)){
    input_files$isoforms_cate$type[j] <- if(input_files$isoforms_cate$category[j] %in% c("FSM", "ISM")){
      "Annotated" 
    } else if (input_files$isoforms_cate$category[j] %in% c("NIC", "NNC")){
      "Novel"
    } else {
      "Others"
    }
  }
  
  return(input_files)
}


## ------------------- rarefaction_distribution

# Aim: plot rarefaction 
# Input:
  # rarefaction files generated from E.Tseng scripts (output from read_rarefaction_files) 
# Output:
  # p1 <- isoform and gene level 
  # p2 <- isoform per category for individual datasets

rarefaction_distribution <- function(all_rarefaction){
  
  # Merge input files pertainng to genes and isoforms 
  all_rarefaction_levels <- bind_rows(all_rarefaction[["genes"]], all_rarefaction[["isoforms"]])
  
  ## Plots
  p1 <- all_rarefaction_levels %>% 
    ggplot(., aes(x = size, y = mean, linetype = type)) + 
    geom_line(size = 1.5) + 
    labs(x ="Number of Subsampled Reads (K)", y = "Number of Genes/Isoforms (K)") + 
    theme_bw() + mytheme + 
    scale_y_continuous(labels = ks) + scale_x_continuous(labels = ks) + 
    #scale_color_manual(values=c(label_colour("mouse")), name = "") + 
    scale_linetype_manual(values=c("dotted","solid")) +
    theme(legend.position = c(0.8, 0.6), legend.spacing.y = unit(-0.1, "cm"),legend.title = element_blank())
  
  
  p2 <- ggplot(all_rarefaction[["isoforms_cate"]], aes(x=size, y=mean, color=category)) + geom_line(aes(linetype = type), size = 1.5) + 
    labs(x = "Number of Subsampled Reads (K)", 
         y = "Number of Isoforms (K)") +
    mytheme + 
    scale_y_continuous(labels = ks, limits = c(0, 20000)) + 
    scale_x_continuous(labels = ks, expand = c(0.15, 0)) + 
    theme(legend.position = "none") +
    scale_linetype_manual(values=c("solid", "longdash","dotted")) +
    geom_dl(aes(label = category),  method = list(dl.trans(x = x + 0.5), "last.bumpup", cex = 1.3, hjust = .1))
  
  return(list(p1,p2))
}
