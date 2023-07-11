#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Compare Iso-Seq whole vs targeted transcriptome  
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## whole_vs_targeted_exp
## whole_vs_targeted_plots
##   
## ---------- Notes -----------------
##

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "compare_datasets/dataset_identifer.R"))
suppressMessages(library("ggplot2"))
suppressMessages(library("wesanderson"))

## ------------------- whole_vs_targeted_exp

# Aim: plot expression of dataset 
# Input:
  # class.files: df = classification file from SQANTI 
  # dataset: str = <"targeted", "whole>, used to determine colour for dataset
# Output: box-plot of expression 

whole_vs_targeted_exp <- function(class.files, dataset){
  
  # colour classifier
  if(dataset == "Targeted"){colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("targeted")))
  }else{colour <-  scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("whole")))}
  
  # plot
  class.files$TargetGene = factor(class.files$TargetGene, levels=c("Target Genes","Not Target Genes"))
  p = ggplot(class.files, aes(x = Matching, y = log10(FL), fill = Matching)) + geom_boxplot() + facet_grid(~TargetGene) + mytheme + 
    labs(x = "", y = paste0("FL Read Counts \n",dataset," Transcriptome (Log10)")) + colour + theme(legend.position = "none")
  
  return(p)
}

## ------------------- whole_vs_targeted_plots

# Aim: generate plots comparing whole vs targeted transcriptome of the same subset samples
# Input:
  # class.files: read in SQANTI classification file
  # targetGene: vec: target genes 
  # wholeSamples: vec: column names of the whole transcriptome sample in the SQANTI classification file
  # targetedSamples: vec: column names of the targeted transcriptome sample in the SQANTI classification file
# Pre-requisite:
  # merged the demultiplexed reads from the whole and targeted transcriptome 
  # align using pbmm2
  # iso-seq3 collapse, sqanti3 QC and filter
# Output:
  # p1: bar_plot of the number of isoforms in both, targeted and whole dataset only across target genes
  # p2: bar-plot of the number of isoforms by structural category in both and unique to targeted dataset
  # matchedSumTargeted: list of the isoforms and sum FL read expression across targeted and whole transcriptome dataset

whole_vs_targeted_plots <- function(class.files, wholeSamples, targetedSamples, targetGene){
  
  # sum the number of FL reads across the targeted and whole samples
  matchedSum = merge(data.frame(class.files %>% select(all_of(targetedSamples)) %>% apply(., 1, sum)),
                     data.frame(class.files %>% select(all_of(wholeSamples)) %>% apply(., 1, sum)),
                     by = 0, all = T)
  colnames(matchedSum) = c("isoform","sumTargeted","sumWhole")
  
  # annotate the isoforms and subset by target genes
  matchedSum = merge(matchedSum, class.files [,c("isoform","associated_gene","structural_category")])
  matchedSumTargeted = data.frame(matchedSum %>% filter(associated_gene %in% targetGene))
  
  # create a column by determining if the isoform is detected in both, whole or targeted
  # if FL read >= 1; then considered detected
  matchedSumTargeted$dataset <- apply(matchedSumTargeted, 1, function(x) identify_dataset_by_counts (x[["sumTargeted"]], x[["sumWhole"]], "Targeted","Whole"))
  
  # plots
  totaln = matchedSumTargeted %>% group_by(associated_gene) %>% tally 
  p1 <- matchedSumTargeted %>% group_by(associated_gene, dataset) %>% tally %>% 
    full_join(., totaln, by = "associated_gene") %>%
    ggplot(., aes(x = reorder(associated_gene,-n.y), y = n.x, fill = dataset)) + geom_bar(stat = "identity") +
    mytheme + labs(x = "Target Genes", y = "Number of isoforms") +
    scale_fill_manual(name = "", values = c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling2")[1],wes_palette("Darjeeling1")[2])) + 
    theme(legend.position = "top") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p2 <- matchedSumTargeted %>% filter(dataset != "Both") %>%
    group_by(structural_category, dataset) %>% tally() %>% 
    ggplot(., aes(x = structural_category, y = n, fill = dataset)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural Category", y = "Number of isoforms") +
    scale_fill_manual(name = "", values = c(wes_palette("Darjeeling2")[1],wes_palette("Darjeeling1")[2])) + 
    theme(legend.position = "top")
  
  return(list(p1,p2,matchedSumTargeted))
}

