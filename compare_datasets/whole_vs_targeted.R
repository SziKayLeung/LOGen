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
## Pre-requisite: 
## 1. generated and read in reults from E.Tseng scripts for calculating on-probe target rate


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
  # cuff_tmap: gffcompare result from comparison of whole vs targeted transcriptome
  # TargetGene: vec: target genes 
  # whole.class.files: df = classification file of whole transcriptome dataset of X samples
  # subsettargeted.class.files: df = classification file of subsetted targeted dataset of X matched samples
# Pre-requisite:
  # run gffcompare on linux 
# Output:
  # p1: bar_plot of the number of isoforms in both, targeted and whole dataset only across target genes
  # p2: bar-plot of the number of isoforms by structural category in both and unique to targeted dataset
  # p3: expression of isoforms unique to whole transcriptome data and commonly identified
  # p4: expression of isoforms unique to targeted transcriptome data and commonly identified

whole_vs_targeted_plots <- function(cuff_tmap, TargetGene, whole.class.files,subsettargeted.class.files){
  
  # identify commonly matched isoforms from gffcompare output
  cuff_tmap_exact = cuff_tmap[cuff_tmap$class_code == "=",]
  whole.class.files = whole.class.files %>% mutate(Matching = ifelse(isoform %in% cuff_tmap_exact$ref_id,"Both","Whole"))
  subsettargeted.class.files = subsettargeted.class.files %>% mutate(Matching = ifelse(isoform %in% cuff_tmap_exact$qry_id,"Both","Targeted")) 
  
  # plots
  cols = c("isoform","associated_gene","Matching","structural_category")
  p1 = rbind(whole.class.files[whole.class.files$Matching != "Both",cols],subsettargeted.class.files[,cols]) %>% 
    filter(associated_gene %in% TargetGene) %>%
    group_by(associated_gene, Matching) %>% tally() %>%
    ggplot(., aes(x = reorder(associated_gene, -n), fill = Matching, y = n)) + geom_bar(stat = "identity") +
    mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(name = "", values = c(label_colour("whole"),label_colour("whole+targeted"),label_colour("targeted")))
  
  p2 = subsettargeted.class.files[subsettargeted.class.files$associated_gene %in% TargetGene,] %>%
    group_by(structural_category, Matching) %>% tally() %>% 
    ggplot(., aes(x = structural_category, y = n, fill = Matching)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural Category", y = "Number of Isoforms \n Targeted Transcriptome") +
    scale_fill_manual(name = "", values = c(label_colour("whole+targeted"),label_colour("targeted"))) + 
    theme(legend.position = "top")
  
  print(rbind(whole.class.files[whole.class.files$Matching != "Both",cols],subsettargeted.class.files[,cols]) %>% 
          filter(associated_gene %in% TargetGene) %>%
          group_by(Matching) %>% tally())
  
  p3 = whole_vs_targeted_exp(whole.class.files,"Whole")
  p4 = whole_vs_targeted_exp(subsettargeted.class.files,"Targeted")
  
  return(list(p1,p2,p3,p4))
}
