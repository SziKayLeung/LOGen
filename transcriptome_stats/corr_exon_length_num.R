#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: correlate length and number of exons with number of isoforms 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## gene_summary_info
## corr_exon_length_num
##   
## ---------- Notes -----------------
## 
## Correlation reported is only of multi-exonic isoforms (mono-exonic isoforms are removed)
## Pre-requisite: 
## 1. SQANTI classification file has a tabulated "FL" column for each isoform with sum reads across each sample
## 2. function: draw_density.R from aesthetics_basics_plots.R


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- gene_summary_info

# Aim: summarise the number of isoforms, minimum number of exons and max etc per gene 
# Input: 
  # class.files = df: classification file generated from SQANTI
# Output: 
  # df of basic information for each gene 

gene_summary_info <- function(class.files){
  
  total_fl <- sum(class.files$FL, na.rm=T)
  
  info <- list(
    # Number of isoforms 
    class.files %>% dplyr::count(associated_gene, name = "Num_of_Isoforms"),
    # Min Exons 
    class.files %>% group_by(associated_gene) %>% dplyr::summarise(Min_exon = min(exons)),
    # Max Exons 
    class.files %>% group_by(associated_gene) %>% dplyr::summarise(Max_exon = max(exons)),
    # Min Length 
    class.files %>% group_by(associated_gene) %>% dplyr::summarise(Min_length = min(length)),
    # Max Length
    class.files %>% group_by(associated_gene) %>% dplyr::summarise(Max_length = max(length)), 
    # Num of FL reads 
    class.files %>% group_by(associated_gene) %>% dplyr::summarise(Total_Reads = sum(FL))
  )
  
  final <- Reduce(function(...) merge(..., by='associated_gene', all.x=TRUE), info) %>%
    filter() 
  
  # TPM by gene (sum of FL reads)
  final <- final %>% mutate("FL_TPM" = round(Total_Reads*(10**6)/total_fl)) %>% mutate("Log10_FL_TPM" = log10(FL_TPM))
  
  return(final)
}


## ------------------- corr_exon_length_num

# Aim: Correlate the number of exons/length with number of isoforms (only multiexonic only)
# Input:
  # class.files = df: classification file generated from SQANTI
# Output:
  # Mulitple plots and .txt file for the output of correlation of stats
  # P1: transcript length vs number of exons for all transcripts 
  # P2: transcrpt length vs number of exons (using only representative transcript per gene)
  # P3: same plot as P2 but filtered at expression level 
  # P4: transcript length vs number of isoforms (using only representative length) 
  # P5: same plot as P4 but filtered at expression level 
  # P6: exon vs number of isoforms (exon defined by representative transcript) 
  # P7: same plot as P7 but filtered at expression level 

corr_exon_length_num <- function(class.files){
  
  # remove mono-exonic transcripts 
  dat1 <- class.files %>% filter(subcategory != "mono-exon")
  
  # gene_summary_info = function to summarise the number of isoforms, etc
  # note TPM calculated with the removal of monoexonic transcripts 
  dat2 <- data.frame(gene_summary_info(dat1))
  
  # Gene expression cut off threshold at Log10_FL_TPM at >2.5 
  dat3 <- dat2 %>% filter(Log10_FL_TPM > 2.5)
  
  
  # Plots 
  p1 <- density_plot(dat1, "length","exons","Transcript length (kb)","Number of exons", 
                     "Transcript length\nvs Number of Exons (all transcripts)")  + labs(title="\n\n\n")
  
  p2 <- density_plot(dat2, "Max_length","Max_exon", "Transcript length (kb)", "Number of exons",
                     "Transcript length \nvs Number of Exons (representative transcript per gene)") + labs(title="") +
    scale_x_continuous(labels = ks)
  
  p3 <- density_plot(dat3, "Max_length","Max_exon", "Transcript length (bases)", "Number of exons",
                     "Transcript Length\nvs Number of Exons (representative transcript per gene):\nFiltered at 25TPM")  + labs(title="\n\n\n")
  
  # Transcript length (max representative)\nvs Number of Isoforms
  p4 <- density_plot(dat2, "Max_length","Num_of_Isoforms", 
                     "Transcript Length (kb)", "Number of Isoforms", "Transcript length (max representative) vs Number of Isoforms") + 
    scale_x_continuous(labels = ks)  + labs(title="")
  
  # Transcript length (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p5 <- density_plot(dat3, "Max_length","Num_of_Isoforms", 
                     "Gene Length (kb)", "Number of Isoforms",
                     " Transcript length (max representative) vs Number of Isoforms:Filtered for high gene expression") + 
    scale_x_continuous(labels = ks) + labs(title="\n\n\n")
  
  # Number of Exons (max representative)\nvs Number of Isoforms
  p6 <- density_plot(dat2, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms") + labs(title="")
  
  # Number of Exons (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p7 <- density_plot(dat3, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms:Filtered for high gene expression") + labs(title="\n\n\n")
  
  
  # linear regression 
  # summary(lm(Num_of_Isoforms ~ Max_length + Max_exon, dat2))
  return(list(p1,p2,p3,p4,p5,p6,p7))
}