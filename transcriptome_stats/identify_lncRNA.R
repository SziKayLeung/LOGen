#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Observe difference between long-non-coding-RNA (lncRNA) and non-long-non-coding-RNA (non lncRNA)
## by expression, length, number of exons ..
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## plot_lncRNA_vs_non
## identify_lncRNA
##   
## ---------- Notes -----------------
## 
## Pre-requisite: 
## 1. Ran alignment of dataset to long-non-coding reference genome, followed by SQANTI
## 2. function: diff_isoform_type_by_features.R from transcriptome_stats


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
source(paste0(LOGEN,"/transcriptome_stats/diff_isoform_type_by_features.R"))


## ------------------- plot_lncRNA_vs_non

# Aim: lncRNA vs non-lncRNA across different features (i.e. expression, exons)
# Input:
  # nonlnc_vec = vec: vector of non long-coding RNA feature (numeric) (i.e expression)
  # lnc_vec = vec: vector of long-coding RNA feature (numeric) (i.e. expression)
  # feature = <expression; exons; orf_length; isoform_num;length> 
# Output depending on the feature:
  # expression: box-plot of the expression between lncRNA and non-lncRNA isoforms
  # exons: box-plot of the number of exons between lncRNA and non-lncRNA isoforms
  # orf_length: density plot of open reading frame length between lncRNA and non-lncRNA isoforms
  # isoform_num: box plot of the number of isoforms between lncRNA and non-lncRNA isoforms
  # length: density plot of the isoform length between lncRNA and non-lncRNA isoforms
  
plot_lncRNA_vs_non <- function(nonlnc_vec, lnc_vec, feature){
  cat("Subsetting by:", feature, "\n")
  
  # combine two vectors into a dataframe, delimited by "Isoform"
  df <- rbind(data.frame("Num" = nonlnc_vec, "Isoform" = "Non-lncRNA"),
        data.frame("Num" = lnc_vec, "Isoform" = "lncRNA"))
  
  if(feature == "expression"){
   
     p <- ggplot(df , aes(y = Num, x = Isoform)) + geom_boxplot() +
      labs(x = "", y = "Expression (Log10 TPM)") +
      scale_fill_manual(values = wes_palette("Darjeeling2")) +
      theme(legend.position = c(.15, 0.85)) 
  
  }else if(feature == "exons"){
    
    p <-  ggplot(df, aes(Num, x = Isoform)) + geom_boxplot() +
      labs(y = "Number of Exons", x = "") +
      mytheme + theme(legend.position = c(.15, 0.85))
    
  }else if(feature == "orf_length"){
    
    p <- ggplot(df, aes(Num, fill = Isoform)) + geom_density(alpha = 0.2) +
      labs(x = "ORF Length (kb)", y = "Density") +
      theme(legend.position = c(.15, 0.85)) +
      scale_x_continuous(labels = ks)
  
  }else if(feature == "isoform_num"){
    
    p <- ggplot(df, aes(y = Num, x = Isoform)) + geom_boxplot() +
      labs(x = "", y = "Number of Isoforms (Log10)") +
      theme(legend.position = c(.25, 0.85)) 
  
  }else if(feature == "length"){
    
    p <- ggplot(df, aes(Num, fill = Isoform)) + geom_density(alpha = 0.2) +
      labs(x = "Transcript Length (kb)", y = "Density") +
      theme(legend.position = c(.45, 0.85)) +
      scale_x_continuous(labels = ks) + theme(legend.position="none") 
    
  }else{
    print("feature term is incorrect: must be <expression; exons; orf_length; isoform_num;length>")
  }
  
  p <- p + mytheme + theme(legend.direction="horizontal")
  return(p)
}


## ------------------- identify_lncRNA

# Aim: wrapper to plot lncRNA vs non lncRNA
# lncRNA genes are called from lnc.class.files, which are then used to subset the original classification file by gene name
# (i.e. able to identify lncRNA genes from alignment to lncRNA reference genome)
# Pre-requisite: annotate isoforms with SQANTI using long non-coding RNA reference genome rather than standard reference
# Input:
  # class.files = df: classification file generated from SQANTI after alignment to standard reference genome
  # lnc.class.files = df: classification file generated from SQANTI after alignment to non-coding reference genome
# Output:
  # p1: Transcripts Length (Annotated Genes only)
  # p2: Number of Exons (Annotated Genes only)
  # p3: Non-Coding Transcripts Length (Annotated Genes only)
  # p4: Coding Transcripts Length (Annotated Genes only)
  # p5: Transcript Expression (Annotated Genes only)
  # p6: ORF Length (Annotated Genes only)
  # p7: Number of Isoforms per Gene

differentiate_lncRNA <- function(class.files, lnc.class.files){
  
  # subset classification files (filtered) to known genes
  # no point working with novel genes
  known.class.files <- subset_known_genes(class.files)

  # vector of the known lncRNA genes from lnc.files
  lncrna_genes <- unique(lnc.class.files[!grepl("novel",lnc.class.files$associated_gene),"associated_gene"])
  
  # subset of lncRNA transcripts from original classification file using genes extracted from lnc.files 
  isoform <- list(
    nonlnc = known.class.files$all[!known.class.files$all$associated_gene %in% lncrna_genes,],
    lnc = known.class.files$all[known.class.files$all$associated_gene %in% lncrna_genes,]
  )
  
  # coding potential 
  coding <- list(
    nonlnc_no = isoform$nonlnc[isoform$nonlnc$coding == "non_coding",],
    nonlnc_yes = isoform$nonlnc[isoform$nonlnc$coding == "coding",],
    lnc_no = isoform$lnc[isoform$lnc$coding == "non_coding",],
    lnc_yes = isoform$lnc[isoform$lnc$coding == "coding",]
  )

  # number of isoforms (n) of lncRNA genes vs nonlncRNA genes for each gene
  isoform_diversity <- lapply(isoform, function(x) x %>% group_by(associated_gene) %>% tally())

  
  # Plots
  p1 <- plot_lncRNA_vs_non(isoform$nonlnc$length, isoform$lnc$length, "length") 
  
  p2 <- plot_lncRNA_vs_non(isoform$nonlnc$exons, isoform$lnc$exons, "exons") 
  
  p3 <- plot_lncRNA_vs_non(coding$nonlnc_no$length, coding$lnc_no$length, "length")
  
  p4 <- plot_lncRNA_vs_non(coding$nonlnc_yes$length, coding$lnc_yes$length, "length")
  
  p5 <- plot_lncRNA_vs_non(isoform$nonlnc$Log_ISOSEQ_TPM, isoform$lnc$Log_ISOSEQ_TPM, "expression") 
  
  p6 <- plot_lncRNA_vs_non(isoform$nonlnc$ORF_length, isoform$lnc$ORF_length, "orf_length") 
  
  p7 <-plot_lncRNA_vs_non(log10(isoform_diversity$nonlnc$n), log10(isoform_diversity$lnc$n),"isoform_num") 
  
  plots <- list(p1,p2,p3,p4,p5,p6,p7)
  return(plots)
}