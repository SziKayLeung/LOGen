#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: characterise extent of NMD from SQANTI classification files, and associated with IR
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## associate_IR_NMD
## plot_IR_NMD_isoform_expression
## 

## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))


## ------------------- associate_IR_NMD

# Aim: plot a venn diagram of the number of genes with isoforms associated with IR and coding, NMD, and IR and NMD
# Function: subset_IR_NMD_knowngenes from diff_isoform_type_by_features.R
# Input:
  # class.files = df: classification file generated from SQANTI (further subsetted)
# Output:
  # venn diagram plot

associate_IR_NMD <- function(class.files){
  
  # subset class.files by IR and NMD categories
  cate <- subset_IR_NMD_knowngenes(class.files)
  
  # plot venn diagram
  myCol <- brewer.pal(3, "Set2")
  p <- venn.diagram(
    x = list(unique(cate$IR_coding$associated_gene), unique(cate$NMD$associated_gene), unique(cate$IR_NMD$associated_gene)),
    category.names = c("IR", "NMD", "IR-NMD"),
    filename = NULL,output=TRUE, print.mode = "raw",fill = myCol, cex = 1.5,fontface = "bold",fontfamily= "CM Roman",
    cat.cex = 1.5,cat.fontfamily = "CM Roman")
  
  return(p)
}


## ------------------- plot_IR_NMD_isoform_expression

# Aim: plot a box-plot of the expression of isoforms associated with IR and NMD
# Function: subset_IR_NMD_knowngenes from diff_isoform_type_by_features.R
# Input:
  # class.files = df: classification file generated from SQANTI (further subsetted)
# Output:
  # box-plot of the expression (Log10TPM) across isoform categories

plot_IR_NMD_isoform_expression <- function(class.files){
  
  # subset class.files by IR and NMD categories
  cate <- subset_IR_NMD_knowngenes(class.files)
  
  # box-plot
  p <- rbind(data.frame("Num" = cate$IR_NMD$Log_ISOSEQ_TPM, "Type" = "IR-NMD"),
              data.frame("Num" = cate$IR_notNMD$Log_ISOSEQ_TPM, "Type" = "IR"),
              data.frame("Num" = cate$notIR_NMD$Log_ISOSEQ_TPM,"Type" = "NMD"),
              data.frame("Num" = cate$notIR_notNMD$Log_ISOSEQ_TPM, "Type" = "Non IR-NMD")) %>%
    ggplot(., aes(y = Num, x = Type)) + geom_boxplot() +
    theme_bw() +
    labs(x = "", y = "Expression (Log10 TPM)") + mytheme +
    theme(legend.position = c(.90, 0.90), legend.title = element_blank()) 
  
  return(p)
}