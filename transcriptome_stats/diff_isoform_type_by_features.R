#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Observe if there is a difference in novel and known isoforms associated to known genes by expression, length, number of exons
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## subset_known_genes
## subset_isoforms_by_cate
## extract_feature
## subset_isoform_feature
## diff_isoform_type_by_features
##   
## ---------- Notes -----------------
## 
## All plots generated here are only of multi-exonic isoforms (mono-exonic isoforms are removed)
## Pre-requisite: 
## 1. SQANTI classification file has "Log_ISOSEQ_TPM" column


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- subset_known_genes

# Aim: subset SQANTI classification files by only known genes (i.e. remove novel genes)
# Input:
  # class.files = df: classification file generated from SQANTI
# Output:
  # list: 1. all isoforms associated with known genes
  #       2. all novel multi-exonic isoforms associated with known genes
  #       3. all known multi-exonic isoforms associated with known genes

subset_known_genes <- function(class.files){
  
  # all isoforms of known genes
  # capitalise associated_gene column to ensure not capturing any novel genes 
  cat("Subsetting file of isoforms associated with known genes:\nall\n")
  all = class.files[!grepl("NOVEL",toupper(class.files$associated_gene)),]
  
  # all coding isoforms of known genes
  cat("only coding\n")
  coding <- all %>% filter(coding == "coding")
  
  # novel multi-exonic isoforms of known genes
  cat("only novel\n")
  novel_transcripts = all[grepl("novel",all$associated_transcript),] %>% filter(subcategory != "mono-exon")
  
  # known multi-exonic isoforms of known genes
  cat("only known\n")
  known_transcripts = all[!grepl("novel",all$associated_transcript),] %>% filter(subcategory != "mono-exon") 
  
  # list output
  output <- list(all, novel_transcripts, known_transcripts, coding)
  names(output) = c("all","novel_multi_isoforms","known_multi_isoforms","coding")
  
  return(output)
}


## ------------------- subset_isoforms_bycate

# Aim: subset SQANTI classification files by isoform structural category
# Input:
  # class.files = df: classification file generated from SQANTI
  # type = str: structural_category in SQANTI classification file <FSM> <NIC> <ISM> <NNC> ...
# Output:
  # classification file subsetted by structural_category

subset_isoforms_bycate <- function(class.files, type){
  output <- class.files %>% filter(structural_category == type)
  return(output)
}


## ------------------- subset_IR_NMD_knowngenes

# Aim: subset SQANTI classification files of known genes (all known and novel isoforms) by intron retention and NMD
# Input:
# class.files = df: classification file generated from SQANTI
# Output:
# list: 1. all isoforms of known genes with intron retention and predicted for NMD
#       2. all isoforms of known genes with intron retention, but not predicted for NMD
#       3. all isoforms of known genes predicted for NMD
#       4. all isoforms of known genes predicted for NMD but not intron retention
#       5. all isoforms of known genes charaterised with intron retention and coding potential 
#       6. all isoforms of known genes not characterised with intron retention or NMD 

subset_IR_NMD_knowngenes <- function(class.files){
  
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  known.class.files <- subset_known_genes(class.files)
  
  # all isoforms of known genes with intron retention and predicted for NMD
  cat("subset isoforms with intron retention and predicted for NMD\n")
  known.IR.nmd <- known.class.files$all %>% filter(subcategory == "intron_retention" & predicted_NMD == "TRUE")
  
  # all isoforms of known genes with intron retention, but not predicted for NMD
  cat("subset isoforms with intron retention but not predicted for NMD\n")
  known.IR.notnmd <- known.class.files$all %>% filter(subcategory == "intron_retention" & predicted_NMD == "FALSE") 
  
  # all isoforms of known genes predicted for NMD
  cat("subset isoforms predicted for NMD\n")
  known.nmd <- known.class.files$all %>% filter(predicted_NMD == "TRUE")
  
  # all isoforms of known genes predicted for NMD but not intron retention
  cat("subset isoforms predicted for NMD, but not characterised with intron retention\n")
  known.notIR.nmd <- known.class.files$all %>% filter(predicted_NMD == "TRUE" & subcategory != "intron_retention")
  
  # all isoforms of known genes characterised with intron retention and coding potential 
  cat("subset isoforms characterised with intron retention and coding potential\n")
  known.IR.coding <- known.class.files$all %>% filter(subcategory == "intron_retention" & coding == "coding")
  
  # all isoforms of known genes not characterised with intron retention or NMD 
  known.notIR.notNMD <- known.class.files$all %>% filter(predicted_NMD == "FALSE" & subcategory != "intron_retention")
  
  # list output
  output <- list(known.IR.nmd, known.IR.notnmd, known.nmd, known.notIR.nmd, known.IR.coding, known.notIR.notNMD)
  names(output) = c("IR_NMD","IR_notNMD","NMD","notIR_NMD","IR_coding","notIR_notNMD")
  
  return(output)
}


## ------------------- extract_feature

# Aim: extract specific column (i.e. length, expression) and the transcripts based on the classification file inputted
# Input:
  # class.files = df: classification file generated from SQANTI (further subsetted)
  # feature = str: column name in classification file of interest (i.e lengths)
  # transcript_type = str: type of transcripts based on input of class.files (i.e. all "known" transcripts)
# Output:
  # dataframe with 2 columns: <feature> <Transcripts>

extract_feature <- function(class.files, feature, transcript_type){
  dat <- class.files %>% mutate(Transcripts = transcript_type) %>% .[,c(feature, "Transcripts")]
  return(dat)
}


## ------------------- subset_isoform_feature

# Aim: generate box-plots of the distribution of <features> (i.e. lengths) between novel and known isoforms
# Note: calls on functions noted above
# Input:
  # class.files = df: classification file generated from SQANTI (further subsetted)
  # col_name_feature = str: column name in SQANTI of interest (i.e. lengths)
  # y_label = str: name for y-axis
  # category = str: <all_novel> <all_split> 
      # all_novel = classify isoforms as either known or novel
      # all_split = classify isoforms by structural category (FSM, ISM..)
# Output: 1 box-plot 

subset_isoform_feature <- function(class.files, col_name_feature, ylabel, category){
  
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  known.class.files <- subset_known_genes(class.files)
  
  # y variable for plot
  y.var <- rlang::sym(quo_name(enquo(col_name_feature)))
  
  if(category == "all_novel") { 
    
    # extract all the known and novel isoforms with the specific feature
    set1 <- extract_feature(known.class.files$known_multi_isoforms, col_name_feature, "Known")
    set2 <- extract_feature(known.class.files$novel_multi_isoforms, col_name_feature, "Novel")
    
    # bind rows and plot
    merge <- bind_rows(list(set1, set2))
    
    p <- merge %>%  
      ggplot(., aes(x = Transcripts, y = !! y.var, fill = Transcripts)) + 
      geom_boxplot() + theme_bw() + 
      mytheme + labs(x = "", y = ylabel) + 
      scale_fill_manual(values = c(label_colour("known"),label_colour("novel"))) +
      theme(legend.position = "none") + 
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
    
  } else { 
    
    set1 <- extract_feature(subset_isoforms_bycate(known.class.files$known_multi_isoforms, "FSM"), col_name_feature, "FSM")
    set2 <- extract_feature(subset_isoforms_bycate(known.class.files$known_multi_isoforms, "ISM"), col_name_feature, "ISM")
    set3 <- extract_feature(subset_isoforms_bycate(known.class.files$novel_multi_isoforms, "NIC"), col_name_feature, "NIC")
    set4 <- extract_feature(subset_isoforms_bycate(known.class.files$novel_multi_isoforms, "NNC"), col_name_feature, "NNC")
    set5 <- extract_feature(subset_isoforms_bycate(known.class.files$novel_multi_isoforms, "Fusion"), col_name_feature, "Fusion")
    
    merge <- bind_rows(list(set1, set2, set3, set4, set5)) 
    merge$Transcripts <- factor(merge$Transcripts, levels = c("FSM", "ISM", "NIC", "NNC","Fusion"))
    
    p <- merge %>%  
      ggplot(., aes(x = Transcripts, y = !! y.var, fill = Transcripts)) + 
      geom_boxplot() + theme_bw() + 
      scale_fill_manual(values = c(alpha(label_colour("known"),0.8),alpha(label_colour("known"),0.3),
                                   alpha(label_colour("novel"),0.8),alpha(label_colour("novel"),0.5),alpha(label_colour("novel"),0.3))) +
      mytheme + labs(x = "", y= ylabel) + 
      theme(legend.position = "none") + 
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
  }
  return(p)
}


## ------------------- diff_isoform_type_by_features

# Aim: wrappyer to generate box-plots of the distribution of 3 features (expression, isoform length, number of exons) between novel and known isoforms
# Note: calls on functions noted above; 
# Prequisite: class.files has "Log_ISOSEQ_TPM", "length", "exons" columns 
# Input:
  # class.files = df: classification file generated from SQANTI 
# Output: 1 box-plot 

diff_isoform_type_by_features <- function(class.files){
  
  # expression
  expression <- subset_isoform_feature(class.files, "Log_ISOSEQ_TPM", "Expression (Log10 TPM)", "all_novel")
  split_expression <- subset_isoform_feature(class.files, "Log_ISOSEQ_TPM", "Expression (Log10 TPM)", "all_split")
  
  # isoform  length
  length <- subset_isoform_feature(class.files, "length", "Transcript Length (kb)", "all_novel") + scale_y_continuous(labels = ks)
  split_length <- subset_isoform_feature(class.files, "length", "Transcript Length (kb)", "all_split") + scale_y_continuous(labels = ks)
  
  # number of exons
  exon <- subset_isoform_feature(class.files, "exons", "Number of Exons", "all_novel")
  split_exon <- subset_isoform_feature(class.files, "exons", "Number of Exons", "all_split")
  
  output <- list(expression,split_expression,length,split_length,exon,split_exon)
  
  return(output)
}