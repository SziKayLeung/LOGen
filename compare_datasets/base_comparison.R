#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: basic functions for dataset comparisons 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## compare_isoform
## corr_gene_exp_toknown
## plot_iso_length_mdatasets 
## rnaseq_support_mdatasets
## common_vs_unique_iso
##


compare_isoform <- function(ref_gencode, class_dat, dataset){
  df = merge(ref_gencode, class_dat %>% group_by(associated_gene) %>% tally(), by = "associated_gene") %>% 
    mutate(Number_of_Gencode_Isoforms = log10(Number_of_Gencode_Isoforms), n = log10(n), MeanRNASeqCounts = log10(MeanRNASeqCounts))
  
  p1 <- density_plot(df, "Number_of_Gencode_Isoforms","n", "Number of Gencode Isoforms", paste("Number of",dataset,"Isoforms"), "")
  p2 <- density_plot(df, "MaxGeneLength","n", "Gene Length", paste("Number of",dataset,"Isoforms"), "")
  p3 <- density_plot(df, "MaxTransLength","n", "Transcript Length", paste("Number of",dataset,"Isoforms"), "")
  p4 <- density_plot(df, "MeanRNASeqCounts","n", "Gene Expression", paste("Number of",dataset,"Isoforms"), "")
  p5 <- density_plot(df, "Maxexons","n", "Number of Exons", paste("Number of",dataset,"Isoforms"), "")
  
  count = 1
  corr_output = data.frame()
  for(x in c("Number_of_Gencode_Isoforms","MaxGeneLength","MaxTransLength","MeanRNASeqCounts","Maxexons")){
    corr_output[count,1] = round(cor(df[[x]],df$n, use = "pairwise.complete.obs"),2)
    count = count + 1
  }
  colnames(corr_output)[1] = dataset
  
  output = list(p1,p2,p3,p4,p5, corr_output)
  names(output) = c("Isodiff","Genelendiff","Translendiff","Geneexpdiff","Exondiff","Corr")
  
  return(output)
}


## ------------------- corr_gene_exp_toknown

# Aim: correlate the gene expression of dataset to the known expression (reference/RNA-Seq)
# Note:   
  # Gene expression is determined from sum of "FL" across associated isoforms for each gene
# Input:
  # ref_gencode = df: contains the MeanRNASeqCounts column of associated_genes of interest
  # class.files = df: SQANTI classification file
  # dataset = str: axis label
# Output:
  # p1: density plot and correlation of the gene expression 

corr_gene_exp_toknown <- function(ref_gencode, class.files, dataset){
  df_gene_exp = class.files %>% group_by(associated_gene) %>% tally(FL)
  
  merged_plot = merge(ref_gencode,df_gene_exp, by = "associated_gene") %>% 
    mutate(logFL = log2(n), MeanRNASeqCounts = log2(MeanRNASeqCounts))
  
  p <- density_plot(merged_plot,"logFL","MeanRNASeqCounts", paste(dataset, "FL Gene Expression (Log2)"),"RNA-Seq Gene Expression (Log2)","")
  return(p)
}


## ------------------- plot_iso_length_mdatasets 

# Aim: plot the distribution of y variable (i.e. transcript length, number of exons) in dataset
# Input:
  # class.files = df: SQANTI classification file 
  # y.var = str: y-variable of interest - ensure variable exists as column name in dataset
  # y.name = str: name of y-variable 
# Output: 
  # p: box-plot of the y-variable distribution

plot_iso_length_mdatasets <- function(class.files, y.var, y.name){
  
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  p <- ggplot(class.files, aes(x = Dataset, y = !! y.var, fill = Dataset)) + geom_boxplot() + 
    mytheme + labs(x = " ", y = y.name) +
    theme(legend.position = "None")
  
  return(p)
}


## ------------------- rnaseq_support_mdatasets

# Aim: plot RNA-Seq support distribution across multiple datasets
# Function:
  # plot_stats_feature_2datasets() from identify_within_cage_SS_peaks.R
# Note:
  # Function relies on "RNASeq_supported" column in classification files (<Yes> <No> based on min_coverage)
  # "RNASeq_supported" column generated from annotate_class_binary()
  # Input should be from a merged dataset using gffcompare for comparison, hence no redundant isoforms
  # Input needs a "Dataset" column for identifier; Both/<dataset1>/<dataset2>
# Input:
  # merged.class.files = df: merged classification file from two datasets
# Output:
  # p1: bar-plot of the number of isoforms with RNA-Seq support in both or unique to dataset
  # p2: bar-plot of the percentage of isoforms with RNA-Seq support in both or unique to dataset
  # p3: bar-plot/stats of number of isoforms unique to dataset

rnaseq_support_mdatasets <- function(merged.class.files){
  
  # tally number and proportions of RNA-Seq support across datasets 
  prop <- merge(merged.class.files  %>% group_by(Dataset,RNASeq_supported) %>% tally(), 
                merged.class.files %>% group_by(Dataset) %>% tally(), by = "Dataset") %>% 
    mutate(perc = n.x/n.y * 100) %>% 
    filter(!is.na(RNASeq_supported))
  
  # remove "NA" RNASeq support which refers to monoexons
  p1 <- ggplot(merged.class.files, aes(x = Dataset, fill = RNASeq_supported)) + 
    geom_bar(stat = "count") +
    mytheme + labs(x = " ", y = "Number of Isoforms") + theme(legend.position = "top") +     
    scale_fill_manual(values = c(label_colour("No"),label_colour("Yes")), name = "RNA-Seq")
  
  p2 <- ggplot(prop, aes(x = Dataset, y = perc, fill = RNASeq_supported)) + 
    geom_bar(stat = "identity") + 
    mytheme + labs(x = " ", y = "Isoforms (%)") + theme(legend.position = "top") +     
    scale_fill_manual(values = c(label_colour("No"),label_colour("Yes")), name = "RNA-Seq")

  p3 <- plot_stats_feature_2datasets(merged.class.files[merged.class.files$Dataset != "Both",],"RNASeq_supported","RNASeq")

  return(list(p1,p2,p3))
}


## ------------------- common_vs_unique_iso

# Aim: plot the expression of commonly and unique isoforms 
# that are supported/not supported by RNA-Seq or within/not within cage-peak 
# Input:
  # merged.class.files = df: classification file of merged dataset using gffcompare for comparison, hence no redundant isoforms
  # merged.class.files.exp = df: df with expression of the isoforms, "FL", and "LogFL" 
  # feature = str: <RNASeq_supported> <within_cage_peak>
  # x_lab = str: label for x-axis
  # dataset = str: dataset name of interest for plotting the FL expression
# Note: 
  # requires <dataset>LogFL, <dataset>_isoform columns in merged.class.files.exp
# Output:
  # box-plot of the isoform expression (Log10FL) of common and unique isoforms that are either supported by RNA-Seq/within CAGE peak

common_vs_unique_iso <- function(merged.class.files, merged.class.files.exp, feature, x_lab, dataset){
  
  # merge class.files and clas.files.exp for the FL abundance and features column
  df = merge(merged.class.files, merged.class.files.exp, by = c("isoform" = paste0(dataset,"_isoform")))
  
  # remove NA
  if(feature == "RNASeq_supported"){df = df %>% filter(!is.na(RNASeq_supported))}
  else if(feature == "within_cage_peak"){df = df %>% filter(!is.na(within_cage_peak))}
  else{print("NULL")}
  
  # plot
  x_var <- rlang::sym(quo_name(enquo(feature)))
  y_var <- rlang::sym(quo_name(paste0(dataset,"_LogFL")))
  
  p <- ggplot(df, aes(x = !! x_var, y = !! y_var, fill = Dataset.x)) + geom_boxplot() +
    mytheme + labs(x = x_lab, y = paste(dataset, "Isoform Expression (Log10FL)")) + theme(legend.position = "top") +
    scale_fill_manual(values = c(label_colour("BothTech"),label_colour(dataset)), name = "Dataset")
  
  return(p)
}