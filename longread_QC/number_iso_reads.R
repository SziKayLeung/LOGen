#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Basic QC of Iso-Seq pipeline 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## number_iso_reads
## isoseq_QC_yield
##   
## ---------- Notes -----------------
##
## Pre-requisite: 
## 1. df of [sample, run_id] 
## 2. ccs and lima ouputs of reads generated from ccs.py and lima.py 
## 3. df of [sample, Total.Bases..GB., Genotype,RIN]

## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))


## ------------------- number_iso_reads

# Aim: huge wrapper function to read in multiple files from CCS, LIMA and REFINE directory, and output 5 plots
# Input:
  # sample =  df of 2 columns - sample & run_id
  # ccs_output = df: output of the number of ccs reads processed
  # lima_output = df: output of the number of lima reads processed
  # refine_dir; cluster_dir = path:  of refine and cluster directory containing ".json" and "cluster_report.csv" files

number_iso_reads <- function(ccs_output, lima_output, refine_dir, cluster_dir){
  
  ## Input CCS, LIMA, REFINE files
  # REFINE summary input
  REFINE_json_list <- list.files(paste0(refine_dir), pattern = "flnc.filter_summary.json", full.names = T)
  REFINE_list <- lapply(REFINE_json_list , function(x) as.data.frame(fromJSON(file = x)))
  names(REFINE_list) <- list.files(paste0(refine_dir), pattern = "flnc.filter_summary.json")
  for(i in 1:length(names(REFINE_list))){
    REFINE_list[[i]]$Sample <- names(REFINE_list)[[i]]
  }
  
  # CLUSTER summary input
  CLUSTER_list_names <- list.files(paste0(cluster_dir), pattern = ".cluster_report.csv", full.names = T)
  CLUSTER_list <- lapply(CLUSTER_list_names, function(x) read.csv(x))
  names(CLUSTER_list) <- list.files(paste0(cluster_dir), pattern = ".cluster_report.csv")
  for(i in 1:length(names(CLUSTER_list))){
    CLUSTER_list[[i]]$Sample <- names(CLUSTER_list)[[i]]
  }
  

  ####################### CCS
  # Extract only values from mix of values and percentage
  CCS_values <- cbind(as.character(ccs_output[,1]), apply(ccs_output[,-1], 2, function(x) word(x, c(1), sep = fixed("("))))
  colnames(CCS_values)[1] <- "Description"
  CCS_values_mod <- as.data.frame(CCS_values) %>% reshape2::melt(., id = "Description") %>%
    mutate(sample = word(variable, c(1), sep = fixed("_")))
  
  
  ####################### LIMA
  # Extract only values from the mix of values and percentage from LIMA output
  LIMA_values <- data.frame(lima_output[,1],lapply(lima_output[,2:ncol(lima_output)], function(x) as.numeric(word(x, c(1),  sep = fixed ('(')) )))
  colnames(LIMA_values)[1] <- "Description"
  LIMA_values_mod <- as.data.frame(LIMA_values) %>% reshape2::melt(., id = "Description") %>%
    mutate(sample = word(variable, c(1), sep = fixed("_")))
  
  
  # Refine into one dataframe and data-wrangled for easy merging
  REFINE <- do.call("rbind", REFINE_list) %>%
    rownames_to_column(., var = "variable") %>%
    mutate(sample = word(variable, c(1), sep = fixed("."))) %>%
    as.data.frame() %>%
    gather(., Description, value, 2:4, factor_key=TRUE) %>%
    .[,c("Description","variable","value","sample")] %>%
    mutate(value = as.numeric(value))
  
  
  CLUSTER <- do.call("rbind",CLUSTER_list) %>% 
    rownames_to_column(., var = "variable") %>%
    mutate(sample = word(variable, c(1), sep = fixed("."))) %>%
    as.data.frame() %>% 
    group_by(cluster_id,sample) %>%
    tally() %>% 
    group_by(sample) %>% tally() %>% mutate(Description = "Clustered_transcripts") %>% mutate(variable = "clustered.csv") %>% 
    .[,c(3,4,2,1)] 
  colnames(CLUSTER)[3] <- "value"
  
  Reads <-   
    CCS_values_mod[CCS_values_mod$Description == "ZMWs input               ",] %>%
    bind_rows(CCS_values_mod[CCS_values_mod$Description == "ZMWs pass filters        ",],)  %>%
    mutate(value = as.numeric(value)) %>%
    bind_rows(REFINE) %>%
    bind_rows(CLUSTER) %>%
    mutate(Description = as.character(Description))
  
  Reads$Description <- revalue(Reads$Description, c("ZMWs input               "="Polymerase Reads", "ZMWs pass filters        "="CCS Reads",
                                                    "num_reads_fl"="FL Reads", "num_reads_flnc"="FLNC reads",
                                                    "num_reads_flnc_polya" = "Poly-A FLNC reads",
                                                    "Clustered_transcripts" = "Transcripts"))
  levels(Reads$Description) <- c("Polymerase Reads","CCS Reads","FL Reads","FLNC reads","Poly-A FLNC reads","Transcripts")
  
  ### To calculate proportions to generate plots
  # total failed ccs reads
  failed_CCS_reads <- as.data.frame(CCS_values) %>% reshape2::melt(., id = "Description") %>% filter(Description %in% c("ZMWs filtered       (C)  "))
  failed_LIMA_reads <- as.data.frame(LIMA_values) %>% reshape2::melt(., id = "Description") %>% filter(Description %in% c("ZMWs below any threshold  (C) "))
  
  
  output <- list(Reads,CCS_values_mod)
  names(output) <- c("Reads","CCS_values")
  return(output)
}


## ------------------- plot_transcripts_postisoseq

# Aim: plot the nummber of full-length transcripts by genotype
# Input: reads_df from wrangled output of number_iso_reads
  # see examples in isoseq_QC_yield, isoseq_QC_yield_batch
# Output: box-plot
plot_transcripts_postisoseq <- function(reads_df){
  
  p <- reads_df %>%
    filter(Description == "Transcripts", !is.na(Genotype)) %>% 
    ggplot(., aes(x = Genotype, y = value, colour = Genotype)) +
    geom_boxplot() + geom_point() +  mytheme + 
    labs(x = "", y = "Number of FL Transcripts (K)") +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) + theme(legend.position = "none") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3,accuracy = 1))
  
  return(p)
}


## ------------------- isoseq_QC_yield

# Aim: read in sequencing data (user-generated) and output from number_of_reads to generate QC related plots
# Input:
  # reads_df = df: number_of_reads function output 1
  # sequenced = df: containing [sample,Total.Bases..GB.,Genotype,RIN]
# Output
  # p1: box-plot - Number of bases by genotype
  # p2: scatter-plot - Number of bases vs RIN
  # p3: line-plot - Number of bases processed through the Iso-Seq pipeline
  # p4: box-plot - Number of transcripts by genotype

iso_QC_yield <- function(reads_df, sequenced){
  
  p1 <- sequenced %>% ggplot(., aes(x = Genotype, y = Total.Bases..GB., fill = Genotype)) + geom_boxplot() +
    geom_jitter(shape=17, position=position_jitter(0)) +
    theme_bw() +
    labs(x = "Genotype", y = "Total Bases (Gb)") +
    mytheme + theme(legend.position = "none") +
    scale_fill_manual(values = c(label_colour("TG"),label_colour("WT")))
  
  p2 <- sequenced %>% ggplot(., aes(x = RIN, y = Total.Bases..GB., color = Genotype)) + geom_point(size = 2, shape = 17) +
    scale_color_manual(values = c(label_colour("TG"),label_colour("WT"))) +
    mytheme + labs(x = "RIN", y = "Total Bases (Gb)") 
  
  p3 <-
    reads_df %>%
    filter(Description != "Transcripts") %>% 
    ggplot(., aes(x = reorder(Description, -value), y = value, group = sample, colour = Genotype)) +
    geom_line() + geom_point() +  mytheme + theme(legend.position = "top") + labs(x = "Reads", y = "Number of Reads (K)") +
    scale_color_manual(values = c(label_colour("WT"),label_colour("TG"))) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    scale_x_discrete(labels=c("Polymerase Reads" = "Polymerase", "CCS Reads" = "CCS","FL Reads" = "FL","FLNC reads" = "FLNC",
                              "Poly-A FLNC reads" = "Poly-A FLNC"))
  
  p4 <- plot_transcripts_postisoseq(reads_df)
  
  ### Correlations 
  # Difference between WT and TG yield
  #with(sequenced, shapiro.test(Total.Bases..GB.[Genotype == "WT"])) # p = 0.4
  #with(sequenced, shapiro.test(Total.Bases..GB.[Genotype == "TG"])) # p = 0.7
  #res.ftest <- var.test(Total.Bases..GB.~ Genotype , data = sequenced)
  #res.ftest # p = 0.29
  res <- t.test(RIN ~ Genotype , data = sequenced, var.equal = TRUE)
  res
  
  # correlation of RIN and total yield
  cor.test(sequenced$RIN,sequenced$Total.Bases..GB.)
  
  res <- t.test(Total.Bases..GB.~ Genotype , data = sequenced, var.equal = TRUE)
  res
  
  return(list(p1,p2,p3,p4))
}


## ------------------- isoseq_QC_yield_batch

# Aim: generate similar plots to isoseq_QC_yield, but for targeted datasets (batched samples)
# Input:
  # reads_df = df: number_of_reads function output 1
  # targetedpheno = df: phenotype file of barcoded samples and batches <Sample> <Phenotype> <Batch> <Barcode>
  # samples = df: list of samples
# Output
  # p1: line-plot - Number of bases processed through the Iso-Seq pipeline
  # p2: box-plot - Number of transcripts by genotype
  # p3: box-plot - Number of PolyA Flnc reads by genotype across batch

iso_QC_yield_batch <- function(Reads, targetedpheno, samples){
  
  Reads_plot <- Reads %>% mutate(Batch = word(Reads$variable, c(3), sep = fixed("_"))) %>% 
    full_join(., targetedpheno, by = c("sample" = "Sample")) %>%
    unite("Batch", Batch.x,Batch.y, na.rm = TRUE) %>%
    mutate(Batch = recode(Batch, "3b" = "3", "3a" = "3 (partial run)")) %>%
    mutate(Genotype = Phenotype)
  
  p1 <- Reads_plot %>% filter(Description != "Transcripts") %>% filter(Batch != "3 (partial run)") %>%
    ggplot(., aes(x = Description, y = value, colour = Batch, group = Batch)) +
    geom_line() + geom_point(size = 3) +  mytheme + theme(legend.position = c(0.8,0.8)) + labs(x = "", y = "Number of Reads (Thousands)") +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_colour_discrete(name = "",labels = c("Batch 1 (n = 6)","Batch 2 (n = 9) ","Batch 3 (n = 9)"))
  
  p2 <-  plot_transcripts_postisoseq(Reads_plot)
  
  p3 <- Reads_plot[Reads_plot$Description == "Poly-A FLNC reads" ,] %>% 
    filter(!is.na(Genotype)) %>%
    ggplot(., aes(x = Batch, y = value)) + geom_boxplot() + geom_point(aes(colour = Genotype),size = 3) + mytheme +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    labs(x = "Batch", y = "Number of Poly-A FLNC Reads (Thousands)") +
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) 
  
  cat("Number of CCS Reads in Batch:\n")
  Reads[Reads$Description == "CCS Reads",]
  
  cat("Poly-A FLNC Reads across all Batches:\n")
  polyA_flnc_reads <- Reads_plot[Reads_plot$Description == "Poly-A FLNC reads","value"]/1000
  cat("Sum number (Thousand):", sum(polyA_flnc_reads),"\n")
  cat("Mean number (Thousand):", mean(polyA_flnc_reads),"\n")
  cat("Min number (Thousand):", min(polyA_flnc_reads),"\n")
  cat("Max number (Thousand):", max(polyA_flnc_reads),"\n")
  for(i in 1:3){cat("Sum number (Thousand) of Poly-A FLNC Reads in Batch",i,":", 
                    sum(Reads_plot[Reads_plot$Description == "Poly-A FLNC reads" & Reads_plot$Batch == i,"value"])/1000,"\n")}
  
  # not normally distributed therefore wilcoxon rank sum test
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "WT"]))
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "TG"]))
  #var.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) #cannot assume variance
  wilcox.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) 
  
  # correlation of FL transcripts and RIN
  transcript_RIN <- merge(Reads_plot %>% filter(Description == "Transcripts"), samples, by.x = "sample", by.y = "Sample.ID", all.x = T)
  #shapiro.test(transcript_RIN$value) # spearman's rank
  #shapiro.test(transcript_RIN$RIN) 
  cor.test(transcript_RIN$value,transcript_RIN$RIN, method = "spearman", exact = FALSE)
  
  return(list(p1,p2,p3))
}
