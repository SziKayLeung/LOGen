#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Determine on-target rate of long-read targeted experiments 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## read_batch_readcounts
## summarise_ont_seq_reads
## batch_summarise_ont_seq_reads
## summarise_ont_demux_reads
## number_ont_reads 
##   
## ---------- Notes -----------------
##
## Pre-requisite: 
## 1. generated Read_Counts.txt files for each ONT processed stage (following bash command)
##


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(wesanderson))


## ------------------- read_batch_readcounts

# Aim: read in ReadCount.txt recording the number of ONT reads at each processed stage
# Pre-requisite:
  # Generate Batch<...>ReadCount.txt in correct folders (note ... should refer to batch number)
# Input:
  # raw_dir = str: path of raw directory containing Merged_ReadCount.txt of the number of basecalled reads  
  # root_dir = str: path of root directory containing 1_nanofilt and 2_demultiplex of number of filtered and demuxed reads
  # bnum1 = start of batch number (i.e. 2)
  # bnum2 = end of batch number (i.e. 3)
# Note: 
  # Function assumes a sequential processing of batch numbers with no skipping (i.e Batch 1, 2, 3...)
# Output: 
  # list of ReadCounts = level1: Basecalled, Filtered, Demultiplexed, level2: Batch ...

read_batch_readcounts <- function(raw_dir, root_dir, bnum1, bnum2){
  
  # Loop through each batch to record the number of reads
  ReadCounts <- list()
  for(i in bnum1:bnum2){
    print(i)
    ReadCounts$Basecalled[[i]] = read.table(paste0(raw_dir,"Batch",i,"_Merged_ReadCount.txt")) %>% mutate(Batch = i)
    ReadCounts$Filtered[[i]] = read.table(paste0(root_dir,"1_nanofilt/Batch",i,"_Filtered_ReadCount.txt")) %>% mutate(Batch = i)
    ReadCounts$Demultiplexed[[i]] = read.table(paste0(root_dir,"2_demultiplex/Batch",i,"_Demultiplex/Batch",i,"_Reads_Count.txt")) %>% mutate(Batch = i)
  }
  
  # remove empty list (generated if bnum1 != 1)
  for(i in 1:length(ReadCounts)){ReadCounts[[i]] <-  Filter(length, ReadCounts[[i]])}
  
  return(ReadCounts)
}


## ------------------- summarise_ont_seq_reads

# Aim: create a table of the number of readcounts per sequencing run 
# Input:
  # sequencing_data = list: output of prepare_summary_file function after parsing through the ONT summary .txt file
  # batchname = str: name of the batch for column
# Output: df

summarise_ont_seq_reads <- function(sequencing_data, batchname){
  dat <- data.frame(ReadCount = c(nrow(sequencing_data$sequencedata),
                                  nrow(sequencing_data$passedSeqs)),
                    Batch = c(str_remove(batchname,"B"),str_remove(batchname,"B")),
                    Type = c("Sequenced","Basecalled"))
  
  return(dat)
}


## ------------------- batch_summarise_ont_seq_reads

# Aim: aggregate the number of readcounts across all the batches
  # Input: 
  # sequencing_data = list: output of prepare_summary_file function after parsing through the ONT summary .txt file
# Pre-requisite:
  # Function: summarise_ont_seq_reads used to loop through each batch
# Output: df

batch_summarise_ont_seq_reads <- function(sequencing_data){
  
  dat <- list()
  for(i in 1:length(sequencing_data)){
    dat[[i]] = summarise_ont_seq_reads(sequencing_data[[i]], names(sequencing_data)[[i]])
  }
  dat <- bind_rows(dat)
  
  return(dat)
}


## ------------------- summarise_ont_demux_reads

# Aim: summarise the number of demultiplexed reads across the batches for downstream processing
# Input:
  # ReadCounts = list: output from read_batch_readcounts function after inputting <ReadCounts.txt> files
  # samples_removed = vector: samples to be removed for downstream processing
  # num_barcodes = int: number of barcodes (assuming 1:num_barcodes)
# Output:
  # dataframe of the merged demultiplexed reads across the batches

summarise_ont_demux_reads <- function(ReadCounts, samples_removed, num_barcodes){
  
  dat <- bind_rows(ReadCounts$Demultiplexed) %>% 
    # datawrangle for plots
    mutate(Sample = word(V2,c(1),sep = fixed("_")),
           Type = as.factor("Demultiplexed"),
           Strand = word(V2,c(2),sep = fixed("_")), 
           Batch = as.factor(Batch),
           BarcodedSample = paste0("Batch",Batch,Sample)) %>% 
    filter(Sample %ni% samples_removed) %>%
    dplyr::rename(File = V2, ReadCount = V1) %>%  
    mutate(Sample = factor(Sample %>% stringr::str_remove("BC"), levels = c(1:num_barcodes)))
  
  return(dat)
}


## ------------------- number_ont_reads

# Aim: wrapper function to read in the number of ONT reads at each processing stage for plotting
# Input:
  # sequencing_data = list: output of summary.txt file processed from prepare_summary_file in ont_runs_qc.R
  # Readcounts = list: output from read_batch_readcounts 
  # samples_removed = vector: list of samples to be removed for downstream processing
  # num_barcoes = int: number of barcodes (assuming 1..x)
  # ONTBarcodedPhenotype = df: samples with barcoded names etc
  # tg4510_samples = df: list of samples
# Output:
  # p1: line graph of number of reads through the ONT pipeline
  # p2: scatter plot of number of reads per barcode by batch and strand
  # p3: boxplot of number of reads per batch and by genotype
  # p4: boxplot of the number of demultiplexed reads between phenotype

number_ont_reads <- function(sequencing_data, ReadCounts, samples_removed, num_barcodes,ONTBarcodedPhenotype,tg4510_samples){
  
  # Summarise the number of sequencing counts in batch 2 and 3 (sequenced and basecalled)
  Sequencing_Counts <- batch_summarise_ont_seq_reads(sequencing_data)
 
  # Concatenate number of plus and minus reads across Batch 2 and Batch 3
  Demultiplexed_Counts <- summarise_ont_demux_reads(ReadCounts, samples_removed, num_barcodes)
  
  # Aggregate plus and minus reads for each sample, merge with phenotype information for barcode
  AllDemultiplexed_Counts = Demultiplexed_Counts %>% group_by(BarcodedSample) %>% tally(ReadCount) %>% merge(.,ONTBarcodedPhenotype) %>% 
    dplyr::rename(ReadCounts = n) %>% mutate(Batch = as.factor(Batch))
  
  # Basecalled Reads; Filtered Reads; Demultiplexed Reads 
  # Plot of stepwise drop in number of reads as processed through the ONT bioinformatics pipeline
  StepwiseReads = bind_rows(Demultiplexed_Counts[,c("ReadCount","Batch","Type")], Sequencing_Counts) %>% 
    mutate(Type = factor(Type, levels = c("Sequenced","Basecalled","Demultiplexed")),
           ReadCount = ReadCount/1000000)
  
  p1 = ggplot(StepwiseReads, aes(x = Type, y = ReadCount, colour = Batch, group = Batch)) + geom_line() + geom_point(size = 3) +
    mytheme + theme(legend.position = c(0.8,0.8)) + labs(x = "", y = "Number of Reads (Million)") +
    scale_colour_manual(name = "",labels = c("Batch 2 (n = 9) ","Batch 3 (n = 9)"), values = c("#00BFC4","#7CAE00")) +
    theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  cat("Mean Demultiplexed reads", mean(AllDemultiplexed_Counts$ReadCounts),"\n")
  cat("Min Demultiplexed reads", min(AllDemultiplexed_Counts$ReadCounts),"\n")
  cat("Max Demultiplexed reads", max(AllDemultiplexed_Counts$ReadCounts),"\n")
  
  # Scatter: Number of reads per barcode by batch and strand
  p2 = ggplot(Demultiplexed_Counts,aes(x = Sample, y = ReadCount, colour = Strand)) + geom_point(size = 2) + 
    facet_grid(~Batch,labeller=labeller(Batch = c("2" = "Batch 2", "3" = "Batch 3"))) + mytheme + 
    labs(x = "Barcode", y = "Number of Demultiplexed Reads (Million)") +
    theme(strip.background = element_blank(), legend.position = "bottom") +
    scale_colour_manual(values = c(label_colour("Plus"),label_colour("Minus")))
  
  # Boxplot: Number of reads per batch and by genotype
  p3 = ggplot(AllDemultiplexed_Counts ,aes(x = Batch, y = ReadCounts,fill = Phenotype)) + geom_boxplot() + 
    geom_point(aes(colour = Phenotype),position = position_jitterdodge(),size = 3) + labs(y = "Number of Demultiplexed Reads (Million)") +
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")))  + mytheme + 
    scale_fill_manual(values = c(alpha(label_colour("TG"),0.4),alpha(label_colour("WT"),0.4)))
  
  # Boxplot: Number of reads across genotype
  p4 = ggplot(AllDemultiplexed_Counts,aes(x = Phenotype, y = ReadCounts)) + geom_boxplot() + 
    geom_point(aes(colour = Phenotype), size = 3) +
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")))  + mytheme + 
    labs(y = "Number of Demultiplexed Reads (Million)") + theme(legend.position = "none")
  
  
  # Wilcoxin rank sum test of genotype against number of filtered reads
  #with(AllDemultiplexed_Counts, shapiro.test(ReadCounts[Phenotype == "WT"]))
  #with(AllDemultiplexed_Counts, shapiro.test(ReadCounts[Phenotype == "TG"]))
  #with(Reads_plot %>% filter(Description == "Transcripts"), shapiro.test(value[Phenotype == "TG"]))
  #var.test(value ~ Phenotype,Reads_plot %>% filter(Description == "Transcripts")) #cannot assume variance
  wilcox.test(ReadCounts ~ Phenotype,AllDemultiplexed_Counts) 
  
  # correlation of reads and RIN
  transcript_RIN <- merge(AllDemultiplexed_Counts,tg4510_samples, by.x = "sample",by.y = "Sample.ID", all.x = T)
  cor.test(transcript_RIN$ReadCounts,transcript_RIN$RIN, method = "spearman", exact = FALSE)
  
  # correlation of reads and Barcode
  cor.test(AllDemultiplexed_Counts$ReadCounts,as.numeric(AllDemultiplexed_Counts$Sample), method = "spearman", exact = FALSE)
  
  AllDemultiplexed_Counts %>% group_by(Phenotype, Batch) %>% tally()
  #wilcox.test(ReadCounts ~ Phenotype,AllDemultiplexed_Counts %>% filter(Batch == 2)) 
  #wilcox.test(ReadCounts ~ Phenotype,AllDemultiplexed_Counts %>% filter(Batch == 3)) 
  
  output = list(AllDemultiplexed_Counts,p1,p2,p3,p4)
  names(output) = c("ReadCounts","p1","p2","p3","p4")
  
  return(output)
}