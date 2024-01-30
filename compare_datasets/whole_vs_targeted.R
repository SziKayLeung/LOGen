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
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/draw_venn.R"))
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
  # p3: venn diagram of the number of isoforms unique and common across targeted and whole dataset
  # p4: density plot of the number of FL reads to isoforms unique to targeted and whole dataset
  # p5: bar-plot of the number of isoforms unique to whole dataset by number of exons and structural category
  # p6: bar-plot of the number of reads of isoforms unique to whole dataset with more than 40 exons
  # p7: bar-plot of the number of isoforms unique to whole dataset across all target genes, coloured by structural category
  # p8: histogram of the sum of FL reads across all samples for isoforms unique to whole dataset
  # matchedSumTargeted: list of the isoforms and sum FL read expression across targeted and whole transcriptome dataset
  # geneRecord: number of isoforms unique/common to targeted vs whole across all target genes 

whole_vs_targeted_plots <- function(classfiles, wholeSamples, targetedSamples, targetGene){
  
  # subset by target genes
  classfiles <- classfiles %>% filter(associated_gene %in% targetGene)
  message("Number of samples in whole dataset: ",  length(wholeSamples))
  message("Number of samples in targeted dataset : ",  length(targetedSamples))
  message("Num of genes: ", length(unique(classfiles$associated_gene)))
  message("Target genes not detected in datasets: ")
  setdiff(targetGene, unique(classfiles$associated_gene))
  message("Num of transcripts: ", length(unique(classfiles$isoform)))
  
  # sum the number of FL reads across the targeted and whole samples
  matchedSum = cbind(data.frame(rowSums(classfiles[targetedSamples])),
                     data.frame(rowSums(classfiles[wholeSamples])))
  colnames(matchedSum) <- c("sumTargeted","sumWhole")
  matchedSum$isoform <- classfiles$isoform
  
  # annotate the isoforms and subset by target genes
  matchedSumTargeted = merge(matchedSum, classfiles [,c("isoform","associated_gene","structural_category")])
  
  # create a column by determining if the isoform is detected in both, whole or targeted
  # if FL read >= 1; then considered detected
  matchedSumTargeted$dataset <- apply(matchedSumTargeted, 1, function(x) identify_dataset_by_counts (x[["sumTargeted"]], x[["sumWhole"]], "Targeted","Whole"))
  
  # list of isoforms by dataset and commonality
  listIso <- list(
    # Targeted (unique + common)
    targeted = matchedSumTargeted[matchedSumTargeted$dataset != "Whole","isoform"],
    # Whole (unique + common)
    whole = matchedSumTargeted[matchedSumTargeted$dataset != "Targeted","isoform"], 
    # Whole Unique
    wholeUnique = matchedSumTargeted[matchedSumTargeted$dataset == "Whole","isoform"]
  )
  
  ## dataframe of the number of isoforms unique, common to whole and targeted dataset 
  geneRecord <- matchedSumTargeted %>% group_by(associated_gene, dataset) %>% tally %>% 
    filter(dataset != "NA") %>%
    spread(., dataset, n) 
  geneRecord[is.na(geneRecord)] <- 0
  geneRecord$Whole <- as.numeric(as.character(geneRecord$Whole))
  # create a column if more transcripts detected in whole dataset vs targeted dataset
  geneRecord <- geneRecord %>% mutate(WholeMore = ifelse(Whole < Targeted, FALSE,TRUE)) %>% 
    mutate(diffWholeTargeted = Whole - Targeted)
  message("Number of genes with more unique transcripts detected in whole than targeted dataset: ",
          nrow(geneRecord %>% filter(WholeMore == TRUE)))
  
  # plots
  totaln = matchedSumTargeted %>% group_by(associated_gene) %>% tally 
  p1 <- matchedSumTargeted %>% group_by(associated_gene, dataset) %>% tally %>% 
    full_join(., totaln, by = "associated_gene") %>% filter(dataset != "NA") %>%
    ggplot(., aes(x = reorder(associated_gene,-n.y), y = n.x, fill = dataset)) + geom_bar(stat = "identity") +
    mytheme + labs(x = "Target Genes", y = "Number of isoforms") +
    scale_fill_manual(name = "", values = c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling2")[1],wes_palette("Darjeeling1")[2])) + 
    theme(legend.position = "top") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(length(TargetGene) > 30){
  p1 <- p1 + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  p2 <- matchedSumTargeted %>% filter(dataset != "Both") %>%
    group_by(structural_category, dataset) %>% tally() %>% 
    filter(dataset != "NA") %>%
    ggplot(., aes(x = structural_category, y = n, fill = dataset)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural Category", y = "Number of isoforms") +
    scale_fill_manual(name = "", values = c(wes_palette("Darjeeling2")[1],wes_palette("Darjeeling1")[2])) + 
    theme(legend.position = "top")

  p3 <- twovenndiagrams(listIso$targeted, listIso$whole, "Targeted","Whole")
  
  
  p4 <- matchedSumTargeted %>% filter(dataset %in% c("Targeted","Whole")) %>% select(sumTargeted, sumWhole, dataset) %>% 
    reshape2::melt(id = "dataset") %>% 
    ggplot(., aes(x = log10(value), fill = dataset)) + geom_density(alpha = 0.2) + 
    mytheme + labs(x = "log10(FL count)", y = "Density") +
    scale_fill_manual(name = "", values = c(wes_palette("Darjeeling2")[1],wes_palette("Darjeeling1")[2])) +
    theme(legend.position = "top")
  
  message("Number of isoforms unique to whole dataset and mono-exonic: ",nrow(classfiles %>% filter(isoform %in% listIso$wholeUnique) %>% filter(exons ==1)))
  p5 <- classfiles %>% filter(isoform %in% listIso$wholeUnique) %>% group_by(structural_category, exons) %>% tally() %>%
    ggplot(., aes(x = exons, y = n, fill = structural_category)) + geom_bar(stat = "identity") + mytheme +
    labs(x = "Number of exons", y = "Number of isoforms", fill = "Structural category") +
    theme(legend.position = c(0.8,0.8))
  
  p6 <- classfiles %>% filter(isoform %in% listIso$wholeUnique) %>% filter(exons > 40) %>% 
    ggplot(., aes(x = isoform, y = nreads, fill = structural_category)) + geom_bar(stat = "identity") + 
    facet_grid(~associated_gene,scales = "free", space='free') + mytheme + theme(legend.position = "top") +
    labs(x = "Transcript", y = "Number of FL reads") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p7 <- classfiles[classfiles$isoform %in% listIso$wholeUnique,] %>% 
    select("isoform","associated_gene","structural_category") %>% 
    group_by(associated_gene, structural_category) %>% tally() %>% 
    ggplot(., aes(y = reorder(associated_gene,n), x = n, fill = structural_category)) + geom_bar(stat = "identity") + theme_classic() + labs(x = "Number of isoforms unique to whole dataset", y = "Genes")
  
  message("Number of isoforms with mean FL read < 10s: ", 
            nrow(classfiles %>% filter(isoform %in% listIso$wholeUnique) %>% 
         mutate(meanreads = nreads/length(wholeSamples)) %>% filter(meanreads < 10)))
  
  p8 <- classfiles[classfiles$isoform %in% listIso$wholeUnique,] %>% 
    select(isoform, contains("Whole")) %>% 
    reshape2::melt(id = "isoform") %>% 
    group_by(isoform) %>%
    dplyr::summarize(mean = mean(value, na.rm=TRUE)) %>% 
    ggplot(., aes(x = mean)) + geom_histogram() + scale_y_log10() + 
    theme_classic() + labs(x = "Mean FL reads across all samples", y = "Number of isoforms")
  
  # total number of FL reads across all isoforms detected across each whole transriptome sample
  allwholetotalFL <- classfiles %>% 
    select(isoform, contains("Whole")) %>% 
    reshape2::melt(id = "isoform") %>% group_by(variable) %>%
    dplyr::summarize(sum = sum(value, na.rm=TRUE)) %>%
    `colnames<-`(c("sample","totalFL"))
  
  # number of FL reads across isoforms unique to whole dataset across each whole transcriptome sample 
  uniquewholetotalFL <- classfiles[classfiles$isoform %in% listIso$wholeUnique,] %>% 
    select(isoform, contains("Whole")) %>% 
    reshape2::melt(id = "isoform") %>% group_by(variable) %>% 
    dplyr::summarize(sum = sum(value, na.rm=TRUE)) %>%
    `colnames<-`(c("sample","uniqueIsoFL"))
  
  merge(allwholetotalFL, uniquewholetotalFL, by = "sample") %>% 
    mutate(perc = uniqueIsoFL/totalFL * 100) %>% 
    ggplot(aes(x = sample, y = perc)) + geom_bar(stat = "identity") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
    labs(y = "FL reads of isoforms unique to whole dataset (%)", x = "Sample")
  
  # write txt file of lists of unique isoforms
  message("Output list of unique isoforms to whole dataset:", paste0(dirnames$output,"uniqueWholeIso.txt"))
  write.table(listIso$wholeUnique,paste0(dirnames$output,"uniqueWholeIso.txt"),quote=F,row.names = F, col.names = F)

  return(list(p1,p2,p3,p4,p5,p6,p7,p8,matchedSumTargeted,geneRecord))
}

