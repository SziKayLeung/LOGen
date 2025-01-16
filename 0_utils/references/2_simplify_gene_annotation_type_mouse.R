#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: To remove duplicated genes from annotation with multiple class
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## --------------------------------
## 
## Problem:
## Observed that a number of genes were from the gtf were classified under different types
## i.e. miRNAX --> classified under miRNA and lncRNA
## thus when merged with other files, creates duplicated entries
## 
## Solution: remove duplicated entries, keeping only class type for every gene
## Manually checked a few, but assumed that "lncRNA" and "pseudogene" are umbrella terms are not specific enough
## Therefore removed these classes from genes with multiple class type
##
## Input: txt file generated from create_gene_annotation_file.sh

## ---------- Packages -----------------

suppressMessages(library("tidyverse"))

## ----- Identify gene symbols with multiple class ---------------------

# read in annotation file
geneType <- read.table("C:/Users/sl693/Dropbox/gencode.vM22.annotation.geneannotation.txt", header = T)

# select on the GeneSymbol and Class and remove redudnant rows
geneType <- distinct(geneType %>% select(GeneSymbol, Class))
geneType <- geneType %>% mutate(GeneSymbolClass = paste0(GeneSymbol, "_", Class))

# Tally the number of times the GeneSymbol occurs, i.e. any genes that are in more than one class
numGeneType <- as.data.frame(table(geneType$GeneSymbol))

# Identify genes with more than one class
Redundant <- geneType[geneType$GeneSymbol %in% numGeneType[numGeneType$Freq > 1,"Var1"],]

# of the genes with more than one class, remove "lncRNA" and "unprocessed_pseudogene" and "TEC"
# after manually checking a few, lncRNA and unprocessed_pseudogene are umbrella terms given to any gene if class is not certain
RedCorrected <- Redundant[!Redundant$Class %in% c("lncRNA", "unprocessed_pseudogene","TEC"), ]

# Manually checked that these class identifications are wrong
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Mir1839_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Grik2_protein_coding")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Mir677_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Mir3068_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Mir1949_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Gm27475_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Gm24105_miRNA")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Zkscan7_protein_coding")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Ighv1-13_IG_V_pseudogene")
RedCorrected <- RedCorrected %>% filter(RedCorrected$GeneSymbolClass != "Ighv5-8_IG_LV_gene")

## ----- Checks -------------------------------

# Check 1: no duplicated gene symbols 
if(nrow(RedCorrected[duplicated(RedCorrected$GeneSymbol),]) == 0){
  message("Passed check 1 for geneType: no duplicated gene symbols")
}else{
  RedCorrected[RedCorrected$GeneSymbol %in% RedCorrected[duplicated(RedCorrected$GeneSymbol),"GeneSymbol"],]  
}

# Check 2: no missing gene symbols since possible that genes with classified as both lncRNA and unprocessed pseudogene
# and now removed
missing <- setdiff(unique(Redundant$GeneSymbol),unique(RedCorrected$GeneSymbol))
if(length(missing) > 0){
  message("Removed these genes as either coded as lncRNA or unprocessed pseudogene:")
  message(missing)
  
  for(i in missing){
    add <- data.frame(i, "lncRNA", paste0(i, "_lncRNA"))
    colnames(add) <- colnames(RedCorrected)
    RedCorrected <- rbind(RedCorrected, add)
  }
  
  reCheck <- setdiff(unique(Redundant$GeneSymbol),unique(RedCorrected$GeneSymbol))
  if(length(reCheck) > 0){
    message("Still missing even after replacing genes back as lncRNA")
  }else{
    message("Added missing genes as lncRNA")
    message("Passed check 2 for geneType: no missing gene symbols after replacement")
  }
}else{
  print("Passed check 2 for geneType: no missing gene symbols after replacement")
}

## ----- Final output file with no duplicated gene symbols -------------------------------

correctedGeneType <- geneType[!geneType$GeneSymbol %in% RedCorrected$GeneSymbol,]
correctedGeneType <- rbind(correctedGeneType, RedCorrected)

if(nrow(correctedGeneType[duplicated(correctedGeneType$GeneSymbol),]) == 0){
  message("Passed check 3 for geneType: no duplicated gene symbols in final list")
}

if(length(setdiff(geneType$GeneSymbol, correctedGeneType$GeneSymbol)) == 0){
  message("Passed check 4 for geneType: no missing gene symbols in final list")
}

write.csv(correctedGeneType, "gencode.vM22.annotation.geneannotation_corrected.csv", row.names = F)
