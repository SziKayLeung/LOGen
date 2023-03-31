#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: generate basic descriptive numbers for target genes after FICLE annotation 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## summarise_gene_stats
## all_summarise_gene_stats
##
## ---------- Notes -----------------
##  
## Pre-requisite: Run FICLE on each target gene to generate (<TargetGene>_Final_Transcript_Classifications.csv)
## (<TargetGene>_Final_Transcript_Classifications.csv) records the number of splicing events across different categories
##


## ------------------- summarise_gene_stats

# Aim: summarise the number of transcripts etc for target gene 
# Pre-requisite: generate output files for target gene using FICLE
# Input:
  # class.files = df: input SQANTI classification file for number of transcripts
  # gene_class = df: output from annotation of target gene using FICLE <associated_gene> <Isoform> <Matching> <.....> 
    # each column records the number of transcripts with said alternative splicing event
  # input_cpat = df: output results from running CPAT for determining coding potential of isoform 
    # <seq_ID> <ID> .... <Coding_prob> <coding_stats>
    # note coding stats <Coding/Non_Coding> determined from the Coding_prob > threshold (0.44 for mouse)
  # input_noORF = df: output results from running CPAT of the isoforms with no open reading frame detected
    # <seq_ID> <coding_stats = No_ORF>
  # gene = str: target gene name 
# Output:
  # df: summary descriptives of the number of transcripts etc for said target gene

summarise_gene_stats <- function(class.files, gene_class, input_cpat, input_noORF, gene){
  
  # filter gene of interest
  gene_class_sqanti = class.files %>% filter(associated_gene == gene)
  
  # extract the coding potential related information from CPAT
  gene_class_sqanti = merge(gene_class_sqanti, rbind(input_cpat[,c("seq_ID","coding_status")], input_noORF[,c("seq_ID","coding_status")]), by.x = "isoform", by.y = "seq_ID")
  
  gene_class = gene_class %>% filter(Isoform %in% gene_class_sqanti$isoform)
  
  # generate output dataframe of the desired summary information
  dat = data.frame(
    "Total Number of Transcripts" = c(nrow(gene_class)),
    "Total Number of Unique Known Transcripts (ISM, FSM)" = c(length(unique(gene_class_sqanti[gene_class_sqanti$associated_transcript != "novel", "associated_transcript"]))),
    "Total Number of Known Transcripts" = c(nrow(gene_class_sqanti %>% filter(associated_transcript != "novel"))),
    "Total Number of Novel Transcripts" = c(nrow(gene_class_sqanti %>% filter(associated_transcript == "novel"))),
    "Total Number of Coding Transcripts" = c(nrow(gene_class_sqanti %>% filter(coding_status == "Coding"))),
    "Total Number of Transcripts with canonical splice sites" = nrow(gene_class_sqanti[gene_class_sqanti$all_canonical == "canonical",]),
    "Total Number of Transcripts with non canonical splice sites" = nrow(gene_class_sqanti[gene_class_sqanti$all_canonical == "non_canonical",]),
    "Total Number of Transcripts with non canonical splice sites and non-coding" = 
      nrow(gene_class_sqanti[gene_class_sqanti$all_canonical != "canonical" & gene_class_sqanti$coding_status == "Non_Coding",]),
    "Number of Transcripts with ES Events" = c(nrow(gene_class %>% filter(ES != 0))),
    "Number of Transcripts with IR Events" = c(nrow(gene_class %>% filter(IR != 0))),
    "Number of Transcripts with A5A3 Events" = c(nrow(gene_class %>% filter(A5A3 != 0))),
    "Number of Transcripts with alternative promoter" = c(nrow(gene_class %>% filter(AF != 0)))
  )
  
  dat = reshape2::melt(dat) 
  colnames(dat) = c("Description",gene)
  
  
  count = 1
  output = data.frame(Features = colnames(gene_class)[3:ncol(gene_class)])
  for(c in 3:ncol(gene_class)){
    output[count,2] = sum(gene_class[[c]])
    count = count + 1
  } 
  colnames(output) = c("Description",gene)
  
  final = rbind(dat,output)
  return(final)
}


## ------------------- all_summarise_gene_stats

# Aim: wrapper to run summarise_gene_stats across all target genes to generate one output table
# Input:
  # Gene_class = lst: stores output information from file (<TargetGene>_Final_Transcript_Classifications.csv) for each target gene after running FICLE 
  # class.files = df: input SQANTI classification file for number of transcripts
  # input_cpat = df: required for summarise_gene_stats
  # input_noORF = df: required for summarise_gene_stats
  # TargetGene = vec: all target gene names
# Output:
  # df: summarised descriptives of all target gene into one table (using information after running FICLE)

all_summarise_gene_stats <- function(Gene_class,class.files,input_cpat,input_noORF,TargetGene){

  # merged list of output of each target gene into one dataframe
  Gene_class_df = bind_rows(Gene_class, .id = "associated_gene") %>% dplyr::rename("Isoform" = "X") %>% subset(., select = -c(isoform))

  # run summarise_gene_stats for each target gene into a list
  TargetGene = sort(TargetGene)   
  Merged_gene_class = lapply(TargetGene, function(x) summarise_gene_stats(class.files,Gene_class_df,input_cpat,input_noORF,x))
  
  # combine list into dataframe
  Merged_gene_class_df = do.call(cbind,Map(cbind,Merged_gene_class))
  Merged_gene_class_df = Merged_gene_class_df  %>% column_to_rownames(., var = "Description") %>% dplyr::select(-contains("Des"))
  colnames(Merged_gene_class_df) = TargetGene
  
  return(Merged_gene_class_df)
}
