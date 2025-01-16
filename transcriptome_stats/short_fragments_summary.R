
#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## subset_class_by_length()
## summarise_subsetted_classLength
## statsReadsByLength()
## ---------------------------------


## ------------------- subset_class_by_length

# Aim: subset unfiltered SQANTI classification transcripts by length 
# Input:
  # lengthCutOff = integer for length threshold i.e. if 200, all transcripts less than 200bp retained 
  # preFilteredClassFile = dataframe of SQANTI classification file before filtering i.e. isoforms and artifacts retained
  # demux = dataframe of expression counts of all isoforms pre-filtered (ensure column with nreads = total number of reads across all samples)
  # ReasonFiltered = dataframe of SQANTI txt of rationale of isoforms discarded during SQANTI filtering
  # geneType = reference annotation file (i.e. gencode.v40.annotation.geneannotation.txt for human); check LOGen/0_utils/references/ 
# Output:
  # df: subsetted classification file by length, and merged with demux, SQANTI reason file and geneType

subset_class_by_length <- function(lengthCutOff, preFilteredClassFile, demux, ReasonFiltered, geneType){
  
  # Convert data frames to data.table for faster processing
  preFilteredClassFile <- as.data.table(preFilteredClassFile)
  demux <- as.data.table(demux, keep.rownames = TRUE) # Keep row names as a column for merging
  ReasonFiltered <- as.data.table(ReasonFiltered)
  geneType <- as.data.table(geneType)
  
  # Add structural_category_exons column
  preFilteredClassFile[, structural_category_exons := paste0(structural_category, "_", exons)]
  
  # Filter rows by length cutoff
  shortdf <- preFilteredClassFile[length <= lengthCutOff, .(
    isoform, chrom, strand, length, exons, structural_category, 
    associated_gene, associated_transcript, filter_result, structural_category_exons
  )]
  
  # Merge with demux based on isoform (row names in demux are now 'rn')
  shortdf <- merge(shortdf, demux, by.x = "isoform", by.y = "rn", all.x = TRUE)
  
  # Merge with geneType for the Class column (i.e the type of RNA)
  shortdf <- merge(shortdf, unique(geneType[, .(GeneSymbol, Class)]), 
                   by.x = "associated_gene", by.y = "GeneSymbol", all.x = TRUE)
  
  # Filter ReasonFiltered by isoform and merge
  filteredlengthCutOff <- ReasonFiltered[isoform %in% shortdf$isoform]
  shortdf <- merge(shortdf, filteredlengthCutOff[, .(isoform, reasons)], 
                   by = "isoform", all.x = TRUE)
  
  return(shortdf)
}


## ------------------- summarise_subsetted_classLength

# Aim: table of summary information of shortdf from subset_class_by_length() at the transcript level
# Input:
  # shortdf = output from subset_class_by_length()
  # preFilteredClassFile = dataframe of SQANTI classification file before filtering i.e. isoforms and artifacts retained
# Output:
  # list of summary stats that can be converted to a dataframe downstream

summarise_subsetted_classLength <- function(shortdf, preFilteredClassFile){
  
  stats <- list()
  
  # number of transcripts before SQANTI filtering
  stats[["All Transcripts before SQANTI filtering"]] <- nrow(preFilteredClassFile)
  
  # number of short transcripts and percentage of all transcripts
  stats[["Short Transcripts before SQANTI filtering"]] <- paste0(nrow(shortdf), " (", round(nrow(shortdf)/nrow(preFilteredClassFile) * 100,2), "%)")
  
  # number of transcripts kept from SQANTI filtering
  shortTRUETranscripts <- nrow(shortdf %>% filter(filter_result == "Isoform")) 
  stats[["Short Transcripts after SQANTI filtering"]] <- paste0(shortTRUETranscripts, " (", round(shortTRUETranscripts/nrow(preFilteredClassFile) * 100,2), "%)")
  
  # number of transcripts kept from SQANTI and mono-exonic filtering
  shortdfSQANTIKept <- shortdf %>% 
    filter(filter_result == "Isoform") %>% 
    filter(!structural_category_exons %in% c("Intergenic_1","Genic_Intron_1","Genic_Genomic_1"))
  
  stats[["Short transcripts after SQANTI and mono-exonic filtering"]] <- 
    paste0(nrow(shortdfSQANTIKept), " (", round(nrow(shortdfSQANTIKept)/nrow(preFilteredClassFile) * 100,2), "%)")
  
  expressionFilter <- shortdfSQANTIKept[shortdfSQANTIKept$nreads >= 2 & shortdfSQANTIKept$nsample >= 2,]
  stats[["Short transcripts after SQANTI, mono-exonic and expression filtering"]] <- 
    paste0(nrow(expressionFilter), " (", round(nrow(expressionFilter)/nrow(preFilteredClassFile) * 100,2), "%)")
  
  return(stats)
}


## ------------------- summarise_subsetted_classLength_read

# Aim: table of summary information of shortdf from subset_class_by_length() at the read level
# Input:
  # shortdf = output from subset_class_by_length()
  # demux = dataframe of expression counts of all isoforms pre-filtered (ensure column with nreads = total number of reads across all samples)
# Output:
  # list of summary stats that can be converted to a dataframe downstream

summarise_subsetted_classLength_read <- function(shortdf, demux){
  
  stats <- list()
  
  stats[["Reads across dataset (all Transcripts)"]] <- sum(demux$nreads)
  
  stats[["Reads for short transcripts before SQANTI filtering"]] <- paste0(sum(shortdf$nreads), " (", round(sum(shortdf$nreads)/sum(demux$nreads) * 100,2), "%)")
  
  final <- shortdf %>% 
    filter(filter_result == "Isoform") %>% 
    filter(!structural_category_exons %in% c("Intergenic_1","Genic_Intron_1","Genic_Genomic_1")) %>%
    filter(nreads >= 2 & nsample >= 2)
  
  stats[["Reads for short transcripts retained after all filtering criteria"]] <- 
    paste0(sum(final$nreads), " (", round(sum(final$nreads)/sum(demux$nreads) * 100,2), "%)")
  
  return(stats)
}


## ------------------- subset_class_length_SQFiltered

# Aim: return dataframe of the transcripts filtered by SQANTI and rationale 
# Input;
  # shortdf = output from subset_class_by_length()
# Output:
  # result_df = count of the number of transcripts classified as "artifact" from SQANTI and the number of reads

subset_class_length_SQFiltered <- function(shortdf){
  
  result_df <- shortdf %>% filter(filter_result == "Artifact") %>% 
    group_by(structural_category, reasons) %>% 
    tally(nreads, name = "total_nreads") %>% 
    as.data.frame() %>% 
    mutate(perc = (total_nreads/sum(total_nreads))  * 100)  
  
  return(result_df)
}


## ------------------- final_kept_shortDataset

# Aim: final dataset of short fragments after applying multiple filtering criteria
  # SQANTI filter of isoform
  # remove monoexonic isoforms classifed as intergenic, genic_intron, genic_genomic
  # expression threshold of minimum 2 reads across any 2 samples
# Input;
  # shortdf = output from subset_class_by_length()
# Output:
  # result_df = count the number of short transcripts retained after applying multiple filtering criteria

final_kept_shortDataset <- function(shortdf){
  
  result_df <- shortdf %>% 
    filter(filter_result == "Isoform") %>% 
    filter(!structural_category_exons %in% c("Intergenic_1","Genic_Intron_1","Genic_Genomic_1")) %>%
    filter(nreads >= 2 & nsample >= 2) %>% 
    mutate(monoExonic = ifelse(exons == 1, "monoExonic", "multiExonic")) %>% 
    group_by(structural_category, Class, monoExonic) %>% tally()   
  
  return(result_df)
}