#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Power calculations  
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## cumulative_added_sample
## plot_cumulative_added_sample
##   
## ---------- Notes -----------------
## 
## Plot the cumulative number of transcripts with added sample 



## ---------- Packages -----------------

suppressMessages(library(RColorBrewer))


## ------------------- cumulative_added_sample

# Aim: generate a table of the cumulative number of transcripts with every added sample by structural_category
# Input: 
  # class.files = df: classification file generated from SQANTI
  # demux = df: demux_fl_count.csv with the first column as "id", and the FL reads of other samples
  # structural_category = <FSM> <NIC> <NNC>
# Output: 
  # table of the cumulative number

cumulative_added_sample <- function(class_file, demux, structural_category){
  
  # subset the class file by structural_category
  df <- class_file[class_file$structural_category == structural_category,]
  
  # select only the number of columns with expression (FL read count)
  df <- df %>% select(isoform, colnames(demux)[-1])
  
  # loop through each column <isoform> <sample 1> <sample 2>
  dfOutput <- data.frame()
  for(i in 2:ncol(df)){
    
    dfOutput[i - 1, 1] <- as.numeric(paste0(i - 1))
    
    # isoform ID that are detected in that specific sample (i.e. expression > 0 FL read count)
    isoformSample <- df[df[[i]] > 0,"isoform"]
    
    if(i == 2){
      # number of isoforms as base in first sample
      dfOutput[i - 1, 2] <- length(isoformSample)
      
      # set that as base isoforms
      common <- isoformSample
    }else{
      
      # for the remaining samples
      # tabulate what is already known, plus additional transcripts
      # update common to include those
      common <- c(common, setdiff(isoformSample,common))
      dfOutput[i - 1, 2] <- length(common) 
    }
    
  }
  
  colnames(dfOutput) <- c("Sample", "number_of_transcripts")
  dfOutput <- dfOutput %>% mutate(structural_category = structural_category)
  return(dfOutput)
}

## ------------------- plot_cumulative_added_sample

plot_cumulative_added_sample <- function(class_file, demux){
  
  cumuAll <- lapply(unique(class_file$structural_category), function(x) cumulative_added_sample(class_file, demux, x))
  
  p <- ggplot(bind_rows(cumuAll), aes(x = Sample, y = number_of_transcripts, colour = structural_category)) + 
    geom_point() +
    labs(x = "Number of samples", y = "Cumulative number of transcripts") + 
    geom_smooth(method = "loess", se = FALSE) +
    scale_colour_brewer(palette = "Set1") +
    theme_classic()
  
  return(p)
}


