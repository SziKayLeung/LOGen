#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot histograms of cage, TSS and TTS peaks across datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## prepare_hist_breaks
## hist_TSS_Cage_TTS_peak_2datasets
## wrapper_hist_peaks_2datasets 
## plot_polyA_motifs
##


## ---------- Packages -----------------

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


## ------------------- prepare_hist_breaks

# Aim: data-wrangle class.files to generate columns for downstream plot
# Input: 
  # class.files = df: SQANTI classification files
  # feature = str: column name of class.files interested for plotting
# Output:
  # class.files with diff2 column for plotting

prepare_hist_breaks <- function(class.files, feature){
  class.files <- class.files %>% .[!is.na(.[[feature]]),] # remove NAs for that specific feature 
  
  # settting the breaks for the histogram
  # threshold is for replacing >200, whic differs for each feature and class.filesaset
  diff_max <- max(max(abs(class.files[[feature]])), max(abs(class.files[[feature]])))
  diff_breaks <- c(-(diff_max+1), seq(-200, 200, by = 20), diff_max+1)
  class.files$diff <- cut(-(class.files[[feature]]), breaks = diff_breaks) 
  threshold <- paste(formatC(diff_breaks[1], format = "g", digits = 3))
  
  # formatting of the x-axis to grab the first number, between the brackets
  class.files$diff2 <- gsub("]", "", paste(class.files$diff))
  class.files$diff2 <- word(class.files$diff2,c(2),  sep = fixed ('('))
  class.files$diff2 <- word(class.files$diff2,c(1),  sep = fixed (','))
  class.files$diff2 <- as.factor(class.files$diff2)
  
  levels(class.files$diff2)[match(threshold,levels(class.files$diff2))] <- ">-220"
  class.files$diff2 <- factor(class.files$diff2,
                              levels = c(">-220","-200","-180","-160","-140","-120",
                                         "-100","-80","-60","-40","-20","0",
                                         "20","40", "60", "80","100",
                                         "120","140","160","180","200" ))
  return(class.files)

}


## ------------------- hist_TSS_Cage_TTS_peak_2datasets

# Aim: plot the histogram of two datasets based on the feature (SQANTI column) of interest
# Input:
  # class.files.1 = df: classification file of dataset 1
  # class.files.2 = df: classification file of dataset 2
  # name1 = str: name of dataset 1 (ensure same name as label_colour)
  # name2 = str: name of dataset 2 (ensure same name as label_colour)
  # feature = str: column name of class.files interested for plotting
# Output:
  # p1 = plot of the histogram of the two datasets of percentage of isoforms over the feature of interest 

hist_TSS_Cage_TTS_peak_2datasets <- function(class.files.1, class.files.2, name1, name2, feature){
  
  # plot 
  x_label <- 
    if(feature == "diff_to_TSS"){paste("Distance to Annotated TSS (bp)")}else{
      if(feature == "dist_to_cage_peak"){paste("Distance to Annotated CAGE Peak (bp)")}else{
        if(feature == "diff_to_TTS"){paste("Distance to Annotated TTS (bp)")}else{
          if(feature == "diff_to_gene_TTS"){paste("Distance to any Annotated TTS (bp)")}else{
            if(feature == "diff_to_gene_TSS"){paste("Distance to any Annotated TSS (bp)")}else{
              paste("NA")
            }
          }
        }
      }
    }
  
  
  df1 <- prepare_hist_breaks(class.files.1, feature)
  df2 <- prepare_hist_breaks(class.files.2, feature)
  
  merge <- rbind(df1[,c("diff2","Dataset")], df2[,c("diff2","Dataset")]) 
  total = merge %>% group_by(Dataset) %>% tally()
  p <- merge(merge %>% group_by(diff2, Dataset) %>% tally(), total, by = "Dataset") %>%
    mutate(perc = n.x / n.y * 100) %>% 
    ggplot(., aes(x=diff2, fill = Dataset, y = perc)) +
    geom_bar(stat = "identity", color="black", linewidth=0.3, 
             position = position_dodge())+
    mytheme +
    scale_fill_manual(values = c(label_colour(name1),label_colour(name2)), labels = c(name1,name2)) +
    labs(x = x_label, y = "Isoforms (%)") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  
  return(p)
  
}


## ------------------- wrapper_hist_peaks_2datasets 

# Aim: wrapper to plot the histogram of the cage peak, TSS and TTS of 2 datasets side-by-side
# Input:
  # class.files.1 = df: classification file of dataset 1
  # class.files.2 = df: classification file of dataset 2
  # name1 = name of dataset 1 (ensure same name as label_colour)
  # name2 = name of dataset 2 (ensure same name as label_colour)
# Output: plot of the histogram of the two datasets of percentage of isoforms with distribution of
  # p1 = cage peaks
  # p2 = TSS 
  # p3 = TTS

wrapper_hist_peaks_2datasets <- function(class.files.1, class.files.2, name1, name2){

  cage_plots_merged <- hist_TSS_Cage_TTS_peak_2datasets(class.files.1, class.files.2, name1, name2, "dist_to_cage_peak")
  tss_plots_merged <- hist_TSS_Cage_TTS_peak_2datasets(class.files.1, class.files.2, name1, name2, "diff_to_gene_TSS")
  tts_plots_merged <- hist_TSS_Cage_TTS_peak_2datasets(class.files.1, class.files.2, name1, name2, "diff_to_TTS")
  
  
  return(list(cage_plots_merged, tss_plots_merged, tts_plots_merged))
}


## ------------------- plot_polyA_motifs

# Aim: plot the number of isoforms with different poly-A motifs (x-variable) across different datasets
# Input:
  # class.files = df: SQANTI classification file with "Dataset" column as identifer 
# Output:
  # p1 = bar-plot of the number of isoforms (log10) vs polyA motifs

plot_polyA_motifs <- function(class.files){
  
  p <- class.files %>% group_by(polyA_motif, Dataset) %>% tally() %>% mutate(log_n = log10(n)) %>% filter(!is.na(polyA_motif)) %>% 
    ggplot(., aes(x = reorder(polyA_motif, -log_n), y = log_n, fill = Dataset)) + geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "PolyA motif", y = "Number of Isoforms (log10)") + mytheme
  
  return(p)
}