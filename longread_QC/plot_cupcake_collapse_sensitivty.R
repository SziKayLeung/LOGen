#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot sensitivity of cupcake collapse and determine threshold for filtering
## processed using SQANTI filtering and applied merge_class_abundance() to generate two columns:
## nsamples, nreads
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## plot_cupcake_collapse_sensitivity
##


# packages
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))


## ------------------- plot_cupcake_collapse_sensitivity

# Aim: scatter-plot of the frequency and percentage of isoforms with sample and read support
# Input:
  # df: sqanti class.files with nreads (total number of reads) and nsamples (tota number of samples) column
  # ntitle: str = title of plot 
# Output:
  # p1: frequency plot of the number of isoforms after filtering by the number of reads and samples
  # p2: percentage plot of the proportion of isoforms retained in the class.files after applying filters
# Notes:
  # nreads and nsamples filter is applied across all the samples
  # i.e. filter of 2 reads and 2 samples ==> isoforms with minimum 2 reads in any 2 samples and more

plot_cupcake_collapse_sensitivity <- function(df, ntitle){
  
  # set the threshold filter of nreads and nsamples
  nread_scale = c(seq(2,10),15,20)
  nsample_scale = c(2,3,5,8,10)
  
  # apply the threshold filter and record the number of isoforms retained
  sensitivity <- data.frame()
  count = 1
  for(i in nsample_scale){
    for(j in nread_scale){
      sensitivity[count,1] = i
      sensitivity[count,2] = j
      sensitivity[count,3] = nrow(df %>% filter(nsamples >= i & nreads >= j))
      count = count + 1  
    }
  }
  colnames(sensitivity) = c("num_sample","num_reads","totalisoremaining")
  
  # the total number of isoforms to calculate proportions
  totaliso = nrow(df)
  sensitivity = sensitivity %>% mutate(perc = totalisoremaining/totaliso)
  
  # plots
  p1 = ggplot(sensitivity, aes(x = num_reads, y = totalisoremaining, colour = as.factor(num_sample))) + geom_point() +
    labs(x = "Number of reads", y = "Number of isoforms retained", title = ntitle) + 
    scale_colour_discrete(name = "Number of samples") + theme_classic() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),labels = label_comma()) + mytheme 
  
  p2 = ggplot(sensitivity, aes(x = num_reads, y = perc, colour = as.factor(num_sample))) + geom_point() +
    labs(x = "Number of reads", y = "Percentage of isoforms retained", title = ntitle) + 
    scale_colour_discrete(name = "Number of samples") + theme_classic() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + mytheme + theme(legend.position = c(0.8,0.8))
  
  return(list(p1,p2))
}
