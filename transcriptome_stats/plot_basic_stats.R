#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose:
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## final_num_iso
## 
##   
## ---------- Notes -----------------
## 


## ------------------- total_num_iso

# Aim: plot the total number of isoforms in finalised dataset by structural category 
# Input: 
  # class.files = df: classification file with only the associated_genes of interest kept 
  # input_title = str: title of the plot


total_num_iso <- function(class.files, input_title){
  
  p <- class.files %>% mutate(structural_category = factor(structural_category, levels = c("FSM","ISM","NIC","NNC","Genic\nGenomic"))) %>%
    group_by(associated_gene, structural_category) %>% tally %>%
    ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = structural_category)) + geom_bar(stat = "identity") +      
    mytheme + labs(x = "", y = "Number of Isoforms", title = input_title) + mytheme +
    scale_fill_manual(name = "Isoform Classification", values = c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                  alpha("#F8766D",0.8),alpha("#F8766D",0.3),
                                                                  alpha("#808080",0.3))) +
    theme(legend.position = c(0.8,0.8)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}


## ------------------- length_exon_description

# Aim: generate two statements about mean, sd and range of isoform length and number of exons
# for paper reporting
# Input:  
  # class.files = df: SQANTI classification file
length_exon_description <- function(class.files){
  print(paste0("Length: mean = ", signif(mean(class.files$length),3),
               " s.d = ", signif(sd(class.files$length),3),
               " range: ", signif(min(class.files$length),3),"-", signif(max(class.files$length),3)))
  
  print(paste0("Exon number: median = ", signif(median(class.files$exons),3),
               " s.d = ", signif(sd(class.files$exons),3),
               " range: ", signif(min(class.files$exons),3),"-", signif(max(class.files$exons),3)))
}


tabulate_structural_cate <- function(class.files){
  
  # group by structural_category and tally with proportion
  output <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100)
  
  return(output)
}