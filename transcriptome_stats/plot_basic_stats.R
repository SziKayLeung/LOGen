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

LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN, "aesthetics_basics_plots/pthemes.R"))
suppressMessages(library(scales))


## ------------------- total_num_iso

# Aim: plot the total number of isoforms per target gene in finalised dataset by structural category 
# Input: 
  # class.files = df: classification file 
  # input_title = str: title of the plot
  # dataset = str:"dataset" then split by dataset (if column exists in clasification file)
  # glimit = int: number of genes to display (top X ranked)

total_num_iso <- function(class.files, input_title, dataset = NA, glimit = NA){
  
  # structural category colours 
  cate_cols <- c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),alpha("#F8766D",0.8),alpha("#F8766D",0.3),alpha("#808080",0.3))
  
  # number of isoforms per gene
  nIso <- class.files %>% group_by(associated_gene) %>% tally %>% arrange(-n) %>% dplyr::rename("totaln" = "n")
  
  if(!is.na(dataset)){
    
      p <- class.files %>% group_by(associated_gene, dataset) %>% tally %>%
        ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = dataset)) +
        scale_fill_discrete(name = "Dataset")
    
  }else{
    
    if(!is.na(glimit)){
      
      p <- class.files %>% group_by(associated_gene, structural_category) %>% tally %>% 
        full_join(., nIso, by = "associated_gene") %>%
        filter(associated_gene %in% nIso$associated_gene[1:glimit]) %>%
        ggplot(.,aes(x = reorder(associated_gene, totaln), y = n, fill = structural_category)) +
        scale_fill_manual(name = "Isoform Classification", values = cate_cols)
      
    }else{
      
      p <- class.files %>% group_by(associated_gene, structural_category) %>% tally %>%
        ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = structural_category)) +
        scale_fill_manual(name = "Isoform Classification", values = cate_cols)
      
    }
    
  }
  
  p <- p + geom_bar(stat = "identity") + mytheme + labs(x = "", y = "Number of isoforms", title = input_title) + mytheme +
    theme(legend.position = c(0.8,0.8)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),labels = label_comma())
  
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


## ------------------- tabulate/plot_structural_cate

# Aim: tabulate the total number of isoforms in finalised dataset by structural category 
# Input:  
  # class.files = df: SQANTI classification file

tabulate_structural_cate <- function(class.files){
  
  # group by structural_category and tally with proportion
  output <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100)
  
  return(output)
}


# Aim: plot the total number of isoforms in finalised dataset by structural category 
# Input:  
  # class.files = df: SQANTI classification file
# Output:
  # p = bar-plot of the number of isoforms by structural category (rotated for paper purposes)

plot_structural_cate <- function(class.files){
  
  # colours
  structuralcolours = data.frame(
    "FSM" = rgb(249,176,161,maxColorValue = 255),
    "ISM" = rgb(228,196,116,maxColorValue = 255),
    "NIC" = rgb(204,220,151,maxColorValue = 255),
    "NNC" = rgb(108,220,164,maxColorValue = 255),
    "Fusion" = rgb(114, 215, 217,maxColorValue = 255),
    "Genic_Genomic" = rgb(226, 226, 243,maxColorValue = 255),
    "Genic_Intron" = rgb(226, 226, 243,maxColorValue = 255),
    "Antisense" = rgb(108,204,251,maxColorValue = 255),
    "Intergenic"= rgb(252,140,220,maxColorValue = 255)
  ) %>% gather(.,category,colours,factor_key=TRUE)
  
  
  p <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100) %>% 
    mutate(structural_category = factor(structural_category, levels = structuralcolours$category)) %>%
    ggplot(., aes(x = structural_category, y = perc, fill = structural_category)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural category", y = "Percentage of isoforms (%)") +
    #scale_fill_discrete(values = c("#F9B0A1", "#E4C474", "#CCDC97", "#6CDCA4", "#72D7D9", "#E2E2F3", "#E2E2F3", "#6CCCFB", "#FC8CDC")) +
    coord_flip() +
    scale_x_discrete(limits = rev(structuralcolours$category)) +
    scale_fill_viridis(discrete=TRUE) +
    theme(legend.position = "None")
  
  return(p)
}
