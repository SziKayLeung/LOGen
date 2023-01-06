#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: plot number/proportion of isoforms within specific feature and test for differences between two datasets (fisher test)
## feature i.e within 50bp of cage peak, within 50bp TTS, within a polyA site
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## plot_features_byisocate
## plot_stats_feature_2datasets



## ---------- Packages -----------------

suppressMessages(library("tibble"))
suppressMessages(library("tidyr"))


## ------------------- plot_features_byisocate

# Aim: bar-plot the proportion of isoforms within features (cage peak, TSS within 50bp), split by structural category
# Input: 
  # class.files = df: SQANTI classification file 
  # feature = str: column of SQANTI classification file of interest for plotting <within_polya_site, within_cage_peak, within_50_TTS>
  # featurename = str: title of the legend
# Pre-requisite: 
  # class.files previously annotated using annotate_class_binary() to generate columns of within_50_features etc
# Output: 
  # p: percentage of isoforms within features vs sturcutral category

plot_features_byisocate <- function(class.files, feature, featurename){
  
  # subset class.files by feature and tally the number of isoforms by structural category
  dat = class.files %>% 
    group_by_(feature, "Dataset", "structural_category") %>% tally() %>%
    left_join(., class.files %>% group_by(Dataset) %>% tally(), by = "Dataset") %>% 
    left_join(., class.files %>% group_by(Dataset, structural_category) %>% tally(), 
              by = c("Dataset","structural_category")) %>%
    mutate(perc = n.x/n.y * 100, perc_lab = paste0(round(n.x/n * 100,2),"%")) %>%
    mutate(structural_category = factor(structural_category, levels = c("NNC","NIC","ISM","FSM")))
  
  if(feature == "within_50_cage"){dat = dat %>%  mutate(perc_lab = ifelse(!grepl("Not",within_50_cage),perc_lab,NA)) }else 
    if (feature == "within_polya_site"){dat = dat %>%  mutate(perc_lab = ifelse(within_polya_site == "Yes",perc_lab,NA))}else
      if (feature == "within_50_TTS"){dat = dat %>%  mutate(perc_lab = ifelse(!grepl("Not",within_50_TTS),perc_lab,NA))}
  
  # plot 
  p = ggplot(dat, aes(x = perc, y = structural_category, fill = forcats::fct_rev(!! enquo(feature)), label = perc_lab)) + 
    geom_bar(stat = "identity") + facet_grid(~Dataset) +
    geom_text(size = 3, hjust = 0.25, colour = wes_palette("Cavalcanti1")[2]) + labs(y = "", x = "Isoforms (%)") + mytheme + 
    theme(legend.position = "top") + guides(fill=guide_legend(title=featurename)) +
    scale_fill_manual(values = c(label_colour("Yes"),label_colour("No")))
  
  return(p)
}


## ------------------- plot_stats_feature_2datasets

# Aim: Plot and run a fisher test on two datasets to test for difference in percentage of isoforms within and not within feature
# Input:
  # class.files = df: SQANTI classification file of merged datasets with "Dataset" column as identifier
  # feature = str: column of SQANTI classification file of interest for plotting <within_polya_site, within_50_cage, within_50_TTS>
  # featurename = str: title of the legend
# Pre-requisite:
  # class.files previously annotated using annotate_class_binary() to generate columns of within_50_features etc
  # class.files of merged > 2 dataset, with "Dataset" column as identifer
# Output:
  # printed fisher.test 
  # p: bar-plot of the two datasets vs percentage of isoforms within and not within feature

plot_stats_feature_2datasets <- function(class.files, feature, featurename){
  
  # remove NA, group by feature and dataset and tally
  df_stat = class.files %>% filter(!is.na(.data[[feature]])) %>%
    group_by_("Dataset",feature) %>% tally() %>% spread(., Dataset, n)  %>%  
    column_to_rownames(., var = feature)
  
  # perform test
  cat("Stats test for:", feature,"\n")
  cat("****** Table for fisher test ******** \n")
  print(df_stat)
  res = fisher.test(df_stat)
  print(res)
  cat("Pvalue: ", res$p.value, "\n")
  
  # generate bar-plot
  p = df_stat %>% mutate_if(is.numeric, funs(./sum(.) * 100)) %>% rownames_to_column("within50") %>% reshape2::melt() %>% 
    ggplot(., aes(x = variable, y = value, fill = forcats::fct_rev(within50))) + geom_bar(stat = "identity") + mytheme + 
    labs(x = "Dataset", y = "Isoforms (%)") + theme(legend.position = "bottom") +
    guides(fill=guide_legend(title=featurename)) +
    scale_fill_manual(values = c(label_colour("Yes"),label_colour("No")))
  
  return(p)
}