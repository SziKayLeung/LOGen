#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Generate summary tables and plots of alternative splicing events for dataset
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## summarise_AS_events
## summarise_AS_events_bygenes_isoforms
## summarise_AS_events_bygenes_isoforms
##   
## ---------- Notes -----------------
## 
## Pre-requisite: 
## 1. Generate output tables from running Alternative-splicing bash scripts


## ------------------- summarise_AS_events

# Aim: extract alternative splicing event related stats into simplified table (overview without isoform level) 
# Input:
  # events_genes_output = path of file tabulating the number of events across each gene (for each dataset)  
    # 4 columns: <associated_gene> <Event> <n> <Sample>
# Output: 
  # list 1: all_stats = total number of splicing events per dataset 
  # list 2: dataset_stats = the number and proportion of each splicing event by dataset 
  # list 3: genes_events_stats = the number of genes with splicing events (cumulative) by dataset 

summarise_AS_events <- function(events_genes_output){
  
  # read in the number of events per gene 
  splicing_events <- read.csv(events_genes_output) %>% mutate(Event = as.character(Event))
  splicing_events$Event[splicing_events$Event == "SE"] <- "ES"
  
  # tally the number of events by sample/dataset
  dataset_tally_events <- splicing_events %>% group_by(Sample) %>% tally(n)  
  
  # tally the number of genes with splicing events by sample/dataset
  dataset_tally_gene <- splicing_events %>% group_by(Sample) %>% dplyr::count(associated_gene) %>% tally()
  
  # Number of splicing events - all (inc Novel Genes)
  splicing_events_tally <- list()
  for(i in 1:length(unique(splicing_events$Sample))){
    splicing_events_tally[[i]] <- splicing_events %>% group_by(Event, Sample) %>% tally(n) %>% 
      left_join(dataset_tally_events, by = "Sample") %>% 
      mutate(perc = n.x/n.y * 100) %>% filter(Sample == unique(splicing_events$Sample)[i]) %>% 
      mutate(Type = "All") %>% 
      mutate(Group = "1") %>% select(!Sample)
  }
  names(splicing_events_tally) <- unique(splicing_events$Sample)
  
  # Number of splicing events across genes (percentage) by sample/dataset
  splicing_events_number <- splicing_events %>% group_by(Sample, associated_gene) %>% tally() %>% 
    `colnames<-`(c("Sample", "associated_gene", "n_events")) %>% 
    group_by(n_events, Sample) %>% 
    tally() %>% left_join(dataset_tally_gene, by = "Sample") %>%
    mutate(perc = n.x/n.y * 100)  %>% 
    `colnames<-`(c("Number_of_Splicing_Events", "Sample", "Genes", "Total_Genes","perc")) %>% 
    as.data.frame(.) 
  
  # output
  output <- list(dataset_tally_events, splicing_events_tally, splicing_events_number)
  names(output) <- c("all_stats","dataset_stats", "genes_events_stats")
  
  return(output)
}


## ------------------- summarise_AS_events_bygenes_isoforms

# Aim: extract alternative splicing event related stats into simplified table at the gene and isoform level
# Input:
  # events_known_genes = path of file tabulating the number and proportion of genes across each event
    # 3 columns: <Event> <n> <perc>
  # events_isoforms = path of file tabulating the number of known and novel isoforms across each event
    # 3 columns: <Event> <known> <novel>
# Output: 
  # df of <Event> <Type:AnnotatedGenes/knownIsoform/novelIsoform> <perc> <Group> 
    # note: percentage relates to the individual group (i.e. at the gene or novel/known isoform level)
  
summarise_AS_events_bygenes_isoforms <- function(events_known_genes, events_isoforms){
  
  # read in input 1: number of splicing events for known genes
  knowngenes_AS <- read.csv(events_known_genes) %>% 
    # create column Group and Type for downstream merging and plotting
    mutate(Group = "1", Type = "Annotated Genes", Event = as.character(Event))
  knowngenes_AS$Event[knowngenes_AS$Event == "SE"] <- "ES"
  
  # read in input 2: number of known and novel isoforms by events
  # gather input to --> <Event=A3/A5..> <Type=known/novel> <n>
  anno_novel_AS <- read.csv(events_isoforms) %>% gather(., Type, n, known:novel, factor_key=TRUE) %>% mutate(Event = as.character(Event))
  anno_novel_AS$Event[anno_novel_AS$Event == "SE"] <- "ES"
  
  # tabulate the total number of known and novel isoforms with AS events
  anno_novel_AS_tally <- anno_novel_AS %>% group_by(Type) %>% tally(n)
  
  # comprehensiely tabulate with percentages of the number of known and novel isoforms with AS events
  anno_novel_AS_tally_present <- anno_novel_AS %>% left_join(anno_novel_AS_tally, by = "Type") %>% mutate(perc = n.x/n.y * 100, Group = "2")
  
  # output
  output <- bind_rows(knowngenes_AS[,c("Event","Type","perc","Group")], anno_novel_AS_tally_present[,c("Event","Type","perc","Group")])
  return(output)
}


## ------------------- AS_genes_events

# Aim: wrapper to plot the number of AS events at the gene and isoform level
# Input:
  # events_genes_output = path of file tabulating the number of events across each gene (for each dataset)    
  # events_known_genes = path of file tabulating the number and proportion of genes across each event
  # events_isoforms = path of file tabulating the number of known and novel isoforms across each event
  # dataset = name of "dataset" to plot specific p2 
# Output: 
  # p1: box-plot of the number of splicing events by annotated genes, known and novel isoforms
  # p2: bar-plot of the number of genes by the number of different types of splicing events

AS_genes_events <- function(events_genes_output, events_isoforms, events_known_genes, dataset){
   
  summarised_bygene_isoform <- summarise_AS_events_bygenes_isoforms(events_known_genes, events_isoforms)
  summarised_stats <- summarise_AS_events(events_genes_output)
  
  # plots
  p1 <-  ggplot(summarised_bygene_isoform, aes(x = Type, y = perc, fill = reorder(Event, -perc))) + geom_bar(stat = "identity") + 
    theme_bw() + labs(y = "Splicing Events (%)", x = "") + mytheme + 
    theme(legend.position = "bottom",legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) +
    facet_grid(~ Group,scales='free') +
    scale_x_discrete(labels = c("known" = "Known \n Isoforms", "novel" =  "Novel \n Isoforms")) +
    theme(strip.background = element_blank(),strip.text.x = element_blank())

  p2 <- summarised_stats$genes_events_stats %>%
   filter(Sample %in% dataset) %>%
    ggplot(., aes(x = Number_of_Splicing_Events, y = perc)) + 
    geom_bar(stat = "identity", position = position_dodge()) + 
    theme_bw() + labs(y = "AS Genes (%)", x = "Number of Splicing Events") + mytheme + 
    theme(legend.position = c(0.85,0.85), legend.title = element_blank()) + 
    scale_x_continuous(breaks = 1:7) + geom_text(aes(label = round(perc,1)),nudge_y = 2)
    
  return(list(p1,p2))
}
