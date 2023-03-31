#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: generate plots of the summary of all the target genes at a global level (IR and ES) after running FICLE
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## dendro_AS_dataset
## plot_dendro_Tgene
## plot_n_events
## plot_ES_Tgene
## plot_IR_Tgene
## plot_A5A3_Tgene
## plot_FirstExon_Tgene
##
## ---------- Notes -----------------
##  
## Pre-requisite: Run FICLE on each target gene to generate (A5A3_tab.csv/Exonskipping_tab.csv/Intron) files across all the genes
##


## ---------- Packages -----------------

suppressMessages(library("ggdendro"))
suppressMessages(library("cowplot"))


## ------------------- input_FICLE_splicing_results_tGene

# Aim: read in individual files (.csv) from target gene stats folder 
# Notes:
  # To be used for downstream functions (plot_ES_Tgene(), plot_IR_TGene(), plot_A5A3_Tgene(), plot_FirstExon_Tgene())
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene name of interest
  # filename = str: suffix of name that is common across the file of interest  
# Output:
  # read in file

input_FICLE_splicing_results_tGene <- function(TG_anno_dir, gene, filename){
  
  file = paste0(TG_anno_dir,"/", gene,"/Stats/", gene, filename)
  cat("Read in:", file, "\n")
  file_df = if(!file.size(file) == 0){read.csv(file)}else{as.data.frame()}
  
  return(file_df)
}



## ------------------- dendro_AS_dataset

# Aim: return a value for the plot_dendro_T_gene() from the input file <gene>_Exonskipping_generaltab.csv
# Note: <gene>_Exonskipping_generaltab.csv generated from FICLE to tabulate the status of each exon in dataset

dendro_AS_dataset <- function(value){
  
  if(value %in% c("FirstExon","Present_FinalExon")){return("Present")
    }else if (value %in% c("No")){return("Present")
      }else if (value %in% c("Yes")){return("ES")
        }else if (value %in% c("IR_FinalExon")){return("IR")
          }else if (value %in% c("NA_FinalExon")){return("Present")
          }else if (is.na(value)){return("Absent")
              }else{value}

}


## ------------------- plot_dendro_Tgene

# Aim: plot the dendogram (overview of the isoforms by exon status) of each gene
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene of interest
# Note:
  # clustering lines are not plotted for mouse Apoe due to too memory intensive with >2000 isoforms
# Output:
  # dendroplot of the gene of interest

plot_dendro_Tgene <- function(TG_anno_dir, gene){
  
  # read in files from FICLE
  # gene_tab = documents for each transcript and exon whether it's skipped, or contains IR
  gene_Exontab = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_Exon_tab.csv")
  tab = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_Exonskipping_generaltab.csv")
  
  # plot the cluster lines of the isoforms (grouped by splicing events) using a dendrogram
  # for mouse Apoe, there are over 2000 isoforms therefore avoid generating the cluster lines as memory intensive
  if(gene != "Apoe"){
    gene_Exontab[gene_Exontab=="1001"]<-3
    gene_Exontab = gene_Exontab %>% column_to_rownames(., var = "X")
    cluster = hclust(dist(gene_Exontab[,2:ncol(gene_Exontab)]))
    hcd = as.dendrogram(cluster)
    ddata <- dendro_data(hcd, type = "rectangle")
    clusterl = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +   theme_dendro() + labs(title = "\n")
  }
  
  # wide to long and classify the splicing event
  tab = tab %>% gather(Gene, Value, -X) 
  for(i in 1:nrow(tab)){tab[["Col"]][i] = dendro_AS_dataset(tab[["Value"]][i])}
  tab$X <- factor(tab$X , levels=cluster$labels[cluster$order])
  
  tab = tab %>% mutate(Col = factor(Col, levels = c("Present","ES","IR","Absent"))) %>%
    mutate(gencode = as.numeric(word(Gene, c(2), sep = "_")))
  
  # plot
  n1 <- length(unique(tab$X))
  n2 <- length(unique(tab$gencode))
  p <- ggplot(tab, aes(x = as.factor(gencode), y = X)) +
    geom_tile(aes(fill = Col)) + labs(x = "") + 
    scale_fill_manual(name = "Classification", 
                      labels = c("Present","ES","IR","Absent"),
                      values = c(alpha(wes_palette("Royal2")[5],0.7),
                                 alpha(wes_palette("Royal1")[2],0.7),
                                 wes_palette("Zissou1")[1],
                                 alpha(wes_palette("Chevalier1")[3])),drop = FALSE) +
    mytheme_font + labs(y = "Transcripts") + 
    theme(legend.position = "top", axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_x_discrete(breaks = seq(1,max(tab$gencode), by = 2)) +
    geom_line(data = data.frame(x = c(0, n2) + 0.5, y = rep(2:n1, each = 2) - 0.5),
              aes(x = x, y = y, group = y), linetype="dotted") 
  
  if(gene != "Apoe"){
    output_hmap = plot_grid(clusterl,p, rel_widths = c(0.3,0.7))
  }else{
    output_hmap = plot_grid(p)
  }
  
  return(output_hmap)
}


## ------------------- plot_n_events

# Aim: plot the number of events for gene (basic function)
# Input:
  # df_tally: df = read in files from FICLE 
  # <gene>_Exonskipping_tab.csv/ <gene>_IntronRetention_tab.csv
  # totaln = num: number of isoforms associated with gene (corresponding to df_tally)
  # as_event = str: <IR/ES>; used to determine the x-axis name and colour
# Output:
  # bar-plot of the number of isoforms with 
  # either the number of exons skipped (as_event == "ES") or the number of IR events (as_event == "IR)

plot_n_events <- function(df_tally, totaln, as_event){
  
  # tally the number of exon_skipping events per transcript
  n_events <- df_tally %>% group_by(transcript_id) %>% tally() %>% dplyr::rename(events = n)
  
  if(as_event == "ES"){
    x_lab = "Number of exons skipped"
    col = alpha(wes_palette("Royal1")[2],0.2)
  }else{
    x_lab = "Number of IR events"
    col = alpha(wes_palette("Zissou1")[[1]],0.2)
  }
  
  p <- n_events %>% group_by(events) %>% tally() %>% mutate(perc = n/totaln * 100) %>%
    mutate(num = "num") %>%
    ggplot(., aes(x = as.factor(events), y = n, fill = num)) + geom_bar(fill="white", stat = "identity", colour = col)  + 
    geom_text(aes(label = paste0(round(perc,1),"%")), position = position_stack(vjust = 0.5)) +
    labs(y = "Number of Isoforms", x = x_lab) + mytheme + theme(legend.position = "none") 
  
  return(p)
}


## ------------------- plot_ES_Tgene 

# Aim: plot the exon-skipping related plots per target gene (using results from FICLE)
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene of interest
  # merged.class.files = df: SQANTI classification file of merged dataset (ensure the isoform ID is the same as the IDs in TG_anno_dir)
# Output:
  # p1: bar-plot of the number of isoforms with specific exon across the gene skipped
  # p2: bar-plot of the number of isoforms with the total number of exons skipped (per isoform)

plot_ES_Tgene <- function(TG_anno_dir, gene, merged.class.files){
  
  cat("********* Processing: ", as.character(gene), "\n")
  
  # number of associated isoforms with gene from merged dataset
  niso = nrow(merged.class.files[merged.class.files$associated_gene == gene,])
  
  # read in file from FICLE
  ES_tally = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_Exonskipping_tab.csv")
  
  cat("Number of transcripts with exon skipping:", length(unique(ES_tally$transcript_id)),"\n")
  
  # tally the number of exon_skipping events per transcript
  nES_events <- ES_tally %>% group_by(transcript_id) %>% tally() %>% dplyr::rename(ES = n)
  
  # plots
  p1 <- ES_tally %>% group_by(ES) %>% tally() %>% 
    mutate(col = as.numeric(word(ES, c(2), sep = fixed("_")))) %>% 
    mutate(col_order = factor(col, levels = sort(as.numeric((col))))) %>%
    ggplot(., aes(x = col_order, y = n)) + geom_bar(stat = "identity") + mytheme + 
    labs(x = "Gencode Exon Skipped", y = "Number of isoforms")
  
  p2 <- plot_n_events(ES_tally,niso,"ES")
  
  return(list(p1,p2))
}


## ------------------- plot_IR_Tgene 

# Aim: plot intron-retention related plots per target gene (using results from FICLE)
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene of interest
  # merged.class.files = df: SQANTI classification file of merged dataset (ensure the isoform ID is the same as the IDs in TG_anno_dir)
# Output:
  # p1: bar-plot of the number of isoforms with the number of intron retention events
  # p2: bar-plot of the number of isoforms with specific exons wtih IR

plot_IR_Tgene <- function(TG_anno_dir, gene, merged.class.files){
  
  cat("********* Processing: ", as.character(gene), "\n")
  
  # number of associated isoforms with gene from merged dataset
  niso = nrow(merged.class.files[merged.class.files$associated_gene == gene,])
  
  # read in file from FICLE
  IR_tally = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_IntronRetention_tab.csv")
  
  totalt = nrow(IR_tally)
  cat("Number of IR events: ", totalt, "\n")
  cat("Number of transcripts with IR: ", length(unique(IR_tally$transcript_id)),"\n")
  
  if(length(unique(IR_tally$transcript_id)) != 0){
    # tally the number of transcripts with intron retention across each exon
    IR_exon = IR_tally %>% group_by(IR) %>% tally() %>% 
      mutate(IR_type = word(IR, c(1), sep = fixed("_")),
             gencode = word(IR, c(3), sep = fixed("_"))) %>% 
      group_by(IR_type, gencode) %>% 
      tally(n)
    #print(IR_exon)
    
    # tally the number of exon_skipping events per transcript
    nIR_events <- IR_tally %>% group_by(transcript_id) %>% tally() %>% dplyr::rename(IR = n)
    
    # plots
    p1 <- plot_n_events(IR_tally,niso,"IR")
    
    p2 <- IR_exon %>% 
      mutate(perc = n/totalt * 100) %>% 
      mutate(xlabel = factor(gencode, levels = sort(as.numeric(unique(gencode))))) %>%
      ggplot(., aes(x = xlabel, y = n, fill = IR_type, group = IR_type)) + geom_bar(stat = "identity") + 
      geom_text(aes(label = paste0(round(perc,0),"%")), position = position_stack(vjust = 0.5)) +
      mytheme + labs(y = "Number of Isoforms", x = "Exon with IR") + 
      scale_fill_manual(name = "IR Classification", label = c("IR", "IR Match"), 
                        values = c(alpha(wes_palette("Zissou1")[[1]],0.2),wes_palette("Zissou1")[[1]])) + 
      theme(legend.position = "top")
  }else{
    p1 <- NULL
    p2 <- NULL
  }
  
  
  return(list(p1,p2))
}


## ------------------- plot_A5A3_Tgene

# Aim: plot A5A3' related plots per target gene (using results from FICLE)
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene of interest
# Output:
  # p1: bar-plot the number of isoforms with the A5' and A3' splice sites by classification

plot_A5A3_Tgene <- function(TG_anno_dir, gene){
  
  # read in file from FICLE
  A5A3_tally = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_A5A3_tab.csv")
  
  A5A3_num = A5A3_tally %>% group_by(cate, gencode_exon) %>% tally() %>% 
    mutate(xlabel = as.numeric(word(gencode_exon, c(2), sep = fixed("_"))))
  #print(A5A3_num)
  
  if(gene == "Vgf"){
    print(length(unique(A5A3_tally[A5A3_tally$gencode_exon == "Gencode_6","transcript_id"])))
  }
  
  p <- A5A3_num %>%
    mutate(cate = factor(cate, levels = c("ExtendedA5","ExtendedA3","TruncatedA5","TruncatedA3"))) %>%
    filter(cate != "NA") %>%
    ggplot(.,aes(x = as.factor(xlabel), y = n, fill = cate)) + geom_bar(stat = "identity") + 
    labs(x = "Exon with alternative sites", y = "Number of Isoforms") + mytheme + 
    scale_fill_manual(name = "Classification", labels = c("A5' Extension", "A3' Extension", "A5' Truncation","A3' Truncation"), 
                      values = c(wes_palette("IsleofDogs1")[1], alpha(wes_palette("IsleofDogs1")[1],0.4),
                                 wes_palette("IsleofDogs2")[3], alpha(wes_palette("IsleofDogs2")[3],0.4))) #+ 
    #theme(legend.position = "top") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    #scale_x_discrete(breaks = plot_break_every_nth(n = 2))
  
  if(gene == "Sorl1"){
    p = p + scale_x_discrete(breaks = plot_break_every_nth(n = 3))
  }
  
  return(p)
}


## ------------------- plot_FirstExon_Tgene

# Aim: identify which exon is commonly used as the first exon of the isoforms in dataset
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # gene = str: gene of interest  
# Output:
  # p: bar-plot of the number of isoforms with the different said exons used as first exon

plot_FirstExon_Tgene <- function(TG_anno_dir, gene){
  
  # read in file from FICLE
  FirstExon_tally = input_FICLE_splicing_results_tGene(TG_anno_dir, gene, "_Exonskipping_generaltab.csv")
  
  dat = melt(FirstExon_tally, id= "X") %>% filter(value == "FirstExon") %>% group_by(variable) %>% tally() %>% 
    mutate(xlabel = as.numeric(word(variable, c(2), sep = fixed("_"))))
  #print(dat)
  
  p <- ggplot(dat, aes(x = as.factor(xlabel), y = n)) + geom_bar(stat = "identity",fill= alpha(wes_palette("Royal2")[5],0.7)) + 
    mytheme + labs(x = "Exon as first exon", y = "Number of Isoforms") + 
    scale_x_discrete(breaks = plot_break_every_nth(n = 2))
  
  ## Picalm ## 
  # Shorter isoforms with alternative first exon, check the number that is detected by both ONT and Iso-Seq
  #table(targeted.class.files[targeted.class.files$isoform %in% gene_Exontab[is.na(gene_Exontab$Gencode_1),"X"],"Dataset"])
  
  return(p)
}