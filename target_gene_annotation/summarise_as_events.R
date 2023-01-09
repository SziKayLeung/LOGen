#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: generate plots of the summary of all the target genes at a global level (IR and ES) after running FICLE
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## input_FICLE_splicing_results
## plot_summarised_AS_events
## bin_num_events
## plot_summarised_ES
## plot_summarised_IR
##
## ---------- Notes -----------------
##  
## Pre-requisite: Run FICLE on each target gene to generate (A5A3_tab.csv/Exonskipping_tab.csv/Intron) files across all the genes
##


## ------------------- input_FICLE_splicing_results

# Aim: read in the list of global files with specified filenames and merge into one dataframe
# Input:
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # filenames = st: <A5A3_tab.csv/Exonskipping_tab.csv/Intron>
# Output:
  # df: all file list merged into one dataframe

input_FICLE_splicing_results <- function(TG_anno_dir, filenames){
  
  # list the files with specific filename pattern in directory and read if file size != 0
  files <- list.files(path = TG_anno_dir, pattern = filenames, recursive = TRUE, full = T)
  print(files)
  tab = lapply(files, function(x)  if(!file.size(x) == 0){read.csv(x)})
  
  # save the files under the names of the gene, which is teh first saved character before "_"
  names(tab) = lapply(list.files(path = TG_anno_dir, pattern = filenames, recursive = TRUE), function(x) word(x, c(1), sep = fixed("/")))
  
  # combind all files
  df = bind_rows(tab, .id = "associated_gene") 
  
  return(df)
}


## ------------------- plot_summarised_AS_events

# Aim: plot the number of AS events across target genes after running FICLE
# Pre-requisite:
  # Run FICLE to generate <TargetGene>_Final_Transcript_Classifications.csv and <TargetGene>_A5A3_tab.csv
# Input:
  # Merged_gene_class_df = df: read and merged <TargetGene>_Final_Transcript_Classifications.csv to record AS events across all target genes
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # ref_gencode = df: output of FICLE reference extract_reference_info.py
    # associated_gene; Number_of_Gencode_Isoforms; Transcript_Length; Number_of_Exons; MaxGeneLength; MaxTransLength; Maxexons; MeanRNASeqCounts
# Output:
  # p1: bar-plot of the percentage of splicing events across target genes 
  # p2: bar-plot of the number of isoforms with different A5 A3 classifications 
  # p3: bar-plot of th percentage of isoforms with A5A3 splice sites in first, internal and last exons

plot_summarised_AS_events <- function(Merged_gene_class_df, TG_anno_dir, ref_gencode){
  
  # only interested in these AS events to plot
  AS = c("AP","ES","IR","A5A3","AT")
  AS_df = subset(Merged_gene_class_df, row.names(Merged_gene_class_df) %in% AS) %>% rownames_to_column(., "AS") %>% reshape2::melt(id="AS") %>% dplyr::rename("associated_gene" = "variable") 
  total_AS_df = AS_df %>% group_by(associated_gene) %>% tally(value)
  
  # p1: bar-plot of the percentage of splicing events across target genes 
  p1 = merge(AS_df,total_AS_df) %>% mutate(perc = value/n * 100) %>%
    mutate(AS = factor(AS, levels = c("AP","ES","IR","A5A3","AT"))) %>%
    ggplot(.,aes(x = associated_gene, y = perc, fill = AS)) + geom_bar(stat = "identity") + 
    scale_fill_manual(values = c(wes_palette("Moonrise1")[3],wes_palette("Royal1")[2],wes_palette("Zissou1")[[1]],
                                 wes_palette("IsleofDogs1")[1],wes_palette("Moonrise1")[4])) + mytheme + 
    labs(x = "", y = "Percentage of Splicing Events (%)") + theme(legend.position = "top") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  total_AS = AS_df %>% group_by(AS) %>% tally(value)
  print(total_AS)
  sum(total_AS$n)
  
  # read in A5A3 files from FICLE 
  A5A3_all <- input_FICLE_splicing_results(TG_anno_dir, "A5A3_tab.csv") %>% dplyr::rename("Isoform" = "transcript_id")
  
  # p2: bar-plot of the number of isoforms with different A5 A3 classifications 
  p2 = A5A3_all %>% group_by(associated_gene,cate) %>% tally() %>% 
    mutate(cate = factor(cate, levels = c("ExtendedA5","ExtendedA3","TruncatedA5","TruncatedA3",
                                          "TruncatedBothA3A5","ExtendedBothA3A5","ExtendedTruncatedA3")))  %>%
    filter(cate != "ExtendedTruncatedA3") %>%
    ggplot(., aes(x = reorder(associated_gene,-n), y = n, fill = cate)) + geom_bar(stat = "identity") + 
    labs(x = "", y = "Number of Isoforms (Thousands)") + mytheme +
    theme(legend.position = c(0.7,0.7)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(name = "Classification", labels = c("A5' Extension", "A3' Extension", "A5' Truncation","A3' Truncation",
                                                          "A5' A3' Extension","A5' A3' Truncation"), 
                      values = c(wes_palette("IsleofDogs1")[1], alpha(wes_palette("IsleofDogs1")[1],0.4),
                                 wes_palette("IsleofDogs2")[3], alpha(wes_palette("IsleofDogs2")[3],0.4),
                                 wes_palette("GrandBudapest2")[1], wes_palette("Moonrise3")[1])) +
    scale_y_continuous(labels = ks)
  
  # input reference number of exons for each target gene
  # ref_maxexon = maximum number of exons
  A5A3_all = A5A3_all %>% mutate(exon = as.factor(paste0(associated_gene,"_", word(gencode_exon,c(2), sep = fixed("_")))))
  ref_maxexon <- ref_gencode %>% mutate(code = paste0(associated_gene,"_", Maxexons), class = "Last") %>% select(associated_gene,code,class)
  ref_firstexon <- ref_gencode %>% mutate(code = paste0(associated_gene,"_1"), class = "First") %>% select(associated_gene, code, class)
  refexons <- rbind(ref_maxexon, ref_firstexon)
  
  dat = merge(A5A3_all, refexons, by.x = "exon", by.y = "code", all.x = T) %>% mutate(class = as.character(class))
  dat$class[is.na(dat$class)] <- "internal"
  
  nA5A3 = dat %>% group_by(associated_gene.x) %>% tally()
  
  # p3: bar-plot of th percentage of isoforms with A5A3 splice sites in first, internal and last exons
  p3 = dat %>% group_by(associated_gene.x, class) %>% tally() %>% 
    left_join(., nA5A3, by = c("associated_gene.x" = "associated_gene.x")) %>% 
    mutate(perc = n.x/n.y * 100) %>%
    ggplot(., aes(x = associated_gene.x, y = perc, fill = class)) + geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    mytheme + labs(x = "", y = "Percentage of isoforms with A5A3 splice sites (%)") + 
    theme(legend.position = "top") +
    scale_fill_manual(name = "Exon", labels = c("First", "Internal", "Last"), 
                      values = c(wes_palette("IsleofDogs1")[2], alpha(wes_palette("IsleofDogs1")[1],0.4),
                                 wes_palette("IsleofDogs2")[4])) 
  
  return(list(p1,p2,p3))
}


## ------------------- bin_num_events

# Aim: Place the number of events into specific bins (used for plot_summarised_ES/IR)
# Input:
  # num = numeric: number of events
# Output:
  # num: returned bin classification

bin_num_events <- function(num){
  
  if(num <= 5){return(num)
  }else if(6 <= num & num <= 10){return("6-10")
  }else if(10 < num & num <= 15){return("10-15")
  }else{return(">15")
  }
  
}


## ------------------- plot_summarised_ES

# Aim: plot the number of exon skipping events across target genes, further divided by constitutive/alternative status
# Input:
  # Gene_class = df: output from annotation of target gene using FICLE <associated_gene> <Isoform> <Matching> <.....> 
    # each column records the number of transcripts with said alternative splicing event
  # class.files = df: SQANTI classification file 
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # ref_gencode = df: output of FICLE reference extract_reference_info.py
    # associated_gene; Number_of_Gencode_Isoforms; Transcript_Length; Number_of_Exons; MaxGeneLength; MaxTransLength; Maxexons; MeanRNASeqCounts
  # ref_altcon = df: output of FICLE reference extract_reference_info.py of exon status as constitutive or alternative
    # exon; exon_status; associated_gene
# Output:
  # p1: bar-plot of the number of unique exons skipped (divided by constitutive/alternative) across target genes
  # p2: bar-plot of the total number of exons skipped (divided by constitutive/alternative) across target genes
  # p3: ratio of the number of exons skipped/number of transcripts per target gene
  # p4: total number of isoforms with exon skipping events (coloured by bins)

plot_summarised_ES <- function(gene_class, class.files, TG_anno_dir, ref_gencode, ref_altcon){
  
  # read in ES files from FICLE
  ES_tab_df <- input_FICLE_splicing_results(TG_anno_dir, "Exonskipping_tab.csv")
 
  # unique number of exon skipped, is it constitutive or alternative 
  ES_gene_count = ES_tab_df %>% group_by(associated_gene,ES) %>% tally() 
  
  # percentage of each target gene with constitutive and alternative exons
  exon_type <- merge(ref_altcon, ES_gene_count, by.x=c("associated_gene", "exon"), by.y=c("associated_gene", "ES")) %>% 
    group_by(associated_gene, exon_status) %>% tally() %>% dplyr::rename("total_ES" = "n") %>% 
    full_join(., ref_altcon %>% group_by(associated_gene) %>% tally() %>% dplyr::rename("total_exon" = "n"), by = "associated_gene") %>% 
    mutate(perc = paste0(round(total_ES/total_exon * 100,0),"%"))
  
  # total number of exons that are constitutive/alternative
  n_exon_type <- merge(ref_altcon, ES_gene_count, by.x=c("associated_gene", "exon"), by.y=c("associated_gene", "ES")) %>% 
    group_by(associated_gene, exon_status) %>% tally(n) 
  
  # total number of transcripts that have constitutive/alternative exons
  exon_type_transcripts <- merge(ref_altcon, ES_gene_count, by.x=c("associated_gene", "exon"), by.y=c("associated_gene", "ES")) %>% 
    group_by(associated_gene, exon_status) %>% tally(n) %>% dplyr::rename("total_ES" = "n") %>%
    full_join(., class.files %>% group_by(associated_gene) %>% tally() %>% dplyr::rename("total_trans" = "n")) %>%
    mutate(ratio = total_ES/total_trans) 
  
  # the number of exon-skipping events per target gene
  # using output from FICLE
  gene_class_df = bind_rows(gene_class, .id = "associated_gene") %>% dplyr::rename("Isoform" = "X") %>% subset(., select = -c(isoform))
  nES = gene_class_df %>% group_by(associated_gene, ES) %>% tally()
  
  # bin the number of AS events 
  for(i in 1:nrow(nES)){nES[["col_group"]][i] = bin_num_events(nES[["ES"]][i])}
  
  # tally the number of transcripts with binned splicing events
  nES = nES %>% mutate(col_group = factor(col_group, levels = c("0","1","2","3","4","5","6-10","10-15",">15"))) %>%
    group_by(associated_gene, col_group) %>% tally(n)
  
  # plots
  p1 <- ggplot(exon_type, aes(x = associated_gene, y = total_ES, fill = exon_status)) + geom_bar(stat = "identity") + 
    geom_text(aes(label = perc), angle = 90, position = position_stack(vjust = .5)) + 
    labs(x = "", y = "Number of unique exons skipped") + mytheme + 
    scale_fill_manual(name = "Classification", labels = c("Alternative","Constitutive"),
                      values = c(alpha(label_colour("isoseq"),0.5),label_colour("rnaseq"))) +
    theme(legend.position = c(0.3,0.85)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  p2 <- ggplot(n_exon_type, aes(x = associated_gene, y = n, fill = exon_status)) + geom_bar(stat = "identity", position = position_dodge()) +
    mytheme + 
    scale_fill_manual(name = "Exon Status", values = c(alpha(label_colour("isoseq"),0.5),label_colour("rnaseq"))) +
    theme(legend.position = "bottom") + labs(x = "", y= "Number of exons skipped") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p3 <- ggplot(exon_type_transcripts, aes(x = associated_gene, y = ratio, fill = exon_status)) + geom_bar(stat = "identity") +
    mytheme + 
    scale_fill_manual(name = "Exon Status", values = c(alpha(label_colour("isoseq"),0.5),label_colour("rnaseq"))) +
    theme(legend.position = "top") + labs(x = "", y= "number of exons skipped/number of trancripts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p4 <- ggplot(nES, aes(x = as.factor(associated_gene), y = n, fill = forcats::fct_rev(col_group))) + geom_bar(stat = "identity") + 
    scale_fill_manual(name = "ES", values = c("#F21A00","#E86F00","#E2B306","#E8C31E","#CAC656","#88BAA3","#5DAABC","#3B9AB2",
                                              alpha(wes_palette("Royal1")[1],0.5))) + 
    mytheme + labs(x = "", y = "Number of Isoforms (Thousands)") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "right") + 
    scale_y_continuous(breaks=number_ticks(5), labels = ks)

  return(list(p1,p2,p3,p4))  
}


## ------------------- plot_summarised_IR

# Aim: plot the number of intron retention events across target genes
# Input:
  # merged.class.files = df: SQANTI classification file of merged dataset (ensure the isoform ID is the same as the IDs in TG_anno_dir)
  # TG_anno_dir = str: path of directory containing output of FICLE 
  # ont_abundance = df: TALON abundance file (ensure the annot_transcript_id is the same as the ONT_isoform in merged.class.files)
# Note:
  # merged.class.files first processed by identify_dataset_specific_name() in dataset_identifier.R
# Output:
  # p1: total number of isoforms with intron retention events per target gene (coloured by bins)
  # p2: box-plot of the transcript expression of transcripts with different number of IR events
  # p3: same as p2 but classified by each target gene

plot_summarised_IR <- function(merged.class.files, TG_anno_dir, ont_abundance){

  # read in IR files from FICLE
  IR_tab_df <- input_FICLE_splicing_results(TG_anno_dir, "Intron")

  # tally the number of intron retention events per transcripts
  nIR_events <- IR_tab_df %>% group_by(associated_gene, transcript_id) %>% tally() %>% dplyr::rename(IR = n)
  nIR <- nIR_events %>% group_by(associated_gene,IR) %>% tally()

  for(i in 1:nrow(nIR)){nIR[["col_group"]][i] = bin_num_events(nIR[["IR"]][i])}

  # keep ONT isoforms only for expression given greater depth
  # need the merged dataframe for the isoform ID, given will be merging later with FICLE output (same ID)
  ont.class.files <- merged.class.files %>% mutate(IR_status = ifelse(.$isoform %in% IR_yes$transcript_id,"Yes","No")) %>% filter(Dataset != "Iso-Seq")
  
  # merge ont.class.files with expression
  ont.class.files.FL <- merge(ont.class.files, ont.abundance, by.x = "ONT_isoform", by.y = "annot_transcript_id", all.x = T)

  # merge the expression with the number of isoforms with IR events
  IR_exp = merge(ont.class.files.FL, nIR_events, by.x = "isoform", by.y = "transcript_id", all.x = T) %>% mutate(IR = ifelse(IR_status == "No", 0, IR))
  
  
  # plots
  p1 <- nIR %>% mutate(col_group = factor(col_group, levels = c("0","1","2","3","4","5","6-10","10-15",">15"))) %>% 
    filter(col_group != "0") %>%
    ggplot(., aes(x = associated_gene, y = n, fill = forcats::fct_rev(col_group))) + geom_bar(stat = "identity") + 
    scale_fill_manual(name = "IR", values = c("#CAC656","#88BAA3","#3B9AB2", alpha(wes_palette("Royal1")[1],0.5))) + 
    mytheme + labs(x = "", y = "Number of Isoforms") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = c(0.8,0.8)) + 
    scale_y_continuous(breaks=number_ticks(10))
  
  p2 <- ggplot(IR_exp, aes(x = IR, y = log10(FL.y), group = IR)) + geom_boxplot() + 
    mytheme + labs(x = "Number of IR events", y = "Transcript expression (log10)")
  
  p3 <- ggplot(IR_exp, aes(x = as.factor(IR), y = log10(FL.y), fill = as.factor(IR))) + geom_boxplot() + facet_grid(~associated_gene.x) +
    labs(x = "Transcripts with Intron Retention", y = "ONT FL read count (Log10)") + mytheme_font +
    scale_fill_discrete(name = "IR events") + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "bottom",
          strip.text.x = element_text(angle = 90))
  
  return(list(p1,p2,p3))
  
}


