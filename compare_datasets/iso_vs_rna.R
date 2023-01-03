#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Compare Iso-Seq vs RNA-Seq defined transcriptome 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## compare_datasets_withincage
## compare_rna_vs_iso_transcriptome
## correlate_rnaseq_isoseq_counts
## 
## ---------- Notes -----------------
## 
## Pre-requisite: 
## 1. Generate RNA-Seq defined transcriptome using stringtie
## 2. Compared RNA-Seq vs Iso-Seq transcriptome using gffcompare 
## 3. function: SQANTI_gene_preparation from sqanti_general.R


## ------------------- compare_datasets_withincage

# Aim: generate a table of the proportion of isoforms within cage peak (50bp) between two datasets
# Input:
  # df1 & 2 = classification file of dataset 1 and 2
  # name1 & 2 = name of dataset 1 and 2
# Output:
  # df of 3 columns: dataset(name 1 & 2); category (within/without); iso_proportion (%) 

compare_datasets_withincage <- function(df1, df2, name1, name2){
  
  within_cage <- list(
    df1 = nrow(df1 %>% filter(abs(dist_to_cage_peak) <= 50 )) / nrow(df1) * 100,
    df2 = nrow(df2 %>% filter(abs(dist_to_cage_peak) <= 50 )) / nrow(df2) * 100
  )
  
  df <- data.frame(dataset = c(name1, name2), 
                   within = c(within_cage$df1, within_cage$df2),
                   without = c(100 - within_cage$df1, 100 - within_cage$df2)) %>% 
    reshape2::melt(id = "dataset") %>%
    mutate(variable = factor(variable, levels = c("without","within"))) %>%
    `colnames<-`(c("dataset", "category", "iso_proportion"))
  
  return(df)
}


## ------------------- compare_rna_vs_iso_transcriptome

# Aim: Generate plots comparing iso-seq vs rna-seq transcriptome
# Input:  
  # class.files = df: classification file generated from SQANTI after alignment to standard reference genome
  # rnaseq.class.files = df: classification file generated from SQANTI after generating stringie transcriptome from rna-seq reads
  # files from Gffcompare: cuffrefmap 
  # files from Gffcompare: cufftmap_input
# Function: 
  # SQANTI_gene_preparation from sqanti_general.R
# Note:
  # Only include isoforms from Iso-Seq that are partially or fully matched to RNA-Seq! Not all isoforms from Iso-Seq dataset!
# Output:
  # p1: venn diagram of the number of isoforms in RNA-Seq vs Iso-Seq transcriptome
  # p2: bar-plot of number of isoform categories per transcriptome
  # p3: box-plot of the isoform length per transcriptome
  # p4: box-plot of the number of exons per transcriptome 
  # p5: bar-plot of the number of isoforms within 50bp CAGE peak
  # p6: bar-plot of the number of isoforms per structural category 

compare_rna_vs_iso_transcriptome <- function(class.files,rnaseq.class.files,cuffrefmap,cufftmap){

  # Replace id of those rnaseq isoforms that fully match with isoseq isoforms with pbid for venn diagram
  rnaseq <- c(as.character(cufftmap[cufftmap$class_code != "=","qry_id"]), as.character(cuffrefmap[cuffrefmap$class_code == "=","ref_id"])) 
  
  # Number of isoforms per dataset 
  subset_cols <- c("isoform", "associated_gene", "novelGene","FSM_class","gene_exp","Sample")
  num_iso <- list(class.files %>% mutate(Sample = "Iso-Seq")  %>%  select(all_of(subset_cols)),
                  rnaseq.class.files %>% mutate(Sample = "RNA-Seq") %>%  select(all_of(subset_cols)))
  isoPerGene <- lapply(num_iso, function(x) SQANTI_gene_preparation(x)) %>% bind_rows()
  
  # Total Number of Genes per Type 
  Total_Num <- isoPerGene %>% group_by(Sample) %>% dplyr::count(Sample)
  
  # Plots
  p1 <- venn.diagram(
    x = list(rnaseq, class.files$isoform), category.names = c("RNA-Seq","Iso-Seq"), filename = NULL, output=TRUE,
    lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC"), main = "\n",
    cex = 1,fontface = "bold",fontfamily = "ArialMT",
    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  
    cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n",
    print.mode = "raw"
  )
  
  p2 <- isoPerGene %>% group_by(Sample) %>% dplyr::count(nIsoCat) %>% 
    full_join(Total_Num,., by = "Sample") %>% mutate(Perc = n.y/n.x * 100) %>% 
    ggplot(., aes(x=nIsoCat, fill=Sample)) +
    geom_bar(stat="identity", aes(y= Perc, group = as.factor(Sample)), 
             color="black", linewidth=0.3, width=0.7, position="dodge") + 
    labs(x ="Number of Isoforms", y = "Genes (%)", fill = "") +
    mytheme + 
    theme(legend.position = c(0.75,0.95))
  
  # Length
  p3 <- bind_rows(class.files,rnaseq.class.files) %>% 
    ggplot(., aes(x = Sample, y = log10(length/1000))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Isoform Length (Log10 kb)")
  
  # Exons
  p4 <- bind_rows(class.files,rnaseq.class.files) %>% 
    ggplot(., aes(x = Sample, y = log10(exons))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Number of Exons (Log10)")
  
  # cage peak
  p5 <- compare_datasets_withincage(class.files, rnaseq.class.files,"Iso-Seq","RNA-Seq") %>%
    ggplot(., aes(x = dataset, y = iso_proportion, fill = category)) + geom_bar(stat = "identity") + mytheme +
    labs(y ="Isoforms within 50bp CAGE (%)", x = "", fill = "") + theme(legend.position = "right") +
    scale_fill_discrete(labels = c("No","Yes"))
  
  # structural_category
  p6 <- bind_rows(tabulate_structural_cate(class.files) %>% mutate(Sample= "Iso-Seq"),
                  tabulate_structural_cate(rnaseq.class.files) %>% mutate(Sample= "RNA-Seq")) %>% 
    ggplot(., aes(x = Sample, y = perc, fill = structural_category)) + geom_bar(stat = "identity") + mytheme +
    labs(x = "", y = "Isoforms (%)", fill = "Structural Category") + theme(legend.position = "left")
  
  return(list(p1,p2,p3,p4,p5,p6))
}


## ------------------- correlate_rnaseq_isoseq_counts

# Aim: correlate the number of reads from RNASeq alignment to Iso-Seq and RNA-Seq defined transcriptome 
# Input: 
  # class.files = df: classification file generated from SQANTI after alignment to standard reference genome
  # files from Gffcompare: cuffrefmap 
  # kallisto_rnaseq = rna-seq defined counts from alignment of short-reads to long-reads using kallisto
  # kallisto_isoseq = iso-seq defined counts from alignment of long-reads to final curated Iso-Seq transcriptome using kallisto
# Note: using only the 23,761 isoforms considered matching from gffcompare (cuffrefmap file from gffcompare output)
# Output: 3 density plots correlating the counts of RNA-Seq aligned to Iso-Seq vs Iso-Seq alone

correlate_rnaseq_isoseq_counts <- function(class.files, cuffrefmap, kallisto_rnaseq, kallisto_isoseq){

  # only consider the isoforms that are "matching" between Iso-Seq and RNA-Seq defined transcriptome (output from gff compare)
  matching <- cuffrefmap[cuffrefmap$class_code == "=",] %>% mutate(qry_id = word(qry_id_list,c(2),sep = fixed("|")))
  
  # number of isoforms matching = 23,761
  # nrow(cuffrefmap[cuffrefmap$class_code == "=",]) 
  
  # ref_id = pacbio id from isoseq-defined transcriptome gtf
  matching_pbisoform <- cuffrefmap[cuffrefmap$class_code == "=","ref_id"]
  
  # qry_id_list = id from rnaseq-defined transcriptome gtf
  matching_rnaseqisoform <- unique(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|")))
  
  # QC no repeated refid (PBid from isoforms)
  #n_occur <- data.frame(table(cuffrefmap[cuffrefmap$class_code == "=","ref_id"]))
  #n_occur[n_occur$Freq > 1,]
  
  # repeated qry_id due to imperfect match from cuffrefmap i.e. two PB isoforms refer to the same RNA-Seq isoform if 5' and 3' end different but internal junction the same
  #n_occur <- data.frame(table(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|"))))
  #n_occur[n_occur$Freq > 1,]
  
  ### Kallisto
  IsoSeq_FL <- class.files[class.files$isoform %in% matching_pbisoform,c("isoform","FL","ISOSEQ_TPM")]
  kallisto_isoseqmatching <- kallisto_isoseq %>% filter(target_id %in% matching_pbisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("PBID", "RNA2Isoseq_Counts","RNA2Isoseq_TPM"))
  RNASeq_Defmatching <- kallisto_rnaseq %>% filter(target_id %in% matching_rnaseqisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("RNASeqID", "RNA2RNAseq_Counts","RNA2RNAseq_TPM"))
  
  
  # Merge all the counts with the isoforms that are considered matching 
  final <- merge(matching,kallisto_isoseqmatching ,by.x = "ref_id", by.y = "PBID", all = T)
  final <- merge(final,RNASeq_Defmatching ,by.x = "qry_id", by.y = "RNASeqID", all = T)
  final <- merge(final,IsoSeq_FL ,by.x = "ref_id", by.y = "isoform", all = T)
  
  # Correlation tests 
  cor.test(final$RNA2Isoseq_TPM,final$RNA2RNAseq_TPM, method = "pearson")
  cor.test(final$RNA2Isoseq_TPM,final$ISOSEQ_TPM, method = "pearson")
  cor.test(final$RNA2RNAseq_TPM,final$ISOSEQ_TPM,method = "pearson")
  
  
  p1 <- density_plot(final,"RNA2Isoseq_TPM","ISOSEQ_TPM", "RNA-Seq Transcript Expression from \n Iso-Seq-defined transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  p2 <- density_plot(final,"RNA2Isoseq_TPM","RNA2RNAseq_TPM", "RNA-Seq Expression from \n Iso-Seq-defined transcriptome (TPM)", "RNA-Seq expression from RNA-Seq transcriptome (TPM)","")
  p3 <- density_plot(final,"RNA2RNAseq_TPM","ISOSEQ_TPM", "RNA-Seq expression from RNA-Seq transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  
  return(list(p1,p2,p3))
}