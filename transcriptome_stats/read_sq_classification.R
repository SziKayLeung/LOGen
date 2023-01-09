#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: prepare SQANTI classification files to be read in for downstream plotting etc
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## gene_summary_info
## corr_exon_length_num
##   
## ---------- Notes -----------------
## 
## Correlation reported is only of multi-exonic isoforms (mono-exonic isoforms are removed)
## Pre-requisite: 
## 1. SQANTI classification file has a tabulated "FL" column for each isoform with sum reads across each sample
## 2. function: draw_density.R from aesthetics_basics_plots.R


## ------------------- annotate_class_binary

# Aim: create extra columns in the classification files for binary purposes
  # i.e. whether within 50bp of cage peak, 50bp of TTS, 50bp of TSS etc
# Note:
  # within_peaks are all assumed within 50bp
  # RNA-Seq support assumed if min_cov > 1, then supported otherwise not supported 
  # within_polyA site assumed if there is a polyA motif
  # 3'fragment assumed if structural category is ISM with 3 prime fragment (assumed partially degraded products)
# Input:
  # class.files = df: SQANTI classification file after running SQANTI_class_preparation()
# Output:
  # class.files = df: additional columns

annotate_class_binary  <- function(class.files){
  
  class.files <- class.files %>% 
    mutate(within_50_cage = ifelse(abs(dist_to_cage_peak) <= 50 & !is.na(dist_to_cage_peak), "Within 50bp","Not within 50bp"),
           within_50_TTS = ifelse(abs(diff_to_TTS) <= 50 & !is.na(diff_to_TTS), "Within 50bp","Not within 50bp"),
           within_50_TSS = ifelse(abs(diff_to_TSS) <= 50 & !is.na(diff_to_TSS), "Within 50bp","Not within 50bp"),
           RNASeq_supported = ifelse(min_cov >= 1, "Supported","Not Supported"),
           within_polya_site = ifelse(is.na(polyA_motif),"No","Yes"),
           ISM_3prime = ifelse(structural_category == "ISM" & subcategory == "3prime_fragment" , "Yes","No"))
  
  return(class.files)
}


## ------------------- SQANTI_class_preparation

# Aim: read and wrangle classification file generated from SQANTI
# Note: function adapted from Liz Tseng (SQANTI_report2.R)
# Input:
  # path.class.file = str: path of classification file generated from SQANTI
  # standard = str: <standard/ns> 
    # "standard" assumes that the quantification of multiple samples are present in file (<sample>.FL)
    # therefore recreate FL column by summing across all .FL 
    # create ISOSEQ_TPM and Log_ISOSEQ_TPM columns
# Output: df 

SQANTI_class_preparation <- function(path.class.file,standard){
  
  cat("Loading classification file:",class.file,"\n")
  data.class = read.table(class.file, header=T, as.is=T, sep="\t")
  rownames(data.class) <- data.class$isoform
  
  xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron")
  xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  
  legendLabelF1 <- levels(as.factor(data.class$coding))
  
  data.class$structural_category = factor(data.class$structural_category,
                                          labels = xaxislabelsF1, 
                                          levels = xaxislevelsF1,
                                          ordered=TRUE)
  
  data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
  data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
  data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
  
  # Label Empty blanks in associated_gene column as "Novel Genes_PB_<isoform_ID>"
  data.class[data.class$associated_gene == "",]
  data.class$associated_gene[data.class$associated_gene == ""] <- paste0("novelGene_PB.", word(data.class$isoform[data.class$associated_gene == ""],c(2), sep = fixed ('.')))
  
  # Create a new attribute called "novelGene"
  data.class$novelGene <- "Annotated Genes"
  data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
  data.class$novelGene = factor(data.class$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  # Create a new attribute called "exonCat"
  data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
  data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
  data.class$exonCat = factor(data.class$exonCat,
                              levels = c("Multi-Exon","Mono-Exon"),
                              ordered=TRUE)
  
  data.class$all_canonical = factor(data.class$all_canonical,
                                    levels = c("canonical","non_canonical"),
                                    ordered=TRUE)
  
  data.class$within_cage_peak = factor(data.class$within_cage_peak)
  data.class$within_cage_peak <- factor(data.class$within_cage_peak, c("True","False"))
  
  if(standard == "standard"){
    # relabel the sum of the samples FL due to demultiplexing 
    # starts with FL refer to columns with samples
    # can append total directly to FL as do not change the order of the rows
    dat <- data.class %>% dplyr::select(starts_with("FL.")) %>% mutate(total =  rowSums(.[1:ncol(.)]))
    data.class$FL <- dat$total
    
    # convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
    total_fl <- sum(data.class$FL, na.rm=T)
    #print(paste0("Total FL counts:", total_fl))
    data.class$ISOSEQ_TPM <- data.class$FL*(10**6)/total_fl
    data.class$Log_ISOSEQ_TPM <- log10(data.class$ISOSEQ_TPM)
  }
  
  # further annotate classification files by within_cage etc.
  data.class <- annotate_class_binary(data.class)

  return(data.class)
}


## ------------------- SQANTI_gene_preparation

# Aim: generate the number of isoforms associated with gene (binned)
  # Note: function adapted from Liz Tseng (SQANTI_report2.R)
# Input:
  # data_class_output_file = df: SQANTI classification file after processing SQANTI_class_preparation()

SQANTI_gene_preparation <- function(data_class_output_file){
  
  # Make "isoPerGene" which is aggregated information by gene
  #  $associatedGene - either the ref gene name or novelGene_<index>
  #  $novelGene      - either "Novel Genes" or "Annotated Genes"
  #  $FSM_class      - "A", "B", or "C"
  #  $geneExp        - gene expression info
  #  $nIso           - number of isoforms associated with this gene
  #  $nIsoCat        - splicing complexity based on number of isoforms
  
  if (!all(is.na(data_class_output_file$gene_exp))){
    isoPerGene = aggregate(data_class_output_file$isoform,
                           by = list("associatedGene" = data_class_output_file$associated_gene,
                                     "novelGene" = data_class_output_file$novelGene,
                                     "FSM_class" = data_class_output_file$FSM_class,
                                     "geneExp"=data_class_output_file$gene_exp),
                           length)
  } else {
    isoPerGene = aggregate(data_class_output_file$isoform,
                           by = list("associatedGene" = data_class_output_file$associated_gene,
                                     "novelGene" = data_class_output_file$novelGene,
                                     "FSM_class" = data_class_output_file$FSM_class),
                           length)
  }
  # assign the last column with the colname "nIso" (number of isoforms)
  colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"
  
  
  isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class, 
                                 levels = c("A", "B", "C"), 
                                 labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"), 
                                 ordered=TRUE)
  
  isoPerGene$novelGene = factor(isoPerGene$novelGene, 
                                levels = c("Annotated Genes", "Novel Genes"), 
                                ordered=TRUE)
  
  # Altered to extend out the graphs
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5,7,9, max(isoPerGene$nIso)+1), labels = c("1", "2-3", "4-5", "6-7","8-9",">=10"))
  Sample_Type <- data_class_output_file$Sample[1]
  isoPerGene$Sample <- Sample_Type
  
  return(isoPerGene)
}


## ------------------- SQANTI_remove_3prime

# Aim: remove ISM 3'fragments that are assumed to be degraded
# Input:
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()

SQANTI_remove_3prime <- function(class.files){
  
  output <- class.files %>% 
    mutate(cate = paste0(structural_category,"_", subcategory)) %>%
    filter(cate!= "ISM_3prime_fragment")
  
  return(output)
}


## ------------------- targeted_remove_3ISM

# Aim: remove ISM 3'fragments that are assumed to be degraded and retain only associated isoforms of target genes
# Input:
  # TargetGenelist = vec: list of target genes of interest  
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()

targeted_remove_3ISM <- function(TargetGenelist, class.files){
  
  class.files <- SQANTI_remove_3prime(class.files) %>% filter(associated_gene %in% TargetGenelist) 
  class.files["FL"] <- rowSums(class.files[, grepl( "FL." , names(class.files))])
  
  return(class.files)
}


## ------------------- subset_class_phenotype

subset_class_phenotype <- function(class.files, phenotype_file, condition){
  
  class.files <- annotate_class_binary(class.files) 
  cols = c("isoform", "min_cov","associated_gene", "exons", "length", "dist_to_cage_peak", 
           "within_50_cage", "dist_to_polya_site","within_50_TTS","within_50_TSS",
           "within_polya_site","polyA_motif","structural_category","subcategory",
           "diff_to_TTS","diff_to_TSS","diff_to_gene_TTS","diff_to_gene_TSS","structural_category")
  
  class.files <- class.files %>% 
    dplyr::select(cols, paste0("FL.", phenotype_file[phenotype_file$Phenotype == condition,"Sample.ID"])) %>% 
    mutate(TotalFL = rowSums(.[paste0("FL.", phenotype_file[phenotype_file$Phenotype == condition,"Sample.ID"])])) %>%
    filter(TotalFL > 0) %>% mutate(Dataset = condition) 
  
  return(class.files)
}