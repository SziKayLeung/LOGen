#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: prepare SQANTI classification files to be read in for downstream plotting etc
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## annotate_class_binary
## SQANTI_class_preparation
## SQANTI_gene_preparation
## SQANTI_remove_3prime
## targeted_remove_3ISM
## subset_class_phenotype
## subset_class_by_sample
## subset_class_by_targets
## quantify_class_abundance
## filter_class_by_counts 
##   
## ---------- Notes -----------------
## 
## Some functions borrowed from SQANTI (SQANTI_report.R)
## ----------------------------------


## ---------- Packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))


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

SQANTI_class_preparation <- function(class.file,standard){
  
  cat("Loading classification file:",class.file,"\n")
  data.class = read.table(class.file, header=T, as.is=T, sep="\t")
  rownames(data.class) <- data.class$isoform
  
  # SQANTI versions (v5.0) change column names for dist_to_cage_peak and within_cage_peak 
  # change columns back for consistency
  if("within_CAGE_peak" %in% colnames(data.class)){
    cat("Processing SQANTI classification file generated from new version (5.0) onwards\n")
    data.class <- data.class %>% dplyr::rename(within_cage_peak = within_CAGE_peak,
                                        dist_to_cage_peak = dist_to_CAGE_peak,
                                        dist_to_polya_site = dist_to_polyA_site)
    
    if(standard != "all"){
      cat("Removing artifacts\n")
      data.class <- data.class %>% filter(filter_result == "Isoform")
    }
    
    data.class$structural_category[data.class$structural_category == "Genic_Genomic"] <- "Genic\nGenomic"
  }
  
  if("FSM" %in% unique(data.class$structural_category)){
    if("Genic_Genomic" %in% unique(data.class$structural_category)){
      xaxislevelsF1 <- c("FSM","ISM","NIC","NNC", "Genic_Genomic","Antisense","Fusion","Intergenic")
      xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic")
    }else{
      xaxislevelsF1 <- c("FSM","ISM","NIC","NNC", "genic","antisense","fusion","intergenic","genic_intron")
    }
  }else{
    xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron")
    xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  }
  
  
  if("dataset" %in% colnames(data.class)){
    data.class <- data.class %>% dplyr::rename(Dataset = dataset)
  }
  
  
  legendLabelF1 <- levels(as.factor(data.class$coding))
  
  data.class$structural_category = factor(data.class$structural_category,
                                          labels = xaxislabelsF1, 
                                          levels = xaxislevelsF1,
                                          ordered=TRUE)
  
  data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
  data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
  data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
  
  # Label Empty blanks in associated_gene column as "Novel Genes_PB_<isoform_ID>"
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
                                     "FSM_class" = data_class_output_file$FSM_class,
                                     "structural_category" = data_class_output_file$structural_category),
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

# Aim: subset classification file based on phenotype 
# Pre-requisites:
  # class.files to have FL count columns for samples 
  # phenotype file to have the sample names (<sample>) under "Sample.ID" column
  # phenotype file to have the condition to subset (under "Phenotype" column)
# Assumptions:
  # 0 FL reads across all samples in condition => not detected in condition
# Input:
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()
  # phenotype_file = df: read in phenotype file; <Sample.ID; condition>
  # condition = str of condition to subset (i.e. AD)
# Output:
  # class.files of the isoforms detected in the samples of interest 

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


## ------------------- subset_class_by_sample

# Aim: subset classification file by sample of interest; only include isoforms detected in sample 
# filter isoforms that have >0 reads for that sample
# Pre-requisites:
  # class.files to have FL count columns for samples
# Input:
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()
  # sample = str: sample name that matches the column name of class.files
# Output:
  # class.files of the isoforms detected in the sample of interest 

subset_class_by_sample <- function(class.files, sample){
  
  # grep the column name that matches with the sample
  col <- colnames(class.files)[[grep(sample, colnames(class.files))]]
  
  # filter the sample abundance > 0 reads for that specific transcript
  # filter_at(vars(1)) = column 1 = sample abundance column
  class.files <- class.files %>% select(all_of(col), associated_gene, isoform) %>% filter_at(vars(1), any_vars(. > 0))
  
  return(class.files)
}


## ------------------- subset_class_by_targets

# Aim: subset classification file by only including isoforms associated with target genes
# filter isoforms associated with exact match to target genes
# Input: 
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()
  # TargetGenes = list: target genes for subsetting (capitals do not make any difference on subset)
# Output:
  # class.files of the isoforms associated with target genes

subset_class_by_targets <- function(class.files, TargetGenes){
  
  ## grep target genes in fusion isoforms
  # subset fusion isoforms
  fusion_isoforms <- class.files %>% filter(structural_category == "fusion")
  
  # generate a list of target fusion isoforms to append
  lst = list()
  
  # loop through each target gene and append the isoform id
  lst <- append(sapply(TargetGenes, function(x) fusion_isoforms[grepl(x, fusion_isoforms$associated_gene),"isoform"]),lst)
  
  # subset the fusion target isoforms of interest
  target_fusion_isoforms <- fusion_isoforms[fusion_isoforms$isoform %in% array(unlist(lst)),]
  
  # subset all other non-fusion isoforms, where associated_gene is exact match
  # capitalise to ensure no inconsistency
  targeted.class.files <- class.files %>% filter(toupper(associated_gene) %in% toupper(TargetGenes))
  
  # rbind both classes of isoforms (non-fusion and fusion)
  targeted.class.files <- rbind(targeted.class.files, target_fusion_isoforms)
  
  return(targeted.class.files)
}


## ------------------- quantify_class_abundance

# Aim: merge the classification file with abundance (if --fl_count not turned on in SQANTI)
# obtain sample level counts for each isoform
# Pre-requisite:
  # isoform column in class.files = row.names in abundance.files
# Input: 
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()
  # abundance.files 
  # sum = <"yes","no"> : if sum then add additional columns of the number of reads and samples detected for each isoform
# Output:
  # class.files with addtional columns of the counts 
# Notes:
  # all the isoforms in the class.files should be retained 
  # NA in count columns suggest no abundance information for isoform (mismatched abundance and sqanti file?)
  # whereas unique isoforms in the abundance files are discarded (from upstream filtering)

quantify_class_abundance <- function(class.files, abundance.files){
  
  # generate two columns recording the 
  # number of samples with reads, and the 
  # total number of reads across all samples
  if("dataset" %in% colnames(abundance.files)){
    abundance.files.counts <- abundance.files %>% 
      # remove columns with sum and the dataset to avoid double counting
      select(-contains("sum_FL"),-dataset) 
  }else{
    abundance.files.counts <- abundance.files
  }
  
  #head(abundance.files.counts)
  abundance.files.counts <- abundance.files.counts %>% filter(rownames(.) != "0") %>% 
    mutate(nsamples = rowSums(.!=0), nreads = rowSums(.)) %>% 
    select(nsamples, nreads)

  # cbind additional two columns
  abundance.files <- abundance.files %>% tibble::rownames_to_column("isoform") %>% filter(isoform != 0)
  abundance.files <- cbind(abundance.files, abundance.files.counts)
  
  # merge class.files and abundance.files
  class.files <- merge(class.files, abundance.files, by = "isoform", all.x = TRUE)
  
  return(class.files)
}


## ------------------- filter_class_by_counts

# Aim: filter SQANTI classification by read counts
# Pre-requisite:
  # SQANTI classification file have cols: <nreads>, <nsamples>, <dataset>
# Input: 
  # class.files = df: SQANTI classification file after processing SQANTI_class_preparation()
  # nread_threshold = number: the minimum number of reads 
  # nsample_threshold = number: the minimum number of samples
# Output:
  # filtered class.files 
# Notes:
  # Isoforms that are detected in "Both" or "All" datasets are retained and not filtered regardless of read counts and samples

filter_class_by_counts <- function(class.files, nread_threshold, nsample_threshold){
  
  cat("Filtering isoforms with less than", nread_threshold, "reads, across", nsample_threshold,"samples\n")
  
  if("dataset" %in% colnames(class.files)){
    cat("Keeping isoforms in both or all datasets\n")
    both.class.files <- class.files %>% filter(dataset %in% c("Both","All"))
    unique.class.files <- class.files %>% filter(!dataset %in% c("Both","All")) %>% filter(nreads >= nread_threshold & nsamples >= nsample_threshold)
    filtered.class.files <- rbind(both.class.files, unique.class.files)
  }else{
    filtered.class.files <- class.files %>% filter(nreads >= nread_threshold & nsamples >= nsample_threshold)
  }
  
  cat("Removed", nrow(class.files) - nrow(filtered.class.files),"isoforms\n")
  
  return(filtered.class.files)
}

plot_sqanti_filtered_reasons <- function(reason.files){
  
  p <- reason.files %>% group_by(structural_category, reasons) %>% tally() %>%
    ggplot(., aes(x = reasons, y = n, fill = structural_category)) + geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Reasons Filtered", y = "Number of Isoforms") 
  
  return(p)
}
