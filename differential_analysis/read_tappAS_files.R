## ---------- Script -----------------
##
## Purpose: read in tappAS-related files
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)



## ------------------- read_dea_files 

read_dea_files <- function(dir){
  f = excel_sheets(dir) %>% purrr::set_names() %>% purrr::map(read_excel, path = dir)
  return(f)
}


## ------------------- input_tappasfiles

input_tappasfiles <- function(tappas_input_dir){
  # read in files generated from TAPPAS
  tappasfiles <- list(paste0(tappas_input_dir,"/Data/gene_matrix.tsv"),
                      paste0(tappas_input_dir,"/Data/gene_transcripts.tsv"),
                      paste0(tappas_input_dir,"/InputData/input_normalized_matrix.tsv"),
                      paste0(tappas_input_dir,"/Data/result_gene_trans.tsv"),
                      paste0(tappas_input_dir,"/Data/transcript_matrix.tsv"),
                      paste0(tappas_input_dir,"/Data/DIU/result_trans.tsv"))
  
  tappasfiles <- lapply(tappasfiles, function(x){print(x); read.table(x, sep = "\t", header = T)})
  names(tappasfiles) <- c("gene_matrix","gene_transcripts","input_normalized_matrix","result_gene_trans",
                          "transcript_matrix","results_DIU")
  
  if (file.exists(paste0(tappas_input_dir,"/Data/DEA/result_trans_CASEvsCONTROL.tsv"))){
    tappasfiles$results_trans <- read.table(paste0(tappas_input_dir,"/Data/DEA/result_trans_CASEvsCONTROL.tsv"),sep = "\t", header = T)
  }else{
    tappasfiles$results_trans <- read.table(paste0(tappas_input_dir,"/Data/DEA/result_trans.tsv"),sep = "\t", header = T)
  }
  
  if (file.exists(paste0(tappas_input_dir,"/Data/DEA/result_trans_CASEvsCONTROL.tsv"))){
    tappasfiles$results_gene <- read.table(paste0(tappas_input_dir,"/Data/DEA/result_gene_CASEvsCONTROL.tsv"),sep = "\t", header = T)
  }else{
    tappasfiles$results_gene <- read.table(paste0(tappas_input_dir,"/Data/DEA/result_gene.tsv"),sep = "\t", header = T)
  }
  
  return(tappasfiles)
}


## ------------------- annotate_tappasfiles 

# annotate_tappasfiles
# prerequisite for plots
# aim1: annotate the tappas Normalised counts file with the isoforms from the class file and join by phenotype
# aim2: deduce the gene expression by the sum of the expression of the filtered isoforms 
# output: norm_transcounts (normalised transcript counts) & GeneExp
annotate_tappasfiles <- function(classification_file, tappas_normalised_expmatrix, phenotype){
  
  #classification_file = class.files
  #tappas_normalised_expmatrix = loaded$iso$input_normalized_matrix
  #phenotype = tappasiso_phenotype
  # Annotate the Normalised counts matrix from tappas with the associated gene and sample phenotypes for plots 
  Norm_transcounts = 
    # annotate the tappas output Normalised counts matrix of transcripts by the associated gene name and transcript name
    merge(classification_file [,c("isoform","associated_gene","associated_transcript","structural_category")], 
          tappas_normalised_expmatrix, by.x = "isoform", by.y = 0) %>% 
    reshape2::melt(id = c("isoform", "associated_gene", "associated_transcript","structural_category")) %>% 
    # rename 
    dplyr::rename(sample = variable) %>% 
    # annotate the samples to the phenotype 
    left_join(., phenotype, by = "sample") %>% 
    # change factor levels for plots
    mutate(structural_category=recode(structural_category, `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",`novel_in_catalog`="NIC"),
           Isoform = paste0(isoform,"_",structural_category))
  
  # Deduce gene expression from the sum of normalised transcript counts 
  GeneExp = Norm_transcounts %>% 
    group_by(associated_gene,sample) %>% 
    dplyr::summarise(Exp = sum(value)) %>%
    left_join(., phenotype, by = "sample")
  
  output <- list(Norm_transcounts,GeneExp)
  names(output) <- c("Norm_transcounts","GeneExp")
  return(output)
}


# filter differential transcript results (NoiSEQ)
diff_results <- function(class.files, loaded_file, level, prob_threshold){
  
  if(level == "transcript"){
    dat <- loaded_file[["results_trans"]] 
  }else{
    dat <- loaded_file [["results_gene"]]
  }
  
  dat <- dat %>% 
    mutate(`1-prob` = 1 - prob) %>% 
    filter(`1-prob` < prob_threshold) %>% 
    dplyr::rename(., control_mean = X1_mean, case_mean = X2_mean)
  
  if(level == "transcript"){
    dat <- merge(class.files[,c("isoform","associated_gene","associated_transcript","structural_category")],
                 dat, by.x = "isoform", by.y = "transcript") %>% 
      select(isoform, associated_gene, prob, `1-prob`, log2FC, control_mean, case_mean)
  }else{
    dat <- dat %>% select(gene, prob, `1-prob`, log2FC, control_mean, case_mean) %>%
      dplyr::rename(associated_gene = gene)
  }
  
  dat <- dat %>% arrange(-desc(`1-prob`))
  
  return(dat)
}


# annotate differential transcript expression results from tappAS
diff_trans_stats <- function(result_trans, class.files){
  
  dat <- merge(class.files[,c("isoform","associated_gene","associated_transcript","structural_category")],
               result_trans, by.x = "isoform", by.y = "transcript")
  
  dat <- dat %>% arrange(-desc(`prob`))
  
  return(dat)
}

