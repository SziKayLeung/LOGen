#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## merge the classification file with abundance file (for counts in the same file) and to only include target genes
## --filter ==> filter classification file to keep only isoforms supported by a minimum number of reads and samples
## Functions: 
##  subset_class_by_targets() 
##  quantify_class_abundance()
##  filter_class_by_counts
## Input:
##  --classfile = SQANTI classification file
##  --expression = expression matrix for merging (isoform column)
##  --target_genes = tsv file of list of genes for subsetting (one column with no header)
##  --filter = TRUE/FALSE (turn on filtering)
##  --nsample, --nreads = filtering threshold
## --------------------------------


## ---------- packages -----------------

#suppressMessages(library("data.table"))
suppressMessages(library("optparse"))

# source general scripts to input SQANTI class.files
LOGEN = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
sink("/dev/null")
suppressMessages(sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv))
sink()

## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-s", "--classfile"), type="character", default=NA, help="SQANTI classification file", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NA, help="expression matrix"),
  make_option(c("-g", "--target_genes"), type="character", default=NA, help="tsv file of list of target genes, no header", metavar="character"),
  make_option(c("-f", "--filter"), action="store_true", default=FALSE, help="Apply filtering thresholds"),
  make_option(c("--nsample"), type="integer", default=5L, help="filtering: minimum number of samples [default %default]", metavar="number"),
  make_option(c("--nreads"), type="integer", default=10L, help="filtering: minimum number of reads [default %default]", metavar="number"),
  make_option(c("-d", "--dir"), type="character", default=NA, help="output directory, optional", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.na(opt$dir)){
  opt$dir = dirname(opt$classfile)
}
message("Directory of output files: ", opt$dir)


## ---------- Input files -----------------

if(!is.na(opt$classfile)){
  message("Reading in:",opt$classfile)
  input.class.files <- SQANTI_class_preparation(opt$classfile,"nstandard")
}else {
  stop("classfile must be provided. See script usage (--help)")
}

if(!is.na(opt$expression)){
  message("Reading in:",opt$expression)
  input.abundance <- read.csv(opt$expression, row.names = 1)
}else {
  stop("expression abundance must be provided. See script usage (--help)")
}

if (!is.na(opt$target_genes)){
  cat("Reading",opt$target_genes,"\n")
  targetgenes = array(read.table(opt$target_genes)[["V1"]])
  message("Filtering by target genes: ", paste0(targetgenes,sep=","))
}

## ---------- function -----------------

# Total number of isoforms and genes prior to filtering for target genes
cat("Total number of isoforms:", nrow(input.class.files),"\n")
cat("Total number of genes:", length(unique(input.class.files$associated_gene)),"\n")

# keep only the target genes
if(!is.na(opt$target_genes)){
  targeted.class.files <- subset_class_by_targets(input.class.files, targetgenes)
  cat("Total number of isoforms:", nrow(targeted.class.files),"\n")
  cat("Total number of genes:", length(unique(targeted.class.files$associated_gene)),"\n")
}

# merge the abundance with the classification file
if(!is.na(opt$target_genes)){
  quantified.class.files <- quantify_class_abundance(targeted.class.files, input.abundance)
}else{
  quantified.class.files <- quantify_class_abundance(input.class.files, input.abundance)
}

# apply filter by count
if(opt$filter){
  filtered.quantified.targeted.class.files <- filter_class_by_counts(quantified.class.files, nread_threshold=opt$nreads, nsample_threshold=opt$nsample)
  if(!is.na(opt$target_genes)){
    filename <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".targetgenes_counts_filtered.txt")
    filenameiso <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".targetgenes_filtered_isoforms.txt")
  }else{
    filename <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".counts_filtered.txt")
    filenameiso <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".filtered_isoforms.txt")
  }
  param_filename <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),"/count_filtered_params.txt")
  #message("Filtering parameters:", file = param_filename)
  #cat("Minimum number of samples:", opt$nsample, "\n", file = param_filename, append = TRUE)
  #cat("Minimum number of reads:", opt$nreads, "\n", file = param_filename, append = TRUE)
  #cat("Target Genes:", paste(targetgenes, collapse = ","), "\n", file = param_filename, append = TRUE)
  #cat("Writing to:", filename,"\n")
  write.table(filtered.quantified.targeted.class.files, filename, quote = F, sep = "\t")
  
  # filename of isoforms kept only
  write.table(filtered.quantified.targeted.class.files$isoform, filenameiso, quote = F, sep = "\t", row.names = F, col.names = F)
}else{
  if(!is.na(opt$target_genes)){
    filename <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".targetgenes_counts.txt")
  }else{
    filename <- paste0(opt$dir,"/",str_remove(basename(opt$classfile),".txt"),".counts.txt")
  }
  message("Writing to:", filename)
  write.table(quantified.class.files, filename, quote = F, sep = "\t")
}
