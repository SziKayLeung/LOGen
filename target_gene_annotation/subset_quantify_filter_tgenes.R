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
suppressMessages(sapply(list.files(path = paste0(LOGEN,"transcriptome_stats"), pattern="*.R", full = T), source,.GlobalEnv))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-s", "--classfile"), type="character", default=NULL, help="SQANTI classification file", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NULL, help="expression matrix"),
  make_option(c("-g", "--target_genes"), type="character", default=NULL, help="tsv file of list of target genes, no header", metavar="character"),
  make_option(c("-f", "--filter"), action="store_true", default=FALSE, help="Apply filtering thresholds"),
  make_option(c("--nsample"), type="integer", default=5L, help="filtering: minimum number of samples [default %default]", metavar="number"),
  make_option(c("--nreads"), type="integer", default=10L, help="filtering: minimum number of reads [default %default]", metavar="number")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ---------- Input files -----------------

cat("Processing",opt$classfile,"\n")
input.class.files <- SQANTI_class_preparation(opt$classfile,"nstandard")

cat("Processing",opt$expression,"\n")
input.abundance <- read.csv(opt$expression, row.names = "id")

cat("Reading",opt$target_genes,"\n")
targetgenes = array(read.table(opt$target_genes)[["V1"]])


## ---------- function -----------------

# Total number of isoforms and genes prior to filtering for target genes
cat("Total number of isoforms:", nrow(input.class.files),"\n")
cat("Total number of genes:", length(unique(input.class.files$associated_gene)),"\n")

# keep only the target genes
cat("Filtering:",targetgenes,"\n")
targeted.class.files <- subset_class_by_targets(input.class.files, targetgenes)
cat("Total number of isoforms:", nrow(targeted.class.files),"\n")
cat("Total number of genes:", length(unique(targeted.class.files$associated_gene)),"\n")

# merge the abundance with the classification file
quantified.targeted.class.files <- quantify_class_abundance(targeted.class.files, input.abundance)

# apply filter by count
if(opt$filter){
  filtered.quantified.targeted.class.files <- filter_class_by_counts(quantified.targeted.class.files, nread_threshold=opt$nreads, nsample_threshold=opt$nsample)
  filename <- paste0(str_remove(opt$classfile,".txt"),".targetgenes_counts_filtered.txt")
  param_filename <- paste0(dirname(opt$classfile),"/count_filtered_params.txt")
  cat("Filtering parameters:\n", file = param_filename)
  cat("Minimum number of samples:", opt$nsample, "\n", file = param_filename, append = TRUE)
  cat("Minimum number of reads:", opt$nreads, "\n", file = param_filename, append = TRUE)
  cat("Writing to:", filename,"\n")
  write.table(filtered.quantified.targeted.class.files, filename, quote = F, sep = "\t")
  
  # filename of isoforms kept only
  filenameiso <- paste0(str_remove(opt$classfile,".txt"),".targetgenes_filtered_isoforms.txt")
  write.table(filtered.quantified.targeted.class.files$isoform, filenameiso, quote = F, sep = "\t", row.names = F, col.names = F)
}else{
  filename <- paste0(str_remove(opt$classfile,".txt"),".targetgenes_counts.txt")
  cat("Writing to:", filename,"\n")
  write.table(quantified.targeted.class.files, filename, quote = F, sep = "\t")
}
