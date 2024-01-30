#!/usr/bin/env Rscript
## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: Create tappAS expression matrix from talon abundance file
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
## 
## Date: 26/01/2022
## ---------- Notes -----------------
##
## 
##   
##
## 

## ---------- Packages and arguments -----------------

suppressMessages(library("optparse"))


# handle command line arguments - don't want to use external package so users don't have to install it
option_list = list(
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="Input TALON expression file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="out", 
              help="output name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


## ---------------------------

# read in abundance files
abundance.files <- read.table(opt$expression, sep = "\t", as.is = T, header = T)

# set row names for expression matrix (ensure unique names)
rownames(abundance.files) <- abundance.files$annot_transcript_id

# subset dataframe to only include abundance columns (column 12 onwards)
exp_matrix <- abundance.files[,12:length(abundance.files)]
cat("\nNumber of samples:", ncol(exp_matrix))

# write output 
write.table(exp_matrix, paste0(opt$output,"/", opt$name,".txt"), quote=F, sep = "\t")
