#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: identify and generate a txt file tabulating isoforms from sqanti classification file 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 


## ---------- packages -----------------

suppressMessages(library("optparse"))
suppressMessages(library("stringr"))

LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "transcriptome_stats/read_sq_classification.R"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-s", "--classfile"), type="character", default=NULL, help="SQANTI classification file", metavar="character"),
  make_option(c("-d", "--output_dir"), type="character", default=NULL, help="Output directory; default is directory of SQANTI classification file", metavar="character"),
  make_option (c("-g","--genes"),default=NULL,help="comma separated list of genes to filter")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$dir)){
  opt$dir = paste0(dirname(opt$classfile),"/")
}
opt$name <- tools::file_path_sans_ext(basename(opt$classfile))

## ---------- Input files -----------------

message(paste0("Processing: ", opt$classfile))
input.class.files <- SQANTI_class_preparation(opt$classfile,"nstandard")


## ---------- Output files -----------------

message("Output:",opt$dir,opt$name,"isoform.txt")
write.table(input.class.files$isoform, paste0(opt$dir,opt$name,"_isoform.txt"),quote=F,row.names=F,col.names=F)



