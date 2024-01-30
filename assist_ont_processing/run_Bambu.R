#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## create a general R script to run Bambu
## --------------------------------


## ---------- packages -----------------

suppressMessages(library("optparse"))
suppressMessages(library("bambu"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-b", "--bam"), type="character", 
              help="aligned bam file", metavar="character"),
  make_option(c("-i", "--input_dir"),  default=NULL, 
              help="directory containing input aligned bam file"),
  make_option(c("-f", "--fa"), type="character",  
              help="genome fasta file", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,   
              help="genome gtf file", metavar="character"),
  make_option(c("-a", "--annotations_RDS"), metavar="character", default=NULL, 
              help="annotations RDS file"),
  make_option(c("-o", "--output_dir"), type="character", 
              help="output directory", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$gtf) & is.null(opt$annotations_RDS)){
  print("Error: either --gtf or --annotations_RDS required")
}

## ---------- input -----------------

if(is.null(opt$input_dir)){
  test.bam <- opt$bam
}else{
  test.bam <- list.files(path = opt$input_dir, pattern = "bam", full = T)
}
message("Processing following bam file:")
for (bam_file in test.bam) {
  message(bam_file)
}

if(!is.null(opt$annotations_RDS)){
  annotations <- readRDS(opt$annotations_RDS)
}else{
  annotations <- prepareAnnotations(opt$gtf)
}

message("Running Bambu")
se <- bambu(reads = test.bam, annotations = annotations, genome = opt$fa)
save(se, file = paste0(opt$output_dir, "/se.RData"))

message("Write Bambu output")
writeBambuOutput(se, path = opt$output_dir)
