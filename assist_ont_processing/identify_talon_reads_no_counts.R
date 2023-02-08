#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## output a list of transcripts that are documented in the TALON gtf but not in the abundance file
## sanity check for TALON pipeline
## Input:
##  --gtf = TALON gtf generated from talon_create_GTF
##  --counts = TALON abundance file
##  --output = output path and name of file
## --------------------------------


## ---------- packages -----------------

suppressMessages(library("data.table"))
suppressMessages(library("optparse"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="TALON gtf", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="TALON abundance file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file and path of missing transcripts in abundance file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ---------- Input files -----------------

cat("Processing",opt$gtf,"\n")
gtf <- fread(opt$gtf)
colnames(gtf) <- c("chr","source","type","start","end","score","strand","phase","attributes")

cat("Processing",opt$counts,"\n")
counts <- read.table(opt$counts, header = T)   


## ---------- function -----------------

# the problem is the attributes column that tends to be a collection
# of the bits of information you're actually interested in
# in order to pull out just the information I want based on the 
# tag name, e.g. "gene_id", I have the following function:
# vased on https://www.biostars.org/p/272889/
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}


## ---------- run function -----------------

# extract transcript_id from gtf
transcripts <- gtf[type == "transcript"]
transcripts$transcript_id <- unlist(lapply(transcripts$attributes, extract_attributes, "transcript_id"))

# find differences between abundance and gtf
missing_abundance <- setdiff(transcripts$transcript_id, counts$annot_transcript_id)
paste("Missing transcripts in abundance, but present in gtf:", length(missing_abundance))

# write output
write.table(missing_abundance,opt$output, quote = FALSE, row.names = FALSE, col.names = FALSE)