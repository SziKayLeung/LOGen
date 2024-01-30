#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## 29/06/2021: Generate RNA-Seq expression file for TappAS after alignment of individual RNA-Seq samples with Kallisto 
## Prequisite: Generate Kallisto files in one directory with alignment to Iso-Seq scaffold 
## Input:
##  --kallisto = directory containing all the subdirectories of results with "abundance.tsv" from running Kallisto 
##  --counts = TALON abundance file
##  --whole = working with whole transcriptome datasets (no subsetting)
##  --targeted = further subsetting of counts abundance file to only include expression of isoforms associated with Target Genes
##  --output = output prefix name
## --------------------------------


## ---------- packages -----------------

suppressMessages(library("stringr"))
suppressMessages(library("tidyr"))
suppressMessages(library("optparse"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-k", "--kallisto"), type="character", default=NULL, 
              help="Path of input kallisto directory containing processed files", metavar="character"),
  make_option(c("--targeted"), action = "store_true", default = FALSE, 
              help="working with targeted transcriptome profiling dataset"),
  make_option(c("--whole"), action = "store_true", default = TRUE, 
              help="working with whole transcriptome dataset"),
  make_option(c("-s", "--sqfile"), type = "character", default=NULL,
              help="SQANTI classification file"),
  make_option(c("-g", "--genes"), type = "character", default=NULL,
              help="List of target genes for subsetting"),
  make_option(c("-o", "--output"), type="character", default="combined", 
              help="utput prefix name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## ---------- Input files -----------------

if(opt$whole & opt$targeted){
  cat("Error! Only one argument --whole or --targeted allowed\n")
  quit(status=1)
}else if(opt$whole & !opt$targeted){
  cat("Processing files from whole transcriptome profiling\n") 
}else{
  cat("Processing files from targeted transcriptome profiling\n") 
}

if(opt$targeted){
  if(opt$sqfile != NULL){
    cat("Processing sqanti classification file:", opt$sqfile)
    class.files = read.table(opt$sqfile, as.is = T, sep = "\t", header = T)
  }else{
    cat("Error! SQANTI classification file required")
    quit(status=1)
  }
  
  if(opt$genes != NULL){
    genes = strsplit(opt$genes, ",")
    cat("Subsetting by target genes", genes, "\n")
  }else{
    cat("Not subsetting by target genes")
  }
}


# list all the files in the rnaseq dir with the pattern "abundance.tsv"
# samples = the folder/sample from which the abundance.tsv came from
cat("Processing files in:", opt$kallisto,"\n")
abundance = lapply(list.files(path = opt$kallisto, recursive = T, include.dirs = T, pattern = "abundance.tsv", full.names = T),
                   function(x) read.table(x,header = T))
samples = word(list.files(path = opt$kallisto, recursive = T, include.dirs = T, pattern = "abundance.tsv"),c(1),sep=fixed("/"))


## ---------- process and combine ----------------

# only extract the relevant columns from each file 
mod_abundance = lapply(abundance, function(x) x[,c("target_id","est_counts")])
names(mod_abundance) = samples

# bind all the files together, with the sample as differentiator 
all = dplyr::bind_rows(mod_abundance, .id = "sample")

# spread from long to wide in the format required for tappAS 
counts = all %>% spread(sample, est_counts)
colnames(counts)[1] <- ""


## ---------- output ---------------- 

output_file = paste0(opt$kallisto,"/", opt$output, ".expression.txt")
cat("Writing output file:", output_file,"\n")

if(opt$whole){
  write.table(counts, output_file, quote = F, row.names = F,sep = "\t")
  
}else{
  
  # subset on the isoforms associated with target genes and the counts 
  TargetIsoforms <- class.files[toupper(class.files$associated_gene) %in% toupper(genes),"isoform"]
  counts = counts[counts[[1]] %in% TargetIsoforms,]
  
  write.table(counts, output_file, quote = F, row.names = F,sep = "\t")
}
