#!/usr/bin/env Rscript
# Szi Kay Leung 29/06/2021: Generate RNA-Seq expression file for TappAS after alignment of individual RNA-Seq samples with Kallisto 
# Rscript script.R <input.dir> <output.file> <type=Whole/Targeted/WholeTargeted> <targeted.class.files>
# Rscript script.R </gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/ALLRNASEQ> <WholeMouseRNASeq_sqantitamafiltered.expression.txt> <Whole>

# Prequisite: Generate Kallisto files in one directory with alignment to Iso-Seq scaffold 
# Different output depending on whether whole or targeted analysis 
  # if targeted analysis, further subsetting of counts abundance file to only include expression of isoforms associated with Target Genes (ease tappAS analysis)

## Arguments ###
# input.dir = Directory containing all the subdirectories of results with "abundance.tsv" from running Kallisto 
# output.file = Name of the output file (note this will be generated in the input dir)
# type = "Whole" or "Targeted"
# targeted.class.files = targeted sqanti classification files for subsetting isoforms associated to the target genes if type == "Targeted" or merged sqanti classification file (whole + targeted transcriptome)

suppressMessages(library("stringr"))
suppressMessages(library("tidyr"))

# arguments
args = commandArgs(trailingOnly=TRUE) 
rnaseq_dir <- args[1] 
output.file <- args[2]
type <- args[3]
class.files <- args[4]

cat("Processing files in:", rnaseq_dir,"\n")
cat("Analysis:", type, "\n")

# list all the files in the rnaseq dir with the pattern "abundance.tsv"
# samples = the folder/sample from which the abundance.tsv came from

abundance = lapply(list.files(path = rnaseq_dir, recursive = T, include.dirs = T, pattern = "abundance.tsv", full.names = T),
                   function(x) read.table(x,header = T))
samples = word(list.files(path = rnaseq_dir, recursive = T, include.dirs = T, pattern = "abundance.tsv"),c(1),sep=fixed("/"))

# only extract the relevant columns from each file 
mod_abundance = lapply(abundance, function(x) x[,c("target_id","est_counts")])
names(mod_abundance) = samples

# bind all the files together, with the sample as differentiator 
all = dplyr::bind_rows(mod_abundance, .id = "sample")

# spread from long to wide in the format required for tappAS 
counts = all %>% spread(sample, est_counts)
colnames(counts)[1] <- ""

# Output 
print(paste0("Writing output file:", rnaseq_dir,"/", output.file))
if(type == "Whole"){
  write.table(counts, paste0(rnaseq_dir,"/",output.file), quote = F, row.names = F,sep = "\t")
    
}else if(type %in% c("Targeted","WholeTargeted")){
  # put "FL" in front of the same columns for same formatting as the phenotype file
  colnames(counts)[2:length(colnames(counts))] = paste0("FL.",colnames(counts)[2:length(colnames(counts))])
  
  # required for subsetting on the isoforms associated with target genes
  TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")
  targeted.class.names.files = read.table(class.files, as.is = T, sep = "\t", header = T)
  
  # subset on the isoforms associated with target genes and the counts 
  TargetIsoforms <- targeted.class.names.files[toupper(targeted.class.names.files$associated_gene) %in% TargetGene,"isoform"]
  counts = counts[counts[[1]] %in% TargetIsoforms,]
  
  write.table(counts, paste0(rnaseq_dir,"/", output.file), quote = F, row.names = F,sep = "\t")
}else{
  print("Targeted or Whole required as 3rd argument")
}



