#!/usr/bin/env Rscript
# Szi Kay Leung 29/06/2021: Whole Transcriptome - Generate files of lists of isoforms present in each group (for downsteam processing - WT2mos vs WT8mos vs TG2mos vs TG8mos)
# List of isoforms present in each group is determined by the counts in the sqanti file 
# Rscript script.R <class.names.files> <output.counts.dir>
# output: 4 files <group>_counts.txt

suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))
suppressMessages(library("stringr"))

# arguments
args = commandArgs(trailingOnly=TRUE) 
class.names.files <- args[1] 
output.counts.dir <- args[2]
#class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq_sqantitamafiltered.classification.txt"

cat("Processing:", class.names.files,"\n")

# read the classification file 
class.files <- read.table(class.names.files, header = T, sep = "\t", as.is = T) 
rownames(class.files) <- class.files$isoform

# Subset the columns from the classification files with the expression 
# Extract the column names, with only the the phenotype and the age (remove FL.<sample>)
counts = class.files %>% select(starts_with("FL."))
colnames(counts) = substr(colnames(counts),8,14)

# Sum the expression counts based on genotype and age
count_simple = as.data.frame(sapply(unique(colnames(counts)), function(x) rowSums(counts[,grepl(x, colnames(counts))])))
count_simple = count_simple %>% rownames_to_column(., var = "transcript_id")

# generate separate files for PBIDs present in group 
# present PBID based on count > 0 
cat("Generating Separate Files for each group in:", class.names.files,"\n")
group_order = c("WT_2mos","TG_2mos","WT_8mos","TG_8mos")
for(group in group_order){
  dat = count_simple[count_simple[[group]] != 0 ,"transcript_id"] %>% as.data.frame()
  write.table(dat,paste0(output.counts.dir,"/",group,"_counts.txt"), row.names = F, quote = F, col.names = F)
  #nrow(count_simple[count_simple$WT_2mos != 0,])
}

# generate separate files for genotype (WT vs TG)
counts = class.files %>% select(starts_with("FL."))
colnames(counts) = substr(colnames(counts),8,9)
count_simple = as.data.frame(sapply(unique(colnames(counts)), function(x) rowSums(counts[,grepl(x, colnames(counts))])))
count_simple = count_simple %>% rownames_to_column(., var = "transcript_id")
for(group in c("WT","TG")){
  dat = count_simple[count_simple[[group]] != 0 ,"transcript_id"] %>% as.data.frame()
  write.table(dat,paste0(output.counts.dir,"/",group,"_counts.txt"), row.names = F, quote = F, col.names = F)
  #nrow(count_simple[count_simple$WT != 0,])
}

cat("Generating Separate Files for each sample in:", class.names.files,"\n")
# generate separate files for PBIDs present in each sample
abundance =  class.files %>% select(starts_with("FL."))
for(col in colnames(abundance)){
  
  file = word(col,c(2), sep = fixed("."))
  cat("Processing",file,"\n")
  dat = rownames(abundance[abundance[[col]] != "0",]) %>% as.data.frame() 
  
  write.table(dat,paste0(output.counts.dir,"/",file,"_counts.txt"), row.names = F, quote = F, col.names = F)
  #nrow(abundance[abundance$FL.S23_WT_8mos != "0",])
}


