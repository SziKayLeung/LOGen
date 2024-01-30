## ---------- Script -----------------
##
## Script name: merge_talon_sqanti_forcounts.R
##
## Purpose of script: Include the counts from the TALON abundance file into the sqanti classification file for downstream processing
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
## Requires a metadata phenotype file (Tab delimitied), with 4 columns: Sample, Condition, ID, SQ_ID
##    SQ_ID has to match the column names of the samples in the TALON abundance file
## Output:
##    Generates a sqanti classification file with the counts from TALON matched by isoform ID


## ---------- packages -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))
suppressMessages(library("optparse"))
suppressMessages(library("tools"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-f", "--file"), type="character", 
              help="Input sqanti classification file", metavar="character"),
  make_option(c("-c", "--counts"), type="character", 
              help="abundance file as csv or txt"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Metadata file required if --talon=TRUE;#Sample#Condition(CasevsControl)#SQ_ID", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output_directory", metavar="character"),
  make_option(c("-t", "--talon"), type="character", default=TRUE, 
              help="TRUE/FALSE if abundance file is generated through TALON", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ---------- Input files -----------------

cat("Processing",opt$file,"\n")
class.files <- read.table(opt$file, as.is = T, sep = "\t", header = T) 

cat("Processing",opt$counts,"\n")
if(file_ext(opt$counts) == "txt"){
  counts <- read.table(opt$counts, as.is = T, header = T, sep = "\t")  
}else{
  counts <- read.csv(opt$counts, header = T)  
}
head(counts)


## ---------- Merge classification and counts file -----------------

if(opt$talon == TRUE){
  cat("Processing files generated from TALON\n")
  if(is.null(opt$metadata)){
    cat("Metadata argument required; -m\n")
    quit(status=1)
  }else{
    cat("Processing",opt$metadata,"\n")
    metadata <- read.table(opt$metadata, as.is = T, header = T, sep = "\t")
    
    # subset the TALON counts file by the counts of the samples of interest
    counts <- counts %>% select(annot_transcript_id, metadata$SQ_ID) %>% rename(isoform = annot_transcript_id)
  }
}

# merge the original classification file with the samples of interest
# only keep the counts from the isoforms in the classification file (as final, filtered file)
class.files.merged <- merge(class.files, counts, by = "isoform", all.x = T)

# write output file
basename = tools::file_path_sans_ext(basename(opt$file))
write.table(class.files.merged, paste0(opt$output_dir,"/",basename,"_counts.txt"), sep = "\t", quote = F, row.names = F)
