#!/usr/bin/env Rscript

## ---------- packages -----------------

suppressMessages(library("data.table"))
suppressMessages(library("optparse"))

LOGEN <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))

## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-c", "--classfile"), type="character", default=NULL, 
              help="SQANTI classification file", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression file csv"),
  make_option(c("-s", "--subset"), type="character", default=NULL, 
              help="Subsetted dataset by column name", metavar="character"),
  make_option(c("--monoexonic"), action="store_true", default=TRUE, help="removing mono-exonic intergenic, genic transcripts"),            
  make_option(c("--nsample"), type="integer", default=2L, help="filtering: minimum number of samples [default %default]", metavar="number"),
  make_option(c("--nreads"), type="integer", default=2L, help="filtering: minimum number of reads [default %default]", metavar="number"),
  make_option(c("--dir"), type="character", default=NULL, help="output directory, default = directory of classification file"),
  make_option(c("--name"), type="character", default=NULL, help="output name, default = name of classification file + subsetted name")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$subset)){
  stop("Need subset parameter")
}

if(is.null(opt$dir)){
  opt$dir = dirname(opt$classfile)
}

if(is.null(opt$name)){
  opt$name = sub("\\.txt$", "", basename(opt$classfile))
}



## ---------- Input files -----------------

classfile = SQANTI_class_preparation(opt$classfile,"ns")

message("Reading",opt$expression)
exp = fread(opt$expression)

if (!any(grepl(opt$subset, colnames(exp)))){
  stop("Subset parameter not in expression columns")
}

if(isTRUE(opt$monoexonic)){
  message("Removing monoexonic intergenic transcripts")
  message("Number of transcripts before filtering mono-exonic: ", nrow(classfile))
  monoexonic <- classfile[classfile$subcategory == "mono-exon" & classfile$structural_category != "full-splice_match","isoform"]
  classfile <- classfile[!classfile$isoform %in% monoexonic,]
  message("Number of transcripts after filtering mono-exonic: ", nrow(classfile))
}else{
  message("Not removing monoexonic transcripts")
}


## ---------- Functions -----------------

# Define a custom function to count non-zero elements in a row
count_non_zero <- function(row) {
  return(sum(row != 0))
}


filter_dataset <- function(subcolnames){
  
  subExp = exp %>% select(contains(subcolnames))
  message("Keeping the following columns for filtering")
  print(colnames(subExp))
  
  nreads <- apply(subExp,1,sum)
  nsample <- apply(subExp, 1, count_non_zero)
  
  subExp$nreads <- nreads
  subExp$nsample <- nsample
  
  subExp$isoform <- exp$id
  message("Filtering by ", opt$nreads, " reads and ", opt$nsample, " samples")
  filterSubExp <- subExp %>% filter(nreads >= opt$nreads, nsample >= opt$nsample)
  
  return(filterSubExp)
}

filterSubExp <- filter_dataset(opt$subset)
subclassfile <- merge(classfile[,1:15],filterSubExp,by="isoform")

message("Output:",opt$dir,"/",opt$name,"_",opt$subset,".txt")
write.table(subclassfile,paste0(opt$dir,"/",opt$name,"_",opt$subset,".txt"),quote=F,sep="\t",row.names=F)
write.table(subclassfile$isoform,paste0(opt$dir,"/",opt$name,"_",opt$subset,"ID.txt"),quote=F,sep="\t",col.names=F,row.names=F)

