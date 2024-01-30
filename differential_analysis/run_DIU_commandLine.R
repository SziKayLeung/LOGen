## ---------- Script -----------------
##
## Purpose: perform differential isoform usage analysis on mouse rTg4510 ONT and Iso-Seq targeted datasets
## Adapted tappAS DIU analysis scripts
## https://github.com/ConesaLab/tappAS/blob/master/scripts/DIU.R
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)


## ---------- packages -----------------

# Specify the libraries to load
libraries_to_load <- c("data.table", "stringr", "dplyr", "ggplot2", "optparse", "here")

# Load or install each library
for (library_name in libraries_to_load) {
  if (!requireNamespace(library_name, quietly = TRUE)) {
    install.packages(library_name, dependencies = TRUE)
  }
  suppressMessages(library(library_name, character.only = TRUE))
}


## ---------- source functions -----------------

LOGEN_ROOT ="/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen"
source(paste0(LOGEN_ROOT, "/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN_ROOT, "/differential_analysis/plot_usage.R"))


## ---------- arguments -----------------

option_list <- list( 
  make_option(c("-c", "--classfile"), type="character", default=NULL, 
              help="Input sqanti classification file", metavar="character"),
  make_option(c("-e", "--exp"), type="character", default=NULL, 
              help="Input normalised expression file"),
  make_option(c("-p", "--phenotype"), type="character", default=NULL, 
              help="Input phenotype file", metavar="character"),
  make_option(c("-f", "--factor"), type="character", default=NULL, 
              help="Input factor file", metavar="character"),
  make_option(c("-t", "--typecharacter"), type="character", default=NULL, 
              help="common character in expression columns to capture from classification file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output_directory", metavar="character"),
  make_option(c("-n", "--output_name"), type="character", default=NULL, 
              help="Output name, default: sqanti prefix", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$output_name)){
  opt$output_name=word(basename(opt$classfile),c(1),sep=fixed("_")) 
}

if(is.null(opt$output_dir)){
  opt$output_dir=paste0(dirname(opt$classfile),"/")
}


## ---------- Input files -----------------

# classification files
message("Read in: ",opt$classfile)
class.files <- read.table(opt$classfile,sep="\t",as.is=T,header=T)
row.names(class.files) <- class.files$isoform
#head(class.files)

# phenotype
message("Read in: ",opt$phenotype)
phenotype <- read.csv(opt$phenotype) 
#print(phenotype)

# DESEQ results
message("Processing: ",opt$exp)
Exp <- list(targeted = as.data.frame(fread(opt$exp)))
Exp$targeted <- merge(Exp$targeted, phenotype, by="sample")
Exp$targeted <- merge(Exp$targeted, class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], by = "isoform", all.x = T)

# expression
ExpMod <- list(
  ontTargTranRaw = class.files %>% dplyr::select(associated_gene, contains(opt$typecharacter)),
  ontTargTranNorm = Exp$targeted %>%
    dplyr::select(sample,normalised_counts, isoform) %>% tidyr::spread(., sample, value = normalised_counts) %>%
    tibble::remove_rownames(.) %>% tibble::column_to_rownames(var="isoform")
)
# +1 to raw data to run spliceVariant.DS and calculate calcNormFactors
ExpMod$ontTargTranRaw <- cbind(ExpMod$ontTargTranRaw[1],ExpMod$ontTargTranRaw[-1] + 1)

# factors
message("Read in: ",opt$factor,"\n")
factorsInput = read.table(opt$factor)
head(factorsInput)


## ---------- Process DIU -----------------

resultsDIU <- runDIU(transMatrixRaw=ExpMod$ontTargTranRaw,transMatrix=ExpMod$ontTargTranNorm,classf=class.files,
                        myfactors=factorsInput,filteringType="FOLD",filterFC=2)

if(!is.null(resultsDIU)){
  message("DIU successful")
  message("Writing output from DIU to:", paste0(opt$output_dir,opt$output_name,"_resultDIU.txt"))
  write.table(resultsDIU$resultDIU,paste0(opt$output_dir,opt$output_name,"_resultDIU.txt"),quote=F,row.names=F,sep="\t")
  
  message("Writing output of kept isoforms to:", paste0(opt$output_dir,opt$output_name,"_keptIso.txt"))
  write.table(resultsDIU$keptIso,paste0(opt$output_dir,opt$output_name,"_keptIso.txt"),sep="\t",quote=F)
  
}else{
  message("DIU did not run")
}

