# SQANTI filtered file
#sqanti.filtered.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
#sqantifil.targeted.classfilesss <- read.table(sqanti.filtered.names.files, header = T)

# SQANTI, TAMA filtered file of targeted transcriptome
#targeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
#targeted.classfiles <- SQANTI_class_preparation(targeted.class.names.files,"standard")


# SQANTI, TAMA filtered file of whole transcriptome
#whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantitamafiltered.classification.txt"
#whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")

collapsed.files <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt", header = T)

sqanti.files <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.txt", header = T) %>% filter(toupper(associated_gene) %in% TargetGene)

sqantifiltered.reasons <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_reasons.txt", sep = ",", header = T) %>% filter(filtered_isoform %in% sqanti.files$isoform)

sqantifiltered.files <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt", header = T) %>% filter(toupper(associated_gene) %in% TargetGene)

tama_filtered <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt", header = T) 

tama_filtered <- tama_filtered %>% filter(toupper(associated_gene) %in% TargetGene)

setdiff(tama_filtered$isoform,collapsed.files$isoform)

nrow(sqanti.files)
nrow(sqantifiltered.files)
nrow(collapsed.files)
nrow(tama_filtered)

tama_retained = setdiff(tama_filtered$isoform,collapsed.files$isoform)
setdiff(collapsed.files$isoform,tama_filtered$isoform)

# all incomplete splice match that are retained by tama but filtered for the collapsed collapsed.filesaset 
sqantifiltered.files[sqantifiltered.files$isoform %in% tama_retained,] %>% group_by(structural_category) %>% tally()

# final checks 
for(trans in tama_filtered$isoform){
  #cat("Processing",trans,"\n")
  ensembl = tama_filtered[tama_filtered$isoform == trans, "associated_transcript"]
  if(!ensembl %in% collapsed.files$associated_transcript){
    print(trans)
  }
}
