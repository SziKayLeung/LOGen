#!/usr/bin/env Rscript
# Szi Kay Leung 
# 09/08/2021: reclassify the merged transcripts from tama merge back to the original transcripts from targeted and whole transcriptome relating to target genes

## Input Files ##
# 1. output file from Tama merge (WholeTargeted)
# 2. SQANTI classification file from targeted transcriptome (file used for tama merge)
# 3. SQANTI classification file from whole transcriptome (file used for tama merge)
# 4. SQANTI classification of merged tama file 

## Output Files ## 
# 1. merged_wholePBID.txt: annotation of merged IDs from tama merge with original source from whole transcriptome
# 2. merged_targetedPBID.txt: annotation of merged IDs from tama merge with original source from targeted transcriptome
# 3. merged_whole_targeted_genesID.txt: inclusion list of transcripts from target genes from tama merge sqanti classification file, for tappAS input 

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")

############# Read Files 
root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/"
whole_root_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/"
#1. output file from TAMA merge 
WholeTargeted <- read.table(paste0(root_dir,"Whole_Targeted/merged_whole_targeted_trans_report.txt"), header = T)

# number of sources from collapsed transcript (by separating all_source_trans by comma and counting occurence )
WholeTargeted = WholeTargeted %>% mutate(source_num = nchar(gsub('[^,]', '', WholeTargeted$all_source_trans))+1)
#WholeTargeted %>% filter(source_num == "2") %>% nrow() == WholeTargeted[WholeTargeted$sources == "Targeted,Whole",] %>% nrow()
#nrow(WholeTargeted[WholeTargeted$sources == "Whole" & WholeTargeted$source_num != "1",]) == 0
#nrow(WholeTargeted[WholeTargeted$sources == "Targeted" & WholeTargeted$source_num != "1",]) == 0

# Create new column with either "Whole", "Targeted" or "Both" depending on the number of sources 
WholeTargeted = WholeTargeted %>% mutate(isoform = ifelse(source_num == "1", word(all_source_trans,c(2), sep = fixed("_")),"Both"))

# 2. SQANTI classification file from targeted transcriptome (file used for tama merge)
targeted.class.names.files = paste0(root_dir,"/COLLAPSE_FILTER/AllMouseTargeted_sqantisubset.classification.txt")
targeted.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

# 3. SQANTI classification file from whole transcriptome (file used for tama merge)
whole.class.names.files = paste0(whole_root_dir,"/COLLAPSE_FILTER/WholeIsoSeq_sqantisubset.classification.txt")
whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

# 4. SQANTI classification of merged tama file 
merged.class.files = read.table(paste0(root_dir,"Whole_Targeted/merged_whole_targeted_classification.txt"), header = T)


############# Annotate merged isoforms with relevant classification files
WholeTargetedBoth <- WholeTargeted %>% filter(source_num == "2")

### Split the merged isoforms by the source transcripts from the two sequencing approaches
# split the all_source_trans with the targeted and whole id into two columns, to extract the relevant source ID depending on the sequencing approach
# note the order is not always the same i.e not always targeted<pbid>,whole<pbid>; therefore need to split into two separate columns and do downstream processing
WholeTargetedBoth_split <- data.frame(do.call('rbind', strsplit(as.character(WholeTargetedBoth$all_source_trans),',',fixed=TRUE))) %>% cbind(WholeTargetedBoth$transcript_id)

# separate_mergedtranscript <Whole/Targeted>
# extract the PBID of merged transcript by source technology from the two columns of the WholeTargetBoth_split
separate_mergedtranscript <- function(type){
  dat = list()
  # loop through each row, and if the type "Targeted" or "Whole" is in column 1, extract the Pb.ID from column 1, otherwise extract from column 2
  for(i in 1:nrow(WholeTargetedBoth_split)){
    dat[[i]] = ifelse(grepl(type,WholeTargetedBoth_split$X1[i]) == "TRUE", 
                                word(WholeTargetedBoth_split$X1[i],c(2),sep = fixed("_")),
                                word(WholeTargetedBoth_split$X2[i],c(2),sep = fixed("_")))
  }
  return(dat)
}

WholeBothPBID = data.frame(transcript_id = WholeTargetedBoth$transcript_id,isoform = separate_mergedtranscript("Whole") %>% do.call(rbind, .), sources = "Both") 
TargetedBothPBID = data.frame(transcript_id = WholeTargetedBoth$transcript_id,isoform = separate_mergedtranscript("Targeted") %>% do.call(rbind, .), sources = "Both")

### Annotate isoforms the source transcripts by associated_gene and transcript using the relevant sqanti classification files (see above)  
Whole_anno = 
  # Filter those transcripts that are only identified in Whole Transcriptome
  WholeTargeted %>% filter(sources == "Whole") %>% select(transcript_id,isoform, sources) %>% 
  # Combined with the merged transcripts with Whole Transcriptome PB.ID 
  rbind(.,WholeBothPBID) %>% 
  # Join with the relevant classification file
  full_join(.,whole.class.files[,c("isoform","associated_transcript","associated_gene")])

# Apply the same steps to filtering transcripts only from the targeted transcriptome
Targeted_anno = WholeTargeted %>% filter(sources == "Targeted") %>% select(transcript_id,isoform, sources) %>% rbind(.,TargetedBothPBID) %>% full_join(.,targeted.class.files[,c("isoform","associated_transcript","associated_gene")])

# transcripts that were detected in both should only be the target genes
#Whole_anno[Whole_anno$sources == "Both","associated_gene"] 
#Targeted_anno[Targeted_anno$sources == "Both","associated_gene]

### Output Files 
# create files for the corresponding tama_merge IDs to the original source PB.ID from targeted and whole transcriptome 
output_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/Whole_Targeted/"
write.table(Whole_anno,paste0(output_dir,"/merged_wholePBID.txt"),quote = F, row.names = F)
write.table(Targeted_anno,paste0(output_dir,"/merged_targetedPBID.txt"),quote = F, row.names = F)


############# Subset Isoforms of associated target genes for tappAS 
# create an inclusion transcript list of the isoforms associated with target genes for input to tappAS for RNA-Seq Input
# Speeds up process as not interested in off-target genes
targetgenemerged_ID = data.frame(merged.class.files[toupper(merged.class.files$associated_gene) %in% TargetGene,"isoform"])
write.table(targetgenemerged_ID,paste0(root_dir,"Whole_Targeted/merged_whole_targeted_genesID.txt"), quote = F, col.names = F, row.names = F)
