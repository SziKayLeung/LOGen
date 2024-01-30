# Szi Kay Leung
# 26/07/2021 
# Aim: Generate collapsed expression file using SQANTI filtered dataset and SQANTI-TAMA filtered datasets 
# Applicable for whole and targeted transcriptome 

### Issue ####
# TAMA filtering removes a lot of partial transcripts and shorter FSM transcripts that may be higher in expression that the longest FSM transcript
# Therefore, combine the counts from FSM and ISM transcripts that are annotated to the same transcript ensembl ID 

#### Aims ####
# 1. Using SQANTI filtered dataset, subset on known isoforms using associated_transcript collumn and collapse FL reads across each sample
  # Therefore collapsing ISM transcripts and other FSM transcripts that are annotated to the same ENSEMBL transcript ID
# 2. Annotate the ENSEMBL transcript ID with most appropriate PB.ID from the SQANTI-TAMA filtered dataset
  # Important as require final annotation for RNA-Seq alignment downstream and for tappas expression file 
  # SQANTI-TAMA filtered dataset are transcripts that have passed TAMA filtering, and are therefore unique 
  # 2a. Annotate the ENSEBML transcript ID with FSM isoform (PB.ID) if available 
  # 2b. Otherwise use ISM isoform 
  # 2c. If no FSM isoforms, but more than two ISM isoforms, use the longest ISM isoform as reference 
  # Note: Only the isoforms that are retained in tama filtering are used for reference and kept for expression
# 3. Keep the novel transcripts as they are and do not collapse reads
# 4. Further, refine transcripts that are filtered by TAMA filtering as many redundant ISM transcripts that are annotated to the same emsembl transcript Id, but are retained due to internal fragment, 3prime fragments (check if these transcripts are also removed by IsoPop)

#### Output ####
# 1 & 2. Long Read expression file for downstream analysis: All the counts are accounted for from SQANTI filtering but only the longest unique isoforms are used as annotation 
  # 1. Used Length (i.e. longest transcript) to identify most appropriate FSM isoform (in conjunction with TAMA filtering) 
  # 2. Used Expression (i.e. most abundant transcript) to identify most appropriate FSM isoform (only SQANTI filtering)
# 3. List of the longest, unique isoforms for further subsetting
# 4. List of the longest, most abundant FSM isoforms for further subsetting


#### Usage ####
# Script <sqanti_filered_inputfile> <tama_filtered_inputfile> <outupt_expression_file>.txt <output_finalisoform_file> <type == "Targeted/Whole">
args = commandArgs(trailingOnly=TRUE) 
sqanti_filtered_inputfile <- args[1] 
tama_filtered_inputfile <- args[2]
output_expression_file <- args[3]
output_finalisoform_file <- args[4]
type <- args[5]

#DiffAnalysis_WKD="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis"
#dataset= "AllMouseTargeted"
#sqname="AllMouseTargeted.collapsed_classification.filtered_lite"
#sqanti_filtered_inputfile <- paste0(DiffAnalysis_WKD,"/SQANTI3/",sqname,"_classification.txt")
#tama_filtered_inputfile <- paste0(DiffAnalysis_WKD,"/SQANTI_TAMA_FILTER/",dataset,"_sqantitamafiltered.classification.txt")
#output_expression_file <- paste0(DiffAnalysis_WKD,"/TAPPAS_INPUT/IsoSeq_Expression/AllMouseTargeted_sqantisubset.expression")
#output_finalisoform_file <- paste0(DiffAnalysis_WKD,"/TAPPAS_INPUT/Retained_collapsed_pbid")
#type <- "Targeted"

# SQANTI filtered file
#sqanti.filtered.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification#.filtered_lite_classification.txt"
#sqantifil.class.files <- read.table(sqanti.filtered.names.files, header = T)

# SQANTI, TAMA filtered file of targeted transcriptome
#targeted.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt"
#tamafil.class.files <- SQANTI_class_preparation(targeted.class.names.files,"standard")

# SQANTI, TAMA filtered file of whole transcriptome
#whole.class.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantitamafiltered.classification.txt"
#whole.class.files <- SQANTI_class_preparation(whole.class.names.files,"standard")

#output_expression_file = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_INPUT/IsoSeq_Expression/AllMouseTargeted_sqantisubset.expression.txt"

#### Library ####
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))


#### Input Files ####
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")


# Input files 
sqantifil.class.files <- SQANTI_class_preparation(sqanti_filtered_inputfile,"standard")


#### Functions & Script ####
#1. Collapse all the counts for known transcripts from the sqanti filtered classification file
# This includes ISM that would have been considered unique by TAMA filtering with 3'fragment or internal fragment however all the junctions all the same suggesting truncation

if(type == "Targeted"){
  cat("Filtering known transcripts by Targeted Genes \n")
  allcounts_known = sqantifil.class.files %>% filter(toupper(associated_gene) %in% TargetGene)
}else{
  allcounts_known = sqantifil.class.files
}

allcounts_known = allcounts_known %>% filter(associated_transcript != "novel") %>% select(associated_transcript,starts_with("FL.")) %>% group_by(associated_transcript) %>% summarise_each(funs(sum))

#sum(sqantifil.class.files[sqantifil.class.files$associated_transcript == "ENSMUST00000022616.13","FL.K17"]) == allcounts_known[allcounts_known$associated_transcript == "ENSMUST00000022616.13","FL.K17"]

## Function to find the most appropriate FSM for annotation 

# Create threshold to use most highly expressed FSM as annotation for RNA-Seq
calc = sqantifil.class.files %>% select(starts_with("FL.")) %>% mutate(n_detected = rowSums(. != 0), total = rowSums(.))
calc_median = sqantifil.class.files %>% select(starts_with("FL.")) %>% rowwise() %>% mutate(med = median(c_across(where(is.numeric)), na.rm = TRUE))
sqantifil.class.files = cbind(sqantifil.class.files, calc[,c("n_detected","total")], calc_median[,c("med")]) 

FSM_select <- function(type){
  
  #2. Annotate the ENSEMBL transcript ID with most appropriate PB.ID 
  if( type == "length_tama"){
    
    tamafil.class.files <- SQANTI_class_preparation(tama_filtered_inputfile,"standard")
    input_class.files  = tamafil.class.files
    
    # ISM that were kept in after tama filter and likely to be unique transcripts 
    # Tama also removed transcripts with different names (all ISM), however likely to be partial transcripts of a different BLAST
    oneFSM = list()                 # list of ENMUST transcript ID that only have one FSM 
    moreFSM = list()                # list of ENMUST transcript ID that with more than one FSM
    noFSM = list()                  # list of ENMUST transcript ID that with no FSM, i.e. only have ISM
    tamaremoved_oneFSM = list()     # list of ENMUST transcript ID that were FSM but filtered out by TAMA collapse
    tamaremoved_noFSM = list()      # list of ENMUST transcript ID that were ISM but filtered out by TAMA collapse
    
    #2a. looping through the known transcript ID to find the structural_category 
    for(trans in allcounts_known$associated_transcript){
      cate = tamafil.class.files[tamafil.class.files$associated_transcript == trans,"structural_category"]
      
      # if transcript is kept from TAMA filtering
      # discern if there is 1 FSM, more than one FSM, or no FSM for that transcript
      if(length(cate) > 0){
        if('FSM' %in% cate){
          if(length(cate[cate == "FSM"]) == 1){oneFSM[[trans]] = trans 
          }else{moreFSM[[trans]] = trans}
        }else{
          noFSM[[trans]] = trans}
        
        # if transcript is not kept from TAMA filtering
        # discern structural category from the SQANTI filtered dataset before filtering from TAMA
      }else{
        cate = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,"structural_category"]
        if(length(cate[cate == "FSM"]) > 1){
          tamaremoved_oneFSM[[trans]] = trans 
        }else{
          tamaremoved_noFSM[[trans]] = trans  
        }
      }
    }
    
    # 2b. For ENSEMBL Transcripts with more than one FSM, keep the longest
    moreFSM_highlyexpressed = list()
    for(trans in row.names(do.call("rbind", moreFSM))){
      # keep the longest ISM transcript as annotation 
      pbid = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,] %>% .[which.max(.$length),"isoform"]
      moreFSM_highlyexpressed[[trans]] = pbid
    }
    
  }else if(type == "length"){
    
    input_class.files = sqantifil.class.files
    
    #2a. Use the most highly expressed FSM transcript as annotation rather than the longest 
    oneFSM = list()                 # list of ENMUST transcript ID that only have one FSM 
    moreFSM = list()                # list of ENMUST transcript ID that with more than one FSM
    noFSM = list()                  # list of ENMUST transcript ID that with no FSM, i.e. only have ISM
    
    for(trans in allcounts_known$associated_transcript){
      cate = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,"structural_category"]
      
      if('FSM' %in% cate){
        if(length(cate[cate == "FSM"]) == 1){oneFSM[[trans]] = trans 
        }else{moreFSM[[trans]] = trans}
      }else{
        noFSM[[trans]] = trans}
    }
    
    # 2b. For ENSEMBL Transcripts with more than one FSM, keep the most abundantly expressed FSM
    moreFSM_highlyexpressed = list()
    for(trans in row.names(do.call("rbind", moreFSM))){
      # keep the longest ISM transcript as annotation 
      pbid = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,] %>% filter(structural_category == "FSM") %>% 
        .[which.max(.$length),"isoform"]
      moreFSM_highlyexpressed[[trans]] = pbid
    }
    
  }else if(type == "expression"){
    
    input_class.files = sqantifil.class.files
    
    #2a. Use the most highly expressed FSM transcript as annotation rather than the longest 
    oneFSM = list()                 # list of ENMUST transcript ID that only have one FSM 
    moreFSM = list()                # list of ENMUST transcript ID that with more than one FSM
    noFSM = list()                  # list of ENMUST transcript ID that with no FSM, i.e. only have ISM
    
    for(trans in allcounts_known$associated_transcript){
      cate = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,"structural_category"]
      
      if('FSM' %in% cate){
        if(length(cate[cate == "FSM"]) == 1){oneFSM[[trans]] = trans 
        }else{moreFSM[[trans]] = trans}
      }else{
        noFSM[[trans]] = trans}
    }
    
    # 2b. For ENSEMBL Transcripts with more than one FSM, keep the most abundantly expressed 
    moreFSM_highlyexpressed = list()
    for(trans in row.names(do.call("rbind", moreFSM))){
      # keep the longest ISM transcript as annotation 
      pbid = sqantifil.class.files[sqantifil.class.files$associated_transcript == trans,] %>% .[.$n_detected >= 3 & .$total >= 5,] %>% .[which.max(.$med),"isoform"]
      moreFSM_highlyexpressed[[trans]] = pbid
    }
    
    
  }else{print("type either as length of expression")}
  
  cat("Number of known transcripts with 1 FSM:", length(oneFSM),"\n")
  cat("Number of known transcripts with more than one FSM:", length(moreFSM),"\n")
  cat("Number of known transcripts with only ISM:", length(noFSM),"\n")
  
  
  # Representative Transcripts for ENMUST where there is only one FSM transcript kept in by TAMA filtering
  oneFSM_known = input_class.files[input_class.files$associated_transcript %in% row.names(do.call("rbind", oneFSM)) & input_class.files$structural_category == "FSM",]
  
  # Representative Transcripts for ENMUST where there is multiple FSM transcripts 
  moreFSM_known = input_class.files[input_class.files$isoform %in% do.call("rbind", moreFSM_highlyexpressed)[,1] & input_class.files$structural_category == "FSM",]
  
  
  # 2b. For ENSEBML Transcripts with no FSM, use either the ISM transcript ID if there is only one, or the longest ISM transcript if more than 1
  # note several combination of noFSM_onetrans + oneFSM_longesttrans < noFSM_known$associated_transcript due to collapsed several ISM transcripts annotated to same transcript
  noFSM_onetrans = list()
  noFSM_longesttrans = list()
  
  # looping through the known transcript ID where there is only ISM transcripts
  for(trans in row.names(do.call("rbind", noFSM))){
    cate = input_class.files[input_class.files$associated_transcript == trans,"structural_category"]
    if(length(cate) == 1){
      pbid = input_class.files[input_class.files$associated_transcript == trans,"isoform"]
      noFSM_onetrans[[trans]] = pbid
    }else{
      # keep the longest ISM transcript as annotation 
      pbid = input_class.files[input_class.files$associated_transcript == trans,] %>% .[which.max(.$length),"isoform"]
      noFSM_longesttrans[[trans]] = pbid
    }
  }
  
  # Representative Transcripts for ENMUST where there is no FSM transcript kept in by TAMA filtering
  noFSM_known = bind_rows(input_class.files[input_class.files$associated_transcript %in% row.names(do.call("rbind", noFSM_onetrans)),],
                          input_class.files[input_class.files$isoform %in% do.call("rbind", noFSM_longesttrans)[,1],])
  
  # Merge the representative PB.ID for each Transcript ENSEMBL id with the counts
  oneFSM_known_collcounts = merge(oneFSM_known[,c("isoform","associated_transcript")],allcounts_known)
  noFSM_known_collcounts = merge(noFSM_known[,c("isoform","associated_transcript")],allcounts_known)
  moreFSM_known_collcounts = merge(moreFSM_known[,c("isoform","associated_transcript")],allcounts_known)
  
  
  # Keep the novel transcripts as they are and do not collapse reads
  if(type == "Targeted"){
    cat("Filtering novel transcripts by Targeted Genes \n")
    allcounts_novel = input_class.files %>% filter(toupper(associated_gene) %in% TargetGene)
  }else{
    allcounts_novel = input_class.files
  }
  
  allcounts_novel = allcounts_novel[allcounts_novel$associated_transcript == "novel",] %>% select(isoform,starts_with("FL."))
  
  # Combine all the transcripts together (known and novel) for final expression file 
  all_counts = bind_rows(oneFSM_known_collcounts %>% select(-associated_transcript),noFSM_known_collcounts %>% select(-associated_transcript),
                         moreFSM_known_collcounts%>% select(-associated_transcript), allcounts_novel)
  colnames(all_counts)[1] = ""
  
  return(all_counts)
  #all_counts_genes = merge(all_counts,sqantifil.class.files[,c("isoform","associated_gene")])
}

all_counts_length = FSM_select("length")
all_counts_expression = FSM_select("expression")

# Generate a list of final unique 1-1 isoform:transcript ensembl ID for downstream filtering
# By Length
cat("Output classification file for longest FSM isoforms:",output_expression_file,"_FSMbylength.txt\n")
write.table(all_counts_length,paste0(output_expression_file,"_bylength.txt"),row.names = F, quote = F, sep = "\t")

cat("Output isoform list for longest FSM isoforms:",output_finalisoform_file,"_FSMbylength.txt\n")
write.table(data.frame(all_counts_length[,1]),paste0(output_finalisoform_file,"_FSMbylength.txt"), quote = F, row.names = F, col.names = F)

# By Expression
cat("Output classification file for most abundant FSM isoforms:",output_expression_file,"_FSMbyexpression.txt\n")
write.table(all_counts_expression,paste0(output_expression_file,"_byexpression.txt"),row.names = F, quote = F, sep = "\t")

cat("Output isoform list for most abundant FSM isoforms:",output_finalisoform_file,"_FSMbylexpression.txt\n")
write.table(data.frame(all_counts_expression[,1]),paste0(output_finalisoform_file,"_FSMbyexpression.txt"), quote = F, row.names = F, col.names = F)