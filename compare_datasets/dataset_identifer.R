#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: identify the dataset (Iso-Seq or ONT) from isoform ID
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## identify_dataset
## identify_dataset_specific_name
##
## ---------- Notes -----------------
##
## Functions are relevant for merged datasets whereby ID has been modified to include both IDs 
## Isoforms with ID:
##  contains "_" = Iso-Seq and ONT derived isoforms
##  starting with PB = Iso-Seq derived isoforms
##  starting with TALON or ENS = ONT-derived isoforms 
## 


## ------------------- identify_dataset

# Aim: identify the source of the dataset based on the isoform ID
# Input:  
  # isoform = str: isoform id
# Output:
  # dataset = str: Iso-Seq only, Both, or ONT only

identify_dataset <- function(isoform){

  if(grepl("_", isoform)){return("Both")
  }else if(grepl("PB", isoform)){return("Iso-Seq")
  }else if (grepl("TALON", isoform)){return("ONT")
  }else if (grepl("ENS", isoform)){return("ONT")
  }else{print("NA")}

}

identify_dataset_by_counts <- function(col1,col2,name1,name2){
  
  col1 = as.numeric(col1)
  col2 = as.numeric(col2)

  if(col1 > 0 & col2 > 0){return("Both")
  }else if(col1 == 0 & col2 > 0){return(name2)
  }else if(col1 > 0 & col2 == 0){return(name1)
  }else{return("NA")}
  
}

## ------------------- identify_dataset_specific_name

# Aim: create a column with the specific ID isoform ID name (extracting from merged dataset)
  # if Dataset == ONT, then ONT_isoform column would be the isoform ID, and the IsoSeq_isoform column would be NA
  # if Dataset == Iso-Seq, then Iso-Seq isoform column would be the isoform ID, and the ONT isoform column would be NA
  # if Dataset == Both, the extract the 1st or 2nd name (separated by _) based on the column 
# Note:
  # Assumes that the ONT isoform name goes before the Iso-Seq isoform name
# Pre-requisite:
  # Run the identify_dataset function to generate a "Dataset" column in the classification file
# Input:
  # Dataset = str: Dataset column in class.files <Iso-Seq/ONT/Both>
  # isoform = str: isoform column ID <TALONXX_PBXX/TALONXX/PBXX>
  # yes_cate = str: Dataset to extract 
  # no_cate = str: Dataset not to extract
# Output: str of the isoform ID 

identify_dataset_specific_name <- function(Dataset, isoform, yes_cate, no_cate){
  
  if(Dataset == yes_cate){return(isoform)
  }else if(Dataset == no_cate){return("NA")
  }else if(Dataset == "Both" & yes_cate == "ONT"){return(word(isoform,c(2), sep = fixed("_")))
  }else if(Dataset == "Both" & yes_cate == "Iso-Seq"){return(word(isoform,c(1), sep = fixed("_")))
  }else{print("ERROR")}

}


## ------------------- wrapper_identify_dataset

# Aim: wrapper to create 3 columns in a merged dataset of ONT and Iso-Seq to determine if the isoforms are unique or common in both
  # and to extract the Iso-Seq and ONT-derived specific isoform if detected in both/unique
# Input:
  # class.files = df: classification file from SQANTI
# Output:
  # class.files = df with 3 additional column, <Dataset, IsoSeq_isoform, ONT_isoform>

wrapper_identify_dataset <- function(class.files){
  
  # identify the dataset based on the isoform ID
  class.files$Dataset <- sapply(class.files$isoform, function(x) identify_dataset(x))
  
  # extract the Iso-Seq isoform or ONT isoform ID based on the "Dataset" column
  class.files$IsoSeq_isoform <- apply(class.files, 1, function(x) identify_dataset_specific_name (x[["Dataset"]], x[["isoform"]], "Iso-Seq","ONT"))
  class.files$ONT_isoform <- apply(class.files, 1, function(x) identify_dataset_specific_name (x[["Dataset"]], x[["isoform"]], "ONT","Iso-Seq"))
  
  return(class.files)
}