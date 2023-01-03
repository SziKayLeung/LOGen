#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose: Add columns to datasets relating to phenotype information for downstream purposes 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## classify_genotype
##


## ------------------- classify_genotype

# Aim: classify the genotype based on the sample ID 
# Input:
  # sample = str: sample id
  # case_names = vec: sample ids of cases 
  # control_names= vec: sample ids of control_names 
# Output:
  # vector of genotypes matching to the input sample ID

classify_genotype <- function(sample, case_names, control_names){
  classified_sample <- if(sample %in% case_names){label_name("case")
  } else if (sample %in% control_names){label_name("control")
  } else {"J20"
  }
  
  classified_sample <- unlist(classified_sample)
  return(classified_sample)
}