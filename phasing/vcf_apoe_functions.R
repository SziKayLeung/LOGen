#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Commands to run genotyping and phasing from deepVariant results
## --------------------------------


read_vcf <- function(inputvcf){
  print(inputvcf)
  vcf <- read.vcfR(inputvcf)
  
  if(nrow(vcf@fix) == 1){
    vcf_header_df <- as.data.frame(t(getFIX(vcf)))
  }else{
    vcf_header_df <- data.frame(getFIX(vcf))
  }
  vcf_dt <- data.frame(vcf@gt)
  vcf <- cbind(vcf_header_df, vcf_dt)
  colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","vcfGT")
  return(vcf)
}

read_gvcf <- function(gvcf){
  gvcf <- read.vcfR(gvcf)
  gvcf_dt <- data.frame(gvcf@gt)
  gvcf_header_df <- data.frame(getFIX(gvcf))
  gvcf <- cbind(gvcf_header_df, gvcf_dt)
  colnames(gvcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","gvcfGT")
  return(gvcf)
}

# obtain the genotype from deepVariant
obtain_genotype_vcf <- function(GT, VAF, DP, REF, ALT){
  
  VAF <- as.numeric(VAF)
  if (is.na(VAF) | is.na(DP)) {
    return("NA")
  }
  
  if(GT == "./."){
    if(DP == "."){
      return(paste0("NA"))
    }
    else if(VAF >= 0.4 & VAF < 0.8){
      return(paste0(REF,"/",ALT))
    }else if(VAF >= 0.8){
      return(paste0(ALT,"/",ALT))
    }else if(VAF < 0.4){
      return(paste0(REF,"/",REF))
    }else{
      return(paste0("NA"))
    }
  }else{
    if(GT == "0/1"){
      return(paste0(REF,"/",ALT))
    }else if(GT == "0/0"){
      return(paste0(REF,"/",REF))
    }else{
      return(paste0(ALT,"/",ALT))
    }
  }
}


## Functions 
genotyping <- function(ref_allele_count, COUNTED, ALT){
  ref_allele_count <- as.numeric(ref_allele_count)
  if(ref_allele_count == 0){
    return(paste0(ALT,"/",ALT))
  }else if(ref_allele_count == 1){
    return(paste0(COUNTED,"/",ALT))
  }else{
    return(paste0(COUNTED,"/",COUNTED))
  }
}

genotyping_num <- function(ref_allele, minor_allele){
  ref_allele_count <- as.numeric(ref_allele)
  minor_allele_count <- as.numeric(minor_allele)
  if(ref_allele_count > 0 & minor_allele_count == 0){
    return(paste0(0))
  }else if(ref_allele_count > 0 & minor_allele_count > 0){
    return(paste0(1))
  }else{
    return(paste0(2))
  }
}

apoeGenotyping <- function(rs429358_genotype, rs7412_genotype){
  
  if(length(rs429358_genotype) == 0 | length(rs7412_genotype) == 0){
    return(NA)
  }
  
  if(rs429358_genotype == "T/T" & rs7412_genotype == "T/T"){
    return("e2/e2")
  }else if(rs429358_genotype == "T/T" & rs7412_genotype %in% c("C/T", "T/C")){
    return("e2/e3")
  }else if(rs429358_genotype %in% c("C/T", "T/C") & rs7412_genotype %in% c("C/T", "T/C")){
    return("e2/e4")
  }else if(rs429358_genotype == "T/T" & rs7412_genotype == "C/C"){
    return("e3/e3")
  }else if(rs429358_genotype %in% c("C/T", "T/C") & rs7412_genotype == "C/C"){
    return("e3/e4")
  }else if(rs429358_genotype == "C/C" & rs7412_genotype == "C/C"){
    return("e4/e4")
  }else{
    return(NA)
  }
}