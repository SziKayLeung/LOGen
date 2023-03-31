
library("wesanderson")
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "differential_analysis/base_DIU.R"))

## ------------------- CalculateIFMean

# Aim: run differntial isoform usage using functions() borrowed from tappAS
# 1. for each isoform, calculate the mean across the group (i.e. WT and TG)
# 2. for each isoform, calculate the isoform fraction by dividing the isoform mean over the sum of all isoform mean
# 3. classify minor isoform if isoform fraction across all groups < 0.5
# Input:  
  # transMatrixRaw = raw expression matrix (col.names = associated_gene;<sample1><sample2>)(raw.names = isoformID)
  # transMatrix = normalised expression matrix (col.names = associated_gene;<sample1><sample2>)(raw.names = isoformID)
  # classf = SQANTI classification file 
  # myfactors = df: row.names = sampleNames as in col.names(transMatrixRaw); col.names = <Replicate>
    # can have other columns, but must have the "Replicate colum"
    # Replicate of groups, 1,1,1,2,2,2 (ordered)
  # filteringType = <PROP,FOLD>
  # filterFC: note Foldchange for FOLD filteringType or proportion for PROP (default: 2 for Foldchange)
# Output:
  # DIU output (p.value and FDR from edgeR output)

runDIU <- function(transMatrixRaw,transMatrix,classf,myfactors,filteringType,filterFC){
  
  # create dataframe: row.names = isoform, col.names = "id" (gene_name)
  genesdf <- classf %>% select(associated_gene) %>% `colnames<-`(c("id"))
  
  # remove row.names and associated_gene column from raw expression value (result as matrix) 
  # while recording the associated_gene names for each row
  # while recording the transcriptID for each row
  isoformID <- rownames(transMatrixRaw)
  genes <- transMatrixRaw$associated_gene
  rownames(transMatrixRaw) <- NULL
  transMatrixRaw <- transMatrixRaw %>% dplyr::select(-associated_gene)
  
  # filter isoforms 
  cat("\nRead ", nrow(transMatrixRaw), " raw counts transcripts expression data rows")
  cat(paste0("\nFiltering new transcript matrix by ",filteringType,"...\n"))
  filterFC <- as.numeric(filterFC)
  # keeptrans = transcripts retained after filtering, represented as transMatrixRaw[index] 
  keeptrans <- minorFoldfilterTappas(transMatrixRaw, genes, filterFC, minorMethod=filteringType)
  transMatrixRaw <- transMatrixRaw[keeptrans,]
  
  # replace the index with original isoformID stored
  row.names(transMatrixRaw) = isoformID[as.numeric(rownames(transMatrixRaw))]
  
  # differential isoform usage
  cat("\nUsing EdgeR\n")
  result = spliceVariant.DS(transMatrixRaw, genesdf, myfactors)
  result$adj_PVALUE = p.adjust(result$PValue, method = "fdr")
  
  # calculate podium change
  transMatrix = transMatrix[rownames(transMatrix) %in% classf$isoform, ]
  pcList = podiumChange(transMatrix, genesdf, myfactors)
  pcdf <- as.data.frame(pcList$podiumChange)
  pcdf["totalChange"] <- pcList$totalChange
  result_diu <- merge(result, pcdf, by=0)
  colnames(result_diu) <- c("Gene","p.value","FDR","podiumChange","totalChange")
  
  # order by p.value
  result_diu <- result_diu %>% arrange(p.value)
  return(result_diu)
}


## ------------------- CalculateIFMean

# Aim: calculate usage of the isoforms with calculated mean across the groups (WT and TG)
# 1. for each isoform, calculate the mean across the group (i.e. WT and TG)
# 2. for each isoform, calculate the isoform fraction by dividing the isoform mean over the sum of all isoform mean
# 3. classify minor isoform if isoform fraction across all groups < 0.5
# Input:  
  # rawExp = normalised expression of only the isoforms of interest (already subsetted by gene)
    # rownames = isoformID (i.e. PB.1.1)
    # colnames = samples (i.e ONT_S19)
  # pheno = phenotype file 
    # sample column: same as rawExp colnames (i.e. ONT_S19)
    # col    column: sample + phenotype for downstream subsetting (ONT_S19_WT)

CalculateIFMean <- function(rawExp, pheno){
  
  # replace column names of rawExp with mathing col (with group) 
  names(rawExp) <- pheno$col[match(names(rawExp), pheno$sample)]
  
  # calculate mean expression of case vs Control 
  # meanisoexp:
    # row.names = isoformID
    # column 1 = group1_mean 
    # column 2 = group2_mean
  group1 = as.character(unique(pheno$group)[1])
  group2 = as.character(unique(pheno$group)[2])
  cat("Group 1:", group1, "\n")
  cat("Group 2:", group2, "\n")
  meanisoexp = cbind(rawExp %>% dplyr::select(contains(group1)) %>% apply(.,1,mean) %>% reshape2::melt(), 
                     rawExp %>% dplyr::select(contains(group2)) %>% apply(.,1,mean) %>% reshape2::melt()) %>% 
    `colnames<-`(c(paste0(group1,"_mean"), paste0(group2, "_mean")))   
  
  # determine the isoform fraction by mean/sum(mean)
  # same format as meanisoexp (only calculated percentage)
  IF = as.data.frame(apply(meanisoexp, 2, function(x) x/sum(x) * 100))
  
  # lowly abundant transcripts = mean expression of group 1 and group 2 < 5
  lowly_abundant = IF[IF[1] < 5 & IF[2] < 5,]
  if(nrow(lowly_abundant) > 0){
    minor_proportions <- data.frame(sum(lowly_abundant[1]),sum(lowly_abundant[2]))
    colnames(minor_proportions) <- colnames(IF)
    # remove abundant isoform proportion
    IF = IF[-which(rownames(IF) %in% rownames(lowly_abundant)),]  
    
    # include the minor proportion sum 
    IF <- rbind(IF, minor_proportions)
    rownames(IF)[nrow(IF)] <- "Minor"
  }

  IF <- IF %>% tibble::rownames_to_column(., var = "isoform")
  
  return(IF)
}


## ------------------- CalculateIFSample

# Aim: calculate usage of the isoforms per sample
# 1. for each isoform and each sample, calculate the fraction (FL/sum(FL) for that sample)
# 2. classify minor isoform in each sample if proporition < 5%
# Input:  
  # rawExp = normalised expression of only the isoforms of interest (already subsetted by gene)
  # rownames = isoformID (i.e. PB.1.1)
  # colnames = samples (i.e ONT_S19)
# Output:
  # df: colnames = sample, isoform, perc,.....,Iso
    # ... = phenotype file columns
    # Iso = <Minor, isoformID> = isoform classification 

CalculateIFSample <- function(rawExp, pheno){
  
  # calculate proportion for each column (i.e each sample)
  IF <- apply(rawExp, 2, function(x) x/sum(x) * 100) %>% 
    reshape2::melt() %>% 
    `colnames<-`(c("isoform", "sample", "perc")) %>% 
    merge(.,pheno,by="sample") %>%
    mutate(Iso = ifelse(perc < 5,"Minor",as.character(isoform)))
  
  return(IF)
}


## ------------------- CalculateIFSample_TimeSeries

# Aim: after calculating usage of the isoforms per sample, dissect for time_series
# 1. for each isoform, sum the percentage across group and time
# 2. filter results by just keeping the top ranked (most abundant) isoform
#    most abundant determined by greatest summed abundace across all the samples
# Input:  
  # CalculateIFSampleOutput = output from CalculateIFSample()
    # ensure has columns: sample,group,time,iso

CalculateIFSample_TimeSeries <- function(IFSampleOutput,rank=3){
  
  IF <- IFSampleOutput %>% group_by(sample,group,time,Iso) %>% tally(perc) %>% mutate(GroupIso = paste0(group, Iso))
  cat("Keeping only the", rank, "most abundant isoforms\n")
  major_isoforms <- as.data.frame(IF %>% group_by(Iso) %>% tally(n)) %>% arrange(-n) %>% .[1:rank,"Iso"]
  IF <- IF %>% filter(Iso %in% major_isoforms)
  IF$Iso <- factor(IF$Iso, levels = major_isoforms)
  
  return(IF)
}


## ------------------- plotIF

plotIF <- function(gene,Exp,pheno,cfiles,design="case_control"){
  
  print(gene)

  # subset expression by all the isoforms associated with the gene 
  iso = subset(cfiles, associated_gene == gene) 
  isoexp = Exp[which(rownames(Exp) %in% iso$isoform),] %>% dplyr::select(-c("associated_gene"))

  if(nrow(isoexp) > 1){
    # groups
    group1 <- levels(pheno$group)[[1]]
    group2 <- levels(pheno$group)[[2]]
    
    ##--- Method 1: determine isoform fraction by mean expression of isoform over sum of mean expression of all isoforms
    IF1 <- CalculateIFMean(isoexp,pheno)
    
    ##--- Method 2: determine isoform fraction for each sample
    IF2 <- CalculateIFSample(isoexp,pheno)
    
    p1 = IF1 %>% reshape2::melt() %>% `colnames<-`(c("Var1", "Var2","value")) %>% 
      mutate(group = factor(word(Var2,c(1),sep = fixed("_")), levels = c(group1,group2))) %>%
      ggplot(., aes(x = reorder(Var1,-value), y = value, fill = group)) + geom_bar(stat = "identity", position = position_dodge()) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = c(0.9,0.9)) +
      labs(x = "Isoform", y = "Isoform Fraction (%)", title = gene) + mythemeNoLegend + 
      scale_fill_manual(values = c(label_colour(group1),
                                   label_colour(group2)), " ",
                        labels = c(label_group(group1),label_group(group2))) 
    
    p2 <- IF2 %>% group_by(sample,group,time,Iso) %>% tally(perc) %>%
      mutate(group = factor(group, levels = c(group1,group2))) %>%
      ggplot(., aes(x = reorder(Iso,-n), y = n, fill = group)) + geom_boxplot() +
      labs(x = "Isoform", y = "Isoform Fraction (%)", title = gene) + 
      scale_fill_manual(values = c(label_colour(group1),label_colour(group2)), " ",
                        labels = c(label_group(group1),label_group(group2)))  + mythemeNoLegend + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = c(0.9,0.9))
    
    if(design == "time_series"){
      IF3 <- CalculateIFSample_TimeSeries(IF2,rank=3)
      IF3$group <- factor(IF3$group, levels = c(group1,group2))
      p3 <- ggplot(IF3, aes(x = as.factor(time), y = n, shape = Iso, colour = group)) + geom_point(size = 3)  +
        stat_summary(data=IF3, aes(x=as.factor(time), y=n, group=GroupIso), fun="mean", geom="line", linetype = "dotted") +
        labs(x = "Age (months)", y = "Isoform Fraction (%)", title = gene) +
        scale_colour_manual(values = c(label_colour(group1),label_colour(group2)), " ",
                            labels = c(label_group(group1),label_group(group2)),
                            guide="none") + mythemeNoLegend + 
        theme(legend.position = "bottom",legend.direction = "vertical", legend.box="horizontal", legend.margin=margin()) + 
        scale_shape_manual(values=c(19, 0, 15),name = "Isoform", guide = guide_legend(title.position = "left"))
      
      output = list(p1,p2,p3)
    }else{
      output = list(p1,p2)
    }
    
  }else{
    p1 = ggplot() + theme_void()
    p2 = ggplot() + theme_void()
    output = list(p1,p2)
    if(design == "time_series"){
      p3 = ggplot() + theme_void()
      output = list(p1,p2,p3)
    }
  }
  
  return(output)
}
