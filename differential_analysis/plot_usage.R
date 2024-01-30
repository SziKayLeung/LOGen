
library("wesanderson")
library("stringr")
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "differential_analysis/base_DIU.R"))
source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))

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
  # resultDIU: DIU output (p.value and FDR from edgeR output)
  # keptIso: expression of kept isoforms after filtering and fed into EdgeR 

runDIU <- function(transMatrixRaw,transMatrix,classf,myfactors,filteringType,filterFC){
  
  # create dataframe: row.names = isoform, col.names = "id" (gene_name)
  genesdf <- classf %>% select(associated_gene) %>% `colnames<-`(c("id"))
  
  # remove row.names and associated_gene column from raw expression value (result as matrix) 
  # while recording the associated_gene names for each row
  # while recording the transcriptID for each row
  isoformID <- rownames(transMatrixRaw)
  genes <- transMatrixRaw$associated_gene
  transMatrixRaw <- transMatrixRaw %>% dplyr::select(-associated_gene)
  
  # filter isoforms 
  cat("\nRead ", nrow(transMatrixRaw), " raw counts transcripts expression data rows")
  cat(paste0("\nFiltering new transcript matrix by ",filteringType,"...\n"))
  filterFC <- as.numeric(filterFC)
  # keeptrans = transcripts retained after filtering, represented as transMatrixRaw[index] 
  keeptrans <- minorFoldfilterTappas(data=transMatrixRaw,gen=genes, minorfilter=filterFC, minorMethod=filteringType)
  transMatrixRaw <- transMatrixRaw[rownames(transMatrixRaw) %in% keeptrans, ]  
  message("Number of kept isoforms: ", length(keeptrans))  
  #transMatrixRaw

  # differential isoform usage
  cat("\nUsing EdgeR\n")
  result = spliceVariant.DS(raw.counts=transMatrixRaw, feature_association=genesdf, factors=myfactors)
  if(!is.null(result)){
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
    output <- list(result_diu, transMatrixRaw)
    names(output) <- c("resultDIU","keptIso")
  }else{
    output <- NULL
  }
  return(output)
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
  rownames(meanisoexp) = rownames(rawExp)
  
  # determine the isoform fraction by mean/sum(mean)
  # same format as meanisoexp (only calculated percentage)
  IF = as.data.frame(apply(meanisoexp, 2, function(x) x/sum(x) * 100))
  
  if(length(unique(row.names(meanisoexp))) == 1){
    IF = IF %>% `colnames<-`(c("perc")) %>% tibble::rownames_to_column(var = "isoform") %>% spread(.,isoform,perc) 
    IF$isoform = unique(row.names(meanisoexp))
    IF = IF %>% select(isoform,contains("mean"))
    lowly_abundant <- data.frame()
  }else{
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
  }
  output <- list(IF, lowly_abundant)
 
  return(output)
}


## ------------------- CalculateIFSample

# Aim: calculate usage of the isoforms per sample
# 1. for each isoform and each sample, calculate the fraction (FL/sum(FL) for that sample)
# 2. classify minor isoform in each sample if proporition < 5% or specify the major isoforms (majorIso option)
# Input:  
  # rawExp = normalised expression of only the isoforms of interest (already subsetted by gene)
    # rownames = isoformID (i.e. PB.1.1)
    # colnames = samples (i.e ONT_S19)
  # pheno = phenotype file with column <sample>
  # majorIso = vector of isoforms (optional)
# Output:
  # df: colnames = sample, isoform, perc,.....,Iso
    # ... = phenotype file columns
    # Iso = <Minor, isoformID> = isoform classification 

CalculateIFSample <- function(rawExp, pheno, majorIso){
  
  # calculate proportion for each column (i.e each sample)
  IF <- apply(rawExp, 2, function(x) x/sum(x) * 100) %>% 
    reshape2::melt() %>% 
    `colnames<-`(c("isoform", "sample", "perc")) %>% 
    merge(.,pheno,by="sample")
  
  if(!is.null(majorIso)){
    cat("Minor isoforms specified\n")
    IF <- IF %>% mutate(Iso = ifelse(isoform %in% majorIso, as.character(isoform), "Minor"))
  }else{
    cat("Minor isoforms determined if less than 5%")
    IF <- IF %>% mutate(Iso = ifelse(perc < 5,"Minor",as.character(isoform)))
  }
  
  return(IF)
}


## ------------------- CalculateIFSample_TimeSeries

# Aim: after calculating usage of the isoforms per sample, dissect for time_series
# 1. for each isoform, sum the percentage across group and time
# 2. filter results by keeping specific isoforms or
#    filter results by just keeping the top ranked (most abundant) isoform (if rank != NULL)
#    most abundant determined by greatest summed abundace across all the samples
# Input:  
  # CalculateIFSampleOutput = output from CalculateIFSample()
    # ensure has columns: sample,group,time,iso

CalculateIFSample_TimeSeries <- function(IFSampleOutput,rank=NULL,isoSpecific=NULL){
  
  # iso = PB.XXX.XX whereas isoform = GeneY.XXX.XX
  IF <- IFSampleOutput %>% group_by(sample,group,time,isoform) %>% tally(perc) %>% mutate(GroupIso=paste0(group,isoform))
  
  if(!is.null(rank) & rank > 0){
    cat("Keeping only the", rank, "most abundant isoforms\n")
    major_isoforms <- as.data.frame(IF %>% group_by(isoform) %>% tally(n)) %>% arrange(-n) %>% .[1:rank,"isoform"]
    IF <- IF %>% filter(isoform %in% major_isoforms)
    IF$isoform <- factor(IF$isoform, levels = major_isoforms)
  }else{
    cat("Keeping specific isoforms")
    IF <- IF %>% filter(isoform %in% isoSpecific)
    IF$isoform <- factor(IF$isoform, levels = isoSpecific)
  }

  return(IF)
}


## ------------------- plotIF

# Aim: Plot the isoform fraction of specific gene
# Input:
  # ExpInput = df: matrix of normalised expression of all associated target reads 
    # use all normalised reads, including pre-filtered expression < 10 reads
    # row.names = isoformID, col.names = sampleID
  # pheno = df: phenotype file 
  # cfiles = df: class.files
  # design = str: <case_control><time_series>
  # majorIso =  vector: list of major isoforms specified in IF 
  # rank = numeric: ranking for CalculateIFSample_TimeSeries
  # isoSpecific = str: isoform ID for CalculateIFSample_TimeSeries

plotIF <- function(gene,ExpInput,pheno,cfiles,design="case_control",majorIso=NULL,rank=3,isoSpecific=NULL,stats=FALSE){
  
  print(gene)

  # subset expression by all the isoforms associated with the gene 
  iso = subset(cfiles, associated_gene == gene) 
  isoexp = ExpInput %>% tibble::rownames_to_column("isoform") %>% 
    filter(isoform %in% iso$isoform) %>% 
    tibble::column_to_rownames(., var = "isoform") #%>% 
    #dplyr::select(-c("associated_gene"))

  if(nrow(isoexp) > 1){
    # groups
    group1 <- levels(pheno$group)[[1]]
    group2 <- levels(pheno$group)[[2]]
    
    ##--- Method 1: determine isoform fraction by mean expression of isoform over sum of mean expression of all isoforms
    IF1 <- CalculateIFMean(rawExp=isoexp,pheno)
    
    ##--- Method 2: determine isoform fraction for each sample
    IF2 <- CalculateIFSample(isoexp,pheno,majorIso)
    if(isTRUE(grepl("PB",IF2$isoform[1]))){
      IF2$renamedIsoform <- paste0("LR.",gsub(paste0("PB.",word(IF2$isoform[1],c(2),sep=fixed("."))), gene, IF2$isoform))
    }else{
      IF2$renamedIsoform <- IF2$isoform
    }

    p1 = IF1[1] %>% reshape2::melt(id="isoform") %>% `colnames<-`(c("Var1", "Var2","value")) %>% 
      select(Var1,Var2,value) %>%
      mutate(group = factor(word(Var2,c(1),sep = fixed("_")), levels = c(group1,group2))) %>%
      ggplot(., aes(x = reorder(Var1,-value), y = value, fill = group)) + geom_bar(stat = "identity", position = position_dodge()) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = c(0.9,0.9)) +
      labs(x = "Isoform", y = "Isoform fraction (%)", title = gene) + mytheme + 
      scale_fill_manual(values = c(label_colour(group1),
                                   label_colour(group2)), " ",
                        labels = c(label_group(group1),label_group(group2))) 
    
    # group by "iso" category (note some would be minor)
    IsoTally <- IF2 %>% group_by(Iso) %>% tally(perc) %>% arrange(-n)
    
    if(rank == 0){
      message("Selecting:", isoSpecific)
      IsoTally <- IsoTally %>% filter(Iso %in% isoSpecific) %>% arrange(-n)

    }else{
      # note if the isoSpecific isoform is not within the top 4, then take only the top however minus the isoSpecific isoforms
      message("Selecting the top ",rank, " isoforms, including:", isoSpecific)
      topRankedTally <- IsoTally %>% filter(Iso != "Minor") %>% top_n(rank - length(isoSpecific))
      IsoTally <- IsoTally %>% filter(Iso %in% isoSpecific) %>% bind_rows(topRankedTally) %>% arrange(-n)
    }
    
    if(isTRUE(grepl("PB",IsoTally$Iso[1]))){
      IsoTally$renamedIsoform <- paste0("LR.", gsub(paste0("PB.",word(IsoTally$Iso[1],c(2),sep=fixed("."))), gene, IsoTally$Iso))
    }else{
      IsoTally$renamedIsoform <- IsoTally$Iso
    }

    p2 <- IF2 %>% filter(isoform %in% IsoTally$Iso) %>% 
      group_by(sample,group,time,renamedIsoform) %>% tally(perc) %>%
      mutate(group = factor(group, levels = c(group1,group2)))
    
    if(length(unique(p2$renamedIsoform)) == 1){
      p2 <- ggplot(p2, aes(x = group, y = n, colour = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(aes(fill = group), size = 2, position = position_jitterdodge()) +
        labs(x = "Genotype", y = "Isoform fraction (%)", title = paste0(p2$renamedIsoform))
      
    }else{
      p2 <- ggplot(p2, aes(x = reorder(renamedIsoform,-n), y = n, colour = group)) + geom_boxplot(outlier.shape = NA) +
        geom_point(aes(fill = group), size = 2, position = position_jitterdodge()) +
        labs(x = "Isoform", y = "Isoform fraction (%)", title = gene)   
    }
    p2 <- p2 + scale_colour_manual(values = c(label_colour(group1),label_colour(group2)), " ",
                          labels = c(label_group(group1),label_group(group2))) + mytheme +
      scale_fill_discrete(guide="none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = c(0.9,0.9))
    
    if(design == "time_series"){
      # remove minor isoforms
      #IF2major <- IF2 %>% filter(Iso != "Minor") 
      IF3 <- CalculateIFSample_TimeSeries(IF2,rank=rank,isoSpecific=unique(IsoTally$Iso))
      IF3$renamedIsoform <- paste0("LR.",gsub(paste0("PB.",word(IF3$isoform[1],c(2),sep=fixed("."))), gene, IF3$isoform))
      IF3$group <- factor(IF3$group, levels = c(group1,group2))
      IF3$isoform <- factor(IF3$isoform, levels = c(unique(IsoTally$Iso)))
      p3 <- ggplot(IF3, aes(x = as.factor(time), y = n, colour = group)) +  
        geom_point(aes(colour = group), size = 2, position = position_jitterdodge(dodge.width=0.2)) +
        stat_summary(data=IF3, aes(x=as.factor(time), y=n, group=GroupIso), fun="mean", geom="line", linetype = "dotted")  +
        labs(x = "Age (months)", y = "Isoform fraction (%)", title = gene) +
        scale_colour_manual(values = c(label_colour(group1),label_colour(group2)), " ",
                            labels = c(label_group(group1),label_group(group2))) + mytheme + 
        theme(legend.justification = c(0, 1), legend.position = "top", legend.margin=margin()) + 
        scale_shape_manual(values=c(19, 0, 15),name = "Isoform", guide = guide_legend(title.position = "left")) 
      
      if(length(unique(IF3$renamedIsoform)) == 1){
        p3 <- p3 + labs(title = paste0(unique(IF3$renamedIsoform)))
      }else{
        p3 <- p3 + facet_grid(~renamedIsoform) + 
          theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))
      }

      output = list(p2,p3)
    }else{
      output = list(p1,p2)
    }
    
  }else{
    output = list(p1,p2)
    if(design == "time_series"){
      p3 = ggplot() + theme_void()
      output = list(p2,p3)
    }
  }
  
  if(isFALSE(stats)){
    return(output) 
  }else{
    return(IF2)
  }
}


## ------------------- plotIFAll

# Aim: Plot the isoform fraction of across all the genes
  # x = target genes  
  # y = isoform fraction 
  # coloured by structural category
  # Isoform fraction calculated from taking the mean of the isoform expression across all the samples
# Input:
  # Exp = matrix of normalised expression of all associated target reads (do not use filtered normalised reads with expression > 10 reads)
    # row.names = isoformID, col.names = sampleID
  # classf = class.files
  # pheno = phenotype file 
  # majorIso = list of isoforms that were retained from tappAS filtering and not removed

plotIFAll <- function(Exp,classf,pheno,majorIso){
  
  # Calculate the mean of normalised expression across all the samples per isoform
  meandf <- data.frame(meanvalues = apply(Exp,1,mean)) %>%
    rownames_to_column("isoform") %>% 
    # annotate isoforms with associated_gene and structural category
    left_join(., classf[,c("isoform","associated_gene","structural_category")], by = "isoform") %>% 
    # classify if isoform is "major" or "minor" from the majorIso argument
    # if the isoform is in the list, keep label, else label as "minor"
    mutate(majorminor = ifelse(isoform %in% majorIso,isoform,"minor")) 
  
  # Group meandf by associated_gene and calculate the sum of mean values for each group
  grouped <- aggregate(meandf$meanvalues, by=list(associated_gene=meandf$associated_gene), FUN=sum)
  
  # Calculate the proportion by merging back, and divide the meanvalues by the grouped values (x)
  merged <- meandf %>% 
    left_join(grouped, by = "associated_gene") %>%
    mutate(perc = meanvalues / x * 100) 
  
  # number of isoforms less than 1% across all target genes 
  minorLessThan1 <- merged %>% filter(perc < 1) %>% group_by(associated_gene) %>% tally()
  message("Median number of isoforms less than 1%")
  median(minorLessThan1$n)
  message("Range number of isoforms less than 1%")
  range(minorLessThan1$n)
  
  # Select major and minor isoforms
  major <- merged %>% filter(majorminor != "minor") %>% select(isoform, associated_gene, structural_category, perc)
  minor <- merged %>% filter(majorminor == "minor") 
  
  # Group the minor isoforms for each associated_gene and sum the percentages
  minorgrouped <- aggregate(minor$perc, by=list(associated_gene=minor$associated_gene), FUN=sum) %>%
    mutate(isoform = "minor", structural_category = "minor") %>% 
    dplyr::rename(perc = x) %>%
    select(isoform, associated_gene, structural_category,perc)
  
  # Tally the number of minor isoforms per associated_gene
  minortally <- minor %>% group_by(associated_gene) %>% tally() %>% mutate(isoform = "minor")
  
  # plot
  p <- rbind(major,minorgrouped) %>%
    full_join(., minortally, by = c("isoform", "associated_gene")) %>%
    ggplot(., aes(x = associated_gene, y = as.numeric(perc), fill = forcats::fct_rev(structural_category))) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    scale_color_manual(values = rep(NA, length(unique(minorgrouped$gene)))) + 
    mythemeNoLegend + labs(x = "Gene", y = "Isoform fraction (%)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(name = "Isoform Classification", values = rev(c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),
                                                                      alpha("#F8766D",0.8),alpha("#F8766D",0.3),
                                                                      alpha("#808080",0.3)))) +
    geom_text(aes(label=n),color="black",size=4,position=position_stack(vjust=0.5)) +
    theme(legend.position = "None")
  
  tab <- rbind(major,minorgrouped) %>%
    full_join(., minortally, by = c("isoform", "associated_gene"))
  
  
  return(list(p, tab))
}
