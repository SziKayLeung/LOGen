#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Purpose:
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
##
## ---------- Functions -----------------
## 
## final_num_iso
## numIso_relationship
##   
## ---------- Notes -----------------
## 

source(paste0(LOGEN_ROOT, "aesthetics_basics_plots/pthemes.R"))
suppressMessages(library(scales))


## ------------------- total_num_iso

# Aim: plot the total number of isoforms per target gene in finalised dataset by structural category 
# Input: 
  # class.files = df: classification file 
  # input_title = str: title of the plot
  # dataset = str:"dataset" then split by dataset (if column exists in clasification file)
  # glimit = int: number of genes to display (top X ranked)

total_num_iso <- function(class.files, input_title, dataset = NA, glimit = NA){
  
  # structural category colours 
  cate_cols <- c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),alpha("#F8766D",0.8),alpha("#F8766D",0.3),alpha("#808080",0.3))
  
  # number of isoforms per gene
  nIso <- class.files %>% group_by(associated_gene) %>% tally %>% arrange(-n) %>% dplyr::rename("totaln" = "n")
  
  if(!is.na(dataset)){
    
      p <- class.files %>% group_by(associated_gene, dataset) %>% tally %>%
        ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = dataset)) +
        scale_fill_discrete(name = "Dataset")
    
  }else{
    
    if(!is.na(glimit)){
      
      p <- class.files %>% group_by(associated_gene, structural_category) %>% tally %>% 
        full_join(., nIso, by = "associated_gene") %>%
        filter(associated_gene %in% nIso$associated_gene[1:glimit]) %>%
        ggplot(.,aes(x = reorder(associated_gene, -totaln), y = n, fill = structural_category)) +
        scale_fill_manual(name = "Isoform classification", values = cate_cols)
      
    }else{
      
      p <- class.files %>% group_by(associated_gene, structural_category) %>% tally %>%
        ggplot(.,aes(x = reorder(associated_gene,-n), y = n, fill = structural_category)) +
        scale_fill_manual(name = "Isoform classification", values = cate_cols)
      
    }
    
  }
  
  p <- p + geom_bar(stat = "identity") + mytheme + labs(x = "Top ranked genes", y = "Number of isoforms", title = input_title) + mytheme +
    theme(legend.position = "top") +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10),labels = label_comma())
  
  return(p)
}


## ---------- numIsoGene -----------------

numIsoGene <- function(class.files, stats=FALSE){
  
  isoPerGene <- SQANTI_gene_preparation(class.files) %>% 
    mutate(isoformType = ifelse(structural_category %in% c("FSM","ISM"),"Known","Novel")) 
  
  sum(isoPerGene$nIso)
  length(unique(isoPerGene$associatedGene))
  nrow(isoPerGene)
  
  isoPerGeneTally <- isoPerGene %>% group_by(isoformType, nIsoCat) %>% dplyr::count(nIso) 
  #isoPerGeneTally <- isoPerGeneTally %>% mutate(isoformType = factor(isoformType, levels = c("Novel","Known")))
  sum(isoPerGeneTally$n)
  if(isTRUE(stats)){
    return(isoPerGeneTally)
  }
  
  p <- isoPerGeneTally %>%
    mutate(totalNumGenes = length(unique(isoPerGene$associatedGene))) %>% 
    mutate(Perc = n/totalNumGenes * 100) %>% 
    ggplot(., aes(x=nIsoCat, y = Perc, fill = isoformType)) +
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_manual(values =  c("#00BFC4","#F8766D")) +
    labs(x ="Number of isoforms", y = "Genes (%)", fill = "Isoform classification", title = "\n") +
    mytheme + 
    theme(legend.position = c(0.8,0.8))
  
  return(p)
}


## ------------------- length_exon_description

# Aim: generate two statements about mean, sd and range of isoform length and number of exons
# for paper reporting
# Input:  
  # class.files = df: SQANTI classification file

length_exon_description <- function(class.files){
  print(paste0("Length: mean = ", signif(mean(class.files$length),3),
               " s.d = ", signif(sd(class.files$length),3),
               " range: ", signif(min(class.files$length),3),"-", signif(max(class.files$length),3)))
  
  print(paste0("Exon number: median = ", signif(median(class.files$exons),3),
               " s.d = ", signif(sd(class.files$exons),3),
               " range: ", signif(min(class.files$exons),3),"-", signif(max(class.files$exons),3)))
}


## ------------------- tabulate/plot_structural_cate

# Aim: tabulate the total number of isoforms in finalised dataset by structural category 
# Input:  
  # class.files = df: SQANTI classification file

tabulate_structural_cate <- function(class.files){
  
  # group by structural_category and tally with proportion
  output <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100)
  
  return(output)
}


# Aim: plot the total number of isoforms in finalised dataset by structural category 
# Input:  
  # class.files = df: SQANTI classification file
# Output:
  # p = bar-plot of the number of isoforms by structural category (rotated for paper purposes)

plot_structural_cate <- function(class.files, column="structural_category",rotate=FALSE){
  
  # colours
  structuralcolours = data.frame(
    "FSM" = rgb(249,176,161,maxColorValue = 255),
    "ISM" = rgb(228,196,116,maxColorValue = 255),
    "NIC" = rgb(204,220,151,maxColorValue = 255),
    "NNC" = rgb(108,220,164,maxColorValue = 255),
    "Fusion" = rgb(114, 215, 217,maxColorValue = 255),
    "Genic_Genomic" = rgb(226, 226, 243,maxColorValue = 255),
    "Genic_Intron" = rgb(226, 226, 243,maxColorValue = 255),
    "Antisense" = rgb(108,204,251,maxColorValue = 255),
    "Intergenic"= rgb(252,140,220,maxColorValue = 255)
  ) %>% gather(.,category,colours,factor_key=TRUE)
  
  cate_cols <- c(alpha("#00BFC4",0.8),alpha("#00BFC4",0.3),alpha("#F8766D",0.8),alpha("#F8766D",0.3))
  
  
  # keep only structural category that are detected
  dat <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100) %>% 
    mutate(structural_category = factor(structural_category, levels = structuralcolours$category)) %>%
    filter(perc > 0) %>% filter(!!sym(column) != "NA")
  
  p <- ggplot(dat, aes(x = !!sym(column), y = perc, fill = structural_category)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "Structural category", y = "Percentage of isoforms (%)", title = NULL) +
    #geom_text(aes(label=n), vjust=-1) + 
    theme(legend.position = "None")
  
  if(isTRUE(rotate)){
    p <- p + coord_flip() + scale_x_discrete(limits = rev(dat[["structural_category"]])) 
  }
  
  if(length(dat[[column]]) > 5){
    p <- p + scale_fill_manual(values = c(cate_cols, wes_palette("GrandBudapest1")[3], alpha("#808080",0.3),
                                          wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[4])) 
  }else{
    p <- p + scale_fill_manual(name = "", values = c(cate_cols,alpha("#808080",0.3))) 
  }
  

  
  return(p)
}

## ------------------- numIso_relationship

# Aim: plot the relationship between the number of isoforms and gene length, exon number and expression
# Pre-requisites:
  # developed for targeted dataset only whereby generated refGencode output from FICLE
# Input:  
  # classfiles = df: SQANTI classification file
  # geneExp = df: gene normalized expression with columm <isoform><sample><normalised_counts> (nb: isoform (XXX) = geneID from pb.XXX.YYY)
  # transExp = df: transcript normalised expression with column <sample><isoform><normalised_counts><associated_gene>
  # refGencode = df: FICLE output with column <associated_gene><MaxGeneLength><Maxexons>
# Output:
  # 4 scatter plots of the features vs number of isoforms
  # statistics from correlation

numIso_relationship <- function(classfiles, geneExp, transExp, refGencode){
  
  # extract geneID from classfiles
  classfiles <- classfiles %>% mutate(gene = word(isoform,c(2),sep=fixed(".")))
  
  # number of isoforms
  TargetIsoNum <- classfiles %>% group_by(associated_gene) %>% tally(name="numIso") %>% arrange(-numIso)
  
  # median gene expression across al the samples
  TargetGeneExp <- aggregate(normalised_counts ~ isoform, geneExp, median)
  
  # merging all files into one dataframe
  dat <- merge(classfiles[,c("associated_gene","gene")], TargetGeneExp, by.x = "gene", by.y = "isoform") %>% distinct()
  dat <- merge(refGencode, dat, by = "associated_gene")
  dat <- merge(TargetIsoNum, dat, by = "associated_gene")
  
  # sum normalized expression of isoforms across all samples
  TargetIsoExp <- aggregate(normalised_counts ~ isoform, transExp, sum)
  TargetIsoExp <- TargetIsoExp %>% 
    merge(., class.files$targ_filtered[,c("isoform","associated_gene")], all.y = T) %>% 
    arrange(-normalised_counts) %>% 
    mutate(associated_gene = factor(associated_gene, levels = TargetIsoNum$associated_gene))
    
  
  # plots
  # p1: scatter plot of the number of isoforms against gene length
  # p2: scatter plot of the number of isoforms against number of exons
  # p3: scatter plot of the number of isoforms against normalized gene expression
  # p4: geom tile of the ranked isoform by abundance (sum across all samples) 
  p1 <- ggplot(dat, aes(x = MaxGeneLength, y = numIso)) + geom_point() +
    labs(x = "Gene length (kb)", y = "Number of isoforms") + mytheme + 
    geom_text_repel(aes(label = associated_gene), vjust = -0.5, max.overlaps = 20) +
    scale_x_continuous(labels = ks)
  
  p2 <- ggplot(dat, aes(x = Maxexons, y = numIso)) + geom_point() + 
    geom_text_repel(aes(label = associated_gene), vjust = -0.5, max.overlaps = 20) +
    labs(x = "Number of exons", y = "Number of isoforms") + mytheme
  
  p3 <- ggplot(dat, aes(x = log10(normalised_counts), y = numIso)) + geom_point() + 
    geom_text_repel(aes(label = associated_gene), vjust = -0.5, max.overlaps = 20) +
    labs(x = "Normalized counts (log10)", y = "Number of isoforms") + mytheme
  
  p4 <- ggplot(TargetIsoExp, aes(x = as.factor(associated_gene), 
                                 y = reorder(isoform,normalised_counts), 
                                 fill = log10(as.numeric(normalised_counts)))) +
    geom_tile() +
    theme(legend.position = "top", axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_fill_viridis(option="viridis", name = "normalized_counts (log10)") + 
    labs(x = "Gene", y = "Isoform (ranked by total abundance)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # correlation
  message("correlation: gene length vs number of isoforms")
  res <- cor.test(dat$MaxGeneLength, dat$numIso)
  print(res)
  print(res$p.value)
  
  message("correlation: number of exons vs number of isoforms")
  res <- cor.test(dat$Maxexons, dat$numIso)
  print(res)
  print(res$p.value)
  
  message("correlation: gene expression vs number of isoforms")
  res <- cor.test(dat$normalised_counts, dat$numIso)
  print(res)
  print(res$p.value)
  
  return(list(p1,p2,p3,p4))
}

