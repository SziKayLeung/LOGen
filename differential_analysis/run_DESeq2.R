suppressMessages(library("DESeq2"))
suppressMessages(library("ggrepel"))


run_DESeq2 <- function(expression, phenotype, exprowname, threshold=10, controlname="Control", interaction="Off",test){
  
  # input phenoptype characteristics as factor
  phenotype$time <- as.factor(phenotype$time)
  phenotype$group <- as.factor(phenotype$group)
  phenotype$group <- relevel(phenotype$group,controlname)
  str(phenotype$group)
  
  # common columns between phenotype and expression
  col_match <- intersect(phenotype$sample, colnames(expression))
  cat("Number of samples:", length(col_match),"\n")
  
  # ensure expression column and phenotype rows are in the same order
  rownames(expression) <- expression[[exprowname]]
  phenotype <- phenotype %>% filter(sample %in% col_match) 
  expression <- expression %>% dplyr::select(phenotype$sample)

  if(nrow(phenotype) != ncol(expression)){
    print("Mismatched number of samples in phenotype and expression file")
  }
  
  cat("Expression file:\n")
  print(head(expression))
  # create DESeq2 data object
  # normalises reads within DESeq2
  # design is linear model – last expression is tested for differential expression
  
  if(interaction == "Off"){
    cat("Modelling interaction effect: off")
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), 
                                  colData = phenotype, 
                                  design = ~ group + time)

  }else{
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(expression), 
                                  colData = phenotype, 
                                  design = ~ group + time + group:time)
    
  }

  
  
  # estimate size factors to account for differences in sequencing depth
  # if all samples have exact same sequencing depth, size factor should be ~ 1
  dds <- estimateSizeFactors(dds)
  
  # pre-filtering the data set
  # remove the rows of the DESeqDataSet that have no counts, or only a single count across all samples
  # arbitrary threshold = 10 (default)
  # minimum 2 samples per group, and expect minimum 2 FL reads per isoform
  cat("No of isoforms before filtering:", nrow(dds),"\n")
  cat("Filtering isoform on count threshold:", threshold,"\n")
  dds <- dds[rowSums(counts(dds)) > threshold, ] 
  cat("No of isoforms after filtering:", nrow(dds),"\n")
  
  #p <- cbind(reshape2::melt(sizeFactors(dds)), reshape2::melt(colSums(counts(dds)))) %>% 
  #  `colnames<-`(c("sizefactors", "nreads")) %>% 
  #  tibble::rownames_to_column("sample") %>% 
  #  mutate(sample = str_remove(sample,"ont_")) %>%
  #  ggplot(., aes(x = sizefactors, y = nreads)) + geom_point() +
  #  geom_label_repel(aes(label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  #  theme_bw() + labs(y = "Number of reads", x = "Size Factors")
  

  # normalization to stabilize variance (regularized logarithm)
  #rld <- rlog(dds, blind = FALSE)
  
  # PCA plot
  #pcaData <- plotPCA(rld, intgroup = c("time", "group"), returnData = TRUE)
  #pcaData$sample <- sapply(pcaData$name, function(x) str_remove(x, "ont_"))
  #percentVar <- round(100 * attr(pcaData, "percentVar"))
  #p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape = group.1)) +
  #  geom_point(size =3) +
  #  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance")) +
  #  coord_fixed() +
  #  geom_label_repel(aes(label = sample), segment.color = 'grey50') 
  
  
  # normalised counts for downstream plotting
  norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))
  
  # run differential analysis
  if(test == "Wald"){
    cat("Running Wald test\n")
    dds_Wald <- DESeq(dds, test="Wald")
    res_Wald <- as.data.frame(results(dds_Wald)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj)
    output <- list(dds_Wald, res_Wald, norm)
    names(output) <- c("dds_Wald", "res_Wald","norm_counts")
  }else if(test == "LRT"){
    cat("Running LRT test\n")
    dds_LRT <- DESeq(dds, reduced=~group + time, test="LRT")
    res_LRT <- as.data.frame(results(dds_LRT),stringsAsFactors=FALSE) %>% arrange(padj)
    output <- list(dds_LRT, res_LRT, norm)
    names(output) <- c("dds_LRT", "res_LRT","norm_counts")
  }else{
    print("test either <Wald/LRT>")
  }

  return(output)
}


# plots
# only keep the isoforms that are in the sqanti classification file and not in the expression file
anno_DESeq2 <- function(deseq_output,class.files,phenotype,controlname="Control"){
  
  phenotype$group <- as.factor(phenotype$group)
  phenotype$group <- relevel(phenotype$group,controlname)
  
  deseq_output$anno_res <-  merge(deseq_output[[2]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
  
  # normalised counts
  deseq_output[[3]] <- merge(deseq_output[[3]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
  deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)

  return(deseq_output)
}
