suppressMessages(library("DESeq2"))
suppressMessages(library("ggrepel"))
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html


# changing the threshold (pre-filtering) impacts the adjusted p-value as increased number of multiple testing
# filter before normalising

run_DESeq2 <- function(test,expression,phenotype,exprowname=NULL,threshold=10,controlname="Control",design="time_series",interaction="Off",groupvar="factor"){
  
  
  # input phenoptype characteristics as factor
  if(groupvar=="factor"){
    cat("Group variable as factor\n")
    phenotype$group <- as.factor(phenotype$group)
    phenotype$group <- relevel(phenotype$group,controlname)
  }else{
    cat("Group variable as continuous variable\n")
    phenotype$group <- as.numeric(phenotype$group)
  }

  if(design == "time_series"){
    phenotype$time <- as.factor(phenotype$time)
  }
  
  # common columns between phenotype and expression
  col_match <- intersect(phenotype$sample, colnames(expression))
  cat("Number of samples:", length(col_match),"\n")
  
  # ensure expression column and phenotype rows are in the same order
  if(!is.null(exprowname)){
    #rownames(expression) <- expression[[exprowname]]
    expression <- expression %>% tibble::column_to_rownames(exprowname)
  }
  phenotype <- phenotype %>% filter(sample %in% col_match) 
  expression <- expression %>% dplyr::select(phenotype$sample)
  rownames(phenotype) <- phenotype$sample
  phenotype <- phenotype %>% dplyr::select(-sample)
  if(all(colnames(expression) == rownames(phenotype))==FALSE){
    print("ERROR: rownames and colnames in expression matrix and phenotype are not in the same order")
  }
  
  if(nrow(phenotype) != ncol(expression)){
    print("Mismatched number of samples in phenotype and expression file")
  }
  
  cat("Expression file:\n")
  print(head(expression))
  # create DESeq2 data object
  # normalises reads within DESeq2
  # counts need to be rounded to the nearest integer
  # design is linear model – last expression is tested for differential expression
  
  
  if(design == "time_series"){
    cat("Design: time_series\n")
    if(interaction == "Off"){
      cat("Modelling interaction effect: off\n")
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                    colData = phenotype,
                                    design = ~ group + time)
      
    }else{
      cat("Modelling interaction effect: on\n")
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                    colData = phenotype,
                                    design = ~ group + time + group:time)
      
    }
  }else{
    cat("Design: case_control\n")
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(expression)),
                                  colData = phenotype,
                                  design = ~ group)
  }
  
  
  
  # estimate size factors to account for differences in sequencing depth
  # if all samples have exact same sequencing depth, size factor should be ~ 1
  dds <- estimateSizeFactors(dds)
  
  # pre-filtering the data set
  # remove the rows of the DESeqDataSet that have no counts, or only a single count across all samples
  # arbitrary threshold = 10 (default)
  # minimum 2 samples per group, and expect minimum 2 FL reads per isoform
  cat("No of isoforms before filtering:", nrow(dds),"\n")
  # normalised counts for downstream plotting
  normAll <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))
  cat("Filtering isoform on count threshold (minimum):", threshold,"\n")
  dds <- dds[rowSums(counts(dds)) >= threshold, ]
  # normalised counts for downstream plotting
  norm <- counts(dds, normalized=TRUE) %>% reshape2::melt() %>% `colnames<-`(c("isoform", "sample", "normalised_counts"))
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
  
  # run differential analysis
  if(test == "Wald"){
    cat("Running Wald test\n")
    dds_output <- DESeq(dds, test="Wald")
    output_names <- c("dds_Wald", "res_Wald","norm_counts", "stats_Wald","norm_counts_all")
  }else if(test == "LRT"){
    cat("Running LRT test\n")
    # reduced model removes the interaction term 
    # significant DE genes will represent those genes that have differences in the effect of genotype over time
    dds_output <- DESeq(dds, reduced=~group+time, test="LRT")
    output_names <- c("dds_LRT", "res_LRT","norm_counts", "stats_LRT", "norm_counts_all")
  }else{
    print("test either <Wald/LRT>")
  }
  
  res <- as.data.frame(results(dds_output)) %>% tibble::rownames_to_column("isoform") %>% arrange(padj)
  stats <- as.data.frame(mcols(dds_output))
  output <- list(dds_output, res, norm, stats, normAll)
  names(output) <- output_names
  
  return(output)
}


# plots
# only keep the isoforms that are in the sqanti classification file and not in the expression file
anno_DESeq2 <- function(deseq_output,class.files,phenotype,controlname="Control",level,sig=0.05){
  
  phenotype$group <- as.factor(phenotype$group)
  phenotype$group <- relevel(phenotype$group,controlname)
  
  if(level == "transcript"){
    print("Processing via transcript")
    deseq_output$anno_res <-  merge(deseq_output[[2]],class.files[,c("isoform","associated_gene","associated_transcript","structural_category","subcategory")], by = "isoform", all.x = T)
    # normalised counts
    deseq_output[[3]] <- merge(deseq_output[[3]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
    deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)
    deseq_output[[5]] <- merge(deseq_output[[5]],class.files[,c("isoform","associated_gene","associated_transcript")], by = "isoform", all.x = T)
    deseq_output[[5]] <- merge(deseq_output[[5]],phenotype, by = "sample", all = T)
    
  }else{
    print("Processing via gene")
    class.files <- class.files %>% mutate(PB_associated_gene = word(isoform,c(2), sep = fixed(".")))
    gene_class_files <- unique(class.files[,c("associated_gene","PB_associated_gene")])
    deseq_output$anno_res <- merge(deseq_output[[2]], gene_class_files, by.x = "isoform", by.y = "PB_associated_gene")
    
    # normalised counts
    deseq_output[[3]] <- merge(deseq_output[[3]],gene_class_files, by.x = "isoform", by.y = "PB_associated_gene", all.x = T)
    deseq_output[[3]] <- merge(deseq_output[[3]],phenotype, by = "sample", all = T)
    deseq_output[[4]] <- merge(deseq_output[[4]], gene_class_files, by.x = 0, by.y = "PB_associated_gene")
  }
  
  deseq_output$anno_res <- deseq_output$anno_res %>% filter(padj < sig) %>% arrange(padj)
  
  return(deseq_output)
}


dissect_DESeq2 <- function(wald_res=wald_res, lrt_res=lrt_res){
  
  sig_iso <- list(
    genotype = setdiff(wald_res$isoform,lrt_res$isoform),
    progressive = lrt_res$isoform
  )
  
  for(i in 1:length(sig_iso)){
    print(paste("Number of significant hits in", names(sig_iso)[[i]], "effect :", length(sig_iso[[i]])))
  }
  
  anno_resTran_split <- list(
    genotype = wald_res %>% filter(isoform %in% sig_iso$genotype),
    progressive = lrt_res %>% filter(isoform %in% sig_iso$progressive)
  )
  
  names(anno_resTran_split) = c("genotype","progressive")
  return(anno_resTran_split)
}

