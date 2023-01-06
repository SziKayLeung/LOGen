# IF_plot 
# gene: input gene 
# gene_transcripts_input = tappas input file of list of geneName and transcript 
# normexp_input = tappas input file of normalised matrix 
# phenotype = phenotype file
IF_plot <- function(gene, gene_transcripts_input, normexp_input, df.phenotype, df.class.files){
  
  print(gene)
  genetrans = gene_transcripts_input
  normexp = normexp_input
  
  # keep only the isoforms that were retained in the classification file (final SQANTI curated)
  # some isoforms in the expression file retained from normalisation and low expression filtering
  # but filtered in sqanti due to intrapriming etc
  normexp = normexp[which(rownames(normexp) %in% df.class.files$isoform),]
  
  # subset expression by all the isoforms associated with the gene 
  iso = subset(genetrans, geneName == gene) 
  isoexp = normexp[which(rownames(normexp) %in% iso$transcript),]
  
  # relabel the expression file column with CASE and CONTROL from the phenotype file
  df.phenotype <- df.phenotype %>% mutate(col = paste0(sample,"_",group))
  names(isoexp) <- df.phenotype$col[match(names(isoexp), df.phenotype$sample)]
  
  # calculate mean expression of case vs Control 
  group1 = as.character(unique(df.phenotype$group)[1])
  group2 = as.character(unique(df.phenotype$group)[2])
  cat("Group 1:", group1, "\n")
  cat("Group 2:", group2, "\n")
  meanisoexp = cbind(isoexp %>% dplyr::select(contains(group1)) %>% apply(.,1,mean) %>% reshape2::melt(), 
                     isoexp %>% dplyr::select(contains(group2)) %>% apply(.,1,mean) %>% reshape2::melt()) %>% 
    `colnames<-`(c(paste0(group1,"_mean"), paste0(group2, "_mean")))   
  
  
  # Method 1: determine isoform fraction by mean expression of isoform over sum of mean expression of all isoforms
  IF1 = as.data.frame(apply(meanisoexp, 2, function(x) x/sum(x) * 100))
  
  # Method 2: determine isoform fraction for each sample
  IF2 <- as.data.frame(apply(isoexp, 2, function(x) x/sum(x) * 100))
  
  # lowly abundant transcripts 
  lowly_abundant = IF1[IF1[1] < 0.5 & IF1[2] < 0.5,]
  
  filter_lowly <- function(IF,type){
    
    # include minor proportions into graph
    if(type == "mean_all"){
      # sum the proportions of the minor isoforms
      minor_proportions <- data.frame(sum(lowly_abundant[1]),sum(lowly_abundant[2]))
      colnames(minor_proportions) <- colnames(IF)
    }else{
      # filter the expression to only include the lowly abundant isoforms
      minor_proportions <- data.frame(IF[which(rownames(IF) %in% rownames(lowly_abundant)),])
      # sum across each column (i.e the sample)
      minor_proportions <- apply(minor_proportions,2,sum)
    }
    
    # filter lowly expressed for plotting
    cat("\nNumber of isoforms pre-filtering for expression for plotting:", nrow(IF))
    IF = IF[-which(rownames(IF) %in% rownames(lowly_abundant)),]  
    
    # include the minor proporition sum into the graph
    IF <- rbind(IF, minor_proportions)
    rownames(IF)[nrow(IF)] <- "The minor isoforms"
    
    IF <- IF %>% rownames_to_column()
    
    cat("\nNumber of isoforms filtered plotting:", nrow(IF))
    
    return(IF)
  }
  
  if(nrow(lowly_abundant) != 0){
    IF1 <- filter_lowly(IF1,"mean_all")
    IF2 <- filter_lowly(IF2,"boxplot")
  }else{
    # for plotting
    IF1 <- IF1 %>% rownames_to_column()
    IF2 <- IF2 %>% rownames_to_column()
  }
  
  
  if(nrow(meanisoexp) != "1"){
    p1 = IF1 %>% reshape2::melt() %>% `colnames<-`(c("Var1", "Var2","value")) %>% 
      mutate(group = factor(word(Var2,c(1),sep = fixed("_")), levels = c(group1,group2))) %>%
      ggplot(., aes(x = Var1, y = value, fill = group)) + geom_bar(stat = "identity", position = position_dodge()) +
      labs(x = "Isoform", y = "Isoform Fraction (%)", title = gene) + mytheme + 
      scale_fill_manual(values = c(label_colour(group1),
                                   label_colour(group2)), " ",
                        labels = c(label_group(group1),label_group(group2))) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
    
    p2 <- reshape2::melt(IF2) %>% `colnames<-`(c("isoform", "sample", "proportion")) %>% 
      mutate(group = as.character(word(sample,c(-1),sep = fixed("_")))) %>%
      mutate(group = factor(group, levels = c(group1,group2)))  %>%
      ggplot(., aes(x = isoform, y = proportion, fill = group)) + 
      geom_boxplot() +
      labs(x = "Isoform", y = "Isoform Fraction (%)", title = gene) + 
      scale_fill_manual(values = c(label_colour(group1),label_colour(group2)), " ",
                        labels = c(label_group(group1),label_group(group2))) + mytheme +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
  }
  else{
    p1 = ggplot() + theme_void()
    p2 = ggplot() + theme_void()
  }
  
  output = list(p1,p2)
  return(p2)
}

#gene = "Trem2"
#genetrans =loaded$targ_ont$gene_transcripts
#normexp = loaded$targ_ont$input_normalized_matrix
#type = "ONT_Targeted"
#phenotype_input = phenotype$targ_ont
#p <- IF_plot_time_series("Trem2",loaded$targ_ont$gene_transcripts,loaded$targ_ont$input_normalized_matrix,phenotype$targ_ont,"ONT_Targeted")
IF_plot_time_series <- function(gene, gene_transcripts_input, normexp_input, phenotype_input, type){
  
  print(gene)
  genetrans = gene_transcripts_input
  normexp = normexp_input
  
  # subset expression by all the isoforms associated with the gene 
  iso = subset(genetrans, geneName == gene) 
  isoexp = normexp[which(rownames(normexp) %in% iso$transcript),]
  
  if(type == "Iso_Targeted"){  
    phenotype_input = phenotype_input %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                                 col = paste0(sample,"_",pheno_group,"_",time))
    names(isoexp) <- phenotype_input$col[match(names(isoexp), phenotype_input$sample)]
    
    # perform prefiltering by fold change 
    #normexp = DIU_time_analysis_filteronly(tappas_input_dir$iso,"3")
    
    
  }else if(type == "ONT_Targeted"){
    phenotype_input = phenotype_input %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                                 col = paste0(sample,"_",pheno_group,"_",time))
    names(isoexp) <- phenotype_input$col[match(names(isoexp), phenotype_input$sample)]
    
    # perform prefiltering by fold change 
    #normexp = DIU_time_analysis_filteronly(tappas_input_dir$ont,"3")
    
  }else{
    # whole transcriptome 
    # given fewer isoforms, no need to perform prefiltering of isoforms
    normexp = normexp_input
  }
  
  # calculate mean expression of WT vs Control 
  meanisoexp = cbind(isoexp %>% dplyr::select(contains("WT")) %>% apply(.,1,mean), 
                     isoexp %>% dplyr::select(contains("TG")) %>% apply(.,1,mean)) %>% `colnames<-`(c("WT_mean", "TG_mean"))   
  
  # determine isoform fraction by mean expression of isoform over sum of mean expression of all isoforms
  IF = apply(meanisoexp, 2, function(x) x/sum(x) * 100)
  
  if(type == "ONT_Targeted"){
    IF = as.data.frame(IF)
    IF = IF[IF$WT_mean > 0.5 & IF$TG_mean > 0.5,] %>% rownames_to_column() 
  }
  
  if(nrow(meanisoexp) != "1"){
    p = IF %>% reshape2::melt() %>% `colnames<-`(c("Var1", "Var2","value")) %>% 
      ggplot(., aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = position_dodge()) +
      labs(x = "Isoform", y = "Isoform Fraction (%)") + mytheme + 
      scale_fill_manual(values = c(label_colour("TG"),label_colour("WT")), "Genotype", labels = c("WT","TG")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "bottom")
    
    
    ### further subset by age 
    if(type == "isoseq"){ages = c("2","8")}else{ages = c("2","4","6","8")}
    meanisoexp_ages = list()
    for(a in 1:length(ages)){
      meanisoexp_ages[[a]] <- cbind(isoexp %>% dplyr::select(contains(paste0("WT_",ages[a]))) %>% apply(.,1,mean), 
                                    isoexp %>% dplyr::select(contains(paste0("TG_",ages[a]))) %>% apply(.,1,mean)) %>% `colnames<-`(c(paste0("WT_mean_",ages[a]),paste0("TG_mean_",ages[a])))
    }
    
    meanisoexp_ages = do.call(cbind, meanisoexp_ages)
    IF_ages = apply(meanisoexp_ages, 2, function(x) x/sum(x) * 100)
    
    IF_ages_plots = IF_ages %>% reshape2::melt() %>% mutate(age = factor(word(Var2, c(3), sep = fixed("_"))), group = factor(word(Var2, c(1), sep = fixed("_")),levels = c("WT","TG"))) 
    
    if(type %in% c("ONT_Targeted","Iso_Targeted")){
      IF_ages_plots = IF_ages_plots %>% mutate(Var1 = ifelse(value > 5, as.character(Var1), "Others"))
    }
    
    p1 = ggplot(IF_ages_plots, aes(x = group, y = value, fill = age)) + geom_bar(stat = "identity", position = position_dodge()) +
      facet_grid(~ Var1) + labs(x = "Genotype", y = "Isoform Fraction (%)") + mytheme +  theme(legend.position = "bottom") 
    
    
    p2 = ggplot(IF_ages_plots, aes(x = group, y = value, fill = reorder(Var1, value))) + geom_bar(stat = "identity") + 
      facet_grid(~ age, labeller=as_labeller(c(`2` = "2 mos", `4` = "4 mos", `6` = "6 mos",`8` = "8 mos"))) +
      labs(x = "Genotype", y = "Isoform Fraction (%)", fill = "Isoforms") + mytheme + 
      theme(legend.position = "bottom", legend.direction="vertical", strip.background = element_blank())
    
    
    if(type == "isoseq"){
      p1 = p1 + scale_fill_manual(values = c(wes_palette("IsleofDogs2")[1], wes_palette("IsleofDogs2")[4]), name = "Age (months)")
    }else{
      p1 = p1 + scale_fill_manual(values = c(wes_palette("IsleofDogs2")[1], wes_palette("IsleofDogs2")[2], wes_palette("IsleofDogs2")[3], 
                                             wes_palette("IsleofDogs2")[4]), name = "Age (months)") + theme(legend.position = "top")
    }
  }
  else{
    p = ggplot() + theme_void()
    p1 = ggplot() + theme_void()
    p2 = ggplot() + theme_void()
  }
  
  #if(type == "ONT_Targeted" & gene %!in% c("Ank1")){
  #  p2 = p2 + scale_fill_manual(values=c(plot_target_diff_col(gene,"ONT"),"Others" = alpha("Grey",0.8)))}
  #if(type == "Iso_Targeted" & gene %!in% c("Ank1")){
  #  p2 = p2 + scale_fill_manual(values=c(plot_target_diff_col(gene,"IsoSeq"),"Others" = alpha("Grey",0.8)))}
  
  output <- list(p,p1,p2)
  names(output) <- c("IFgen","IFage1","IFage2")
  return(output)
}  