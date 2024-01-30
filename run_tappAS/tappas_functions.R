### TappAS Differential Analysis #################################################################
input_tappasfiles <- function(tappas_input_dir){
  # read in files generated from TAPPAS
  tappasfiles <- list.files(path = tappas_input_dir, pattern = ".tsv", full.names = T)
  tappasfiles <- lapply(tappasfiles, function(x) read.table(x, sep = "\t", header = T))
  names(tappasfiles) <- list.files(path = tappas_input_dir, pattern = ".tsv")
  
  # InputExpression Matrix from TAPPAS of which transcripts are filtered during normalisation 
  # match the filtered transcripts with the associated gene uing the isoform id from classification file
  tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv <- 
    merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], 
          tappasfiles$tappAS_Transcripts_InputExpressionMatrix.tsv, by.x = "isoform", by.y = "Id")
  
  return(tappasfiles)
}

tappas_removediso <- function(filteredtappasfile){
  # number of transcripts that are filtered for statistical purposes
  dat = filteredtappasfile
  
  # tally of the number of transcripts filtered per gene 
  dat2 = dat %>% group_by(associated_gene, structural_category, Filtered) %>% tally()
  
  # plot the number of transcripts by gene only 
  p1 <- ggplot(dat2, aes(x = associated_gene, y = n, fill = Filtered)) + 
    geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms") + mytheme + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(labels = c("Retained","Removed due to low coverage"), values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[5])) + 
    theme(legend.position = c(0.85,0.8))
  
  p2 <- dat2 %>% filter(Filtered != "NO") %>% 
    ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + 
    geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms Removed") + 
    mytheme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = c(0.8,0.8))
  
  return(list(p1,p2))
  
}

# tappas_resultsanno
# prerequisite for plots
# aim1: annotate the tappas normalised expression file with the isoforms from the class file and join by phenotype
# aim2: deduce the gene expression by the sum of the expression of the filtered isoforms 
# output: norm_transcounts (normalised transcript counts) & GeneExp
tappas_resultsanno <- function(classification_file, tappas_normalised_expmatrix,phenotype){
  
  # Annotate the normalised expression matrix from tappas with the associated gene and sample phenotypes for plots 
  Norm_transcounts = 
    # annotate the tappas output normalised expression matrix of transcripts by the associated gene name and transcript name
    merge(classification_file [,c("isoform","associated_gene","associated_transcript","structural_category")], 
          tappas_normalised_expmatrix, by.x = "isoform", by.y = 0) %>% reshape2::melt() %>% 
    # annotate the samples to the phenotype 
    left_join(., phenotype, by = c("variable" = "sample")) %>% 
    # change factor levels for plots
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG")),
           structural_category=recode(structural_category, `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",`novel_in_catalog`="NIC"),
           Isoform = paste0(isoform,"_",structural_category))
  
  # Deduce gene expression from the sum of normalised transcript counts 
  GeneExp = Norm_transcounts %>% group_by(associated_gene,variable) %>% dplyr::summarise(Exp = sum(value)) %>%
    left_join(., phenotype, by = c("variable" = "sample"))
  
  output <- list(Norm_transcounts,GeneExp)
  names(output) <- c("Norm_transcounts","GeneExp")
  return(output)
}


# plot gene expression or by isoformID 
plot_mergedexp <- function(InputGene,IsoformID,GeneExp,Norm_transcounts){
  if (InputGene != "NA"){
    df <- GeneExp %>% filter(associated_gene == InputGene) 
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- InputGene
  }else if(IsoformID != "NA"){
    df <- Norm_transcounts  %>% filter(isoform == IsoformID) %>% left_join(., phenotype, by = c("variable" = "sample"))
    colnames(df)[6] <- "Exp"
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- paste0(df$associated_gene,": ",df$associated_transcript)
  }else{
    print("2 variables required")
  }
  
  p <- ggplot(df, aes(x = time, y = Exp, colour = group)) + geom_point() + 
    stat_summary(data=df, aes(x=time, y=Exp, group=group), fun="mean", geom="line", linetype = "dotted") + 
    labs(y = "Normalised Gene Expression", x = "Age (months)", title = paste0(plot_title,"\n\n")) + mytheme +
    # colours in the opposite direction of the levels 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT"))) + 
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  
  # subtitles
  if(IsoformID != "NA"){
    p <- p + labs(title = plot_title, subtitle = paste0(df$isoform,"\n\n"), y = "Normalised Isoform Expression") + theme(plot.subtitle = element_text(hjust = 0.5, size = 12,face = "italic"))
  }
  
  return(p)
}

# group_plots 
group_plots <- function(genegroup,plottype){
  myplots <- list()
  myplots2 <- list()
  
  if(length(genegroup) > 6){
    for(i in 1:6){myplots[[i]] <- plottype[[genegroup[[i]]]]}
    output1 = plot_grid(plotlist=myplots,labels = paste(letters[1:6]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
    count = 1
    for(i in 7:length(genegroup)){myplots2[[count]] <- plottype[[genegroup[[i]]]]; count = count + 1}
    output2 = plot_grid(plotlist=myplots2,labels = paste(letters[6:length(genegroup)]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
    output = list(output1,output2)
  }else{
    for(i in 1:length(genegroup)){myplots[[i]] <- plottype[[genegroup[[i]]]]}
    output = plot_grid(plotlist=myplots,labels = paste(letters[1:length(genegroup)]),label_size = 30, label_fontfamily = "CM Roman", nrow = 3,ncol = 2, scale = 0.9)
  }
  return(output)
}

group_plots_rnavsiso <- function(gene, plottype1, plottype2){
  myplots <- list()
  myplots[[1]] <- plottype1[[gene]] 
  myplots[[2]] <- plottype2[[gene]]
  output = plot_grid(plotlist=myplots,labels = c("a","b"),label_size = 30, label_fontfamily = "CM Roman",scale = 0.9)
  return(output)
}

# Differentially expressed genes 
# Input: tappassig with the sheet names "WholeIso_Geneexp" and "WholeRNA_Geneexp"
# Plots: 
# P1: venn diagram of genes that are differentially expressed between RNA+RNA(Isabel),Iso+RNA,Iso+Iso
tappas_genesig <- function(){
  
  # plot of the significant genes from Targeted Transcriptome using either Iso-Seq or RNA-Seq as expression input 
  p1 <- bind_rows(tappassig$TargetedIso_Genexp %>% mutate(type = "TargetedIso"),
                  tappassig$TargetedRNA_Genexp %>% mutate(type = "TargetedRNA")) %>% 
    ggplot(., aes(x = `...1`, y = `R-squared`, fill = type)) + geom_bar(position = position_dodge(preserve = "single"), stat = "identity") + mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(legend.position = "bottom") + labs(x = "", y =  expression(paste(R^2))) + 
    geom_hline(yintercept=0.5,linetype="dashed") + scale_fill_manual(name = "Transcript Expression Input", labels = c("Iso-Seq","RNA-Seq"), values = c(label_colour("isoseq"),label_colour("rnaseq")))
  
  return(p1)
}

plot_transexp <- function(InputGene,Norm_transcounts,type,name){
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  df$time <- as.factor(df$time)
  
  #p <-
  #  ggplot(df, aes(x = reorder(isoform,-value), y = value, colour = group, shape = time)) + #geom_boxplot() + 
  #  geom_jitter(size = 3, position = position_jitterdodge()) +
  #   stat_summary(data=df, aes(x=group, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
  #  mytheme + labs(x = "", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
  #  theme(strip.background = element_blank(), 
  #        plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
  #        panel.spacing = unit(2, "lines"), 
  #        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #  legend_theme +
  #  facet_grid(~structural_category,  scales = "free", space = "free") +
  #  scale_y_continuous(trans='log10') +
  #  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  p <- df %>% mutate(grouping = paste0(group,"_",time)) %>% group_by(grouping, isoform) %>% dplyr::summarise(structural_category,time,group, mean_exp = mean(value), .groups = 'drop') %>% mutate(group = factor(group, levels = c("TG","WT"))) %>% 
    ggplot(., aes(x = reorder(isoform,-mean_exp), y = mean_exp)) + geom_point(aes(colour = group, shape = time),size = 3) + mytheme + 
    labs(x = "", y = "Mean Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    theme(strip.background = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    legend_theme +
    facet_grid(~structural_category,  scales = "free", space = "free") +
    scale_y_continuous(trans='log10') +
    guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))
  
  
  if(type == "isoseq"){
    p <- p + scale_shape_manual(name = "Age", values=c(1, 16), label = c("2 months", "8 months")) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  } else {
    p <- p + scale_shape_manual(name = "Age (months)", values=c(1, 16, 2, 17)) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  }
  
  return(p)
}

plot_transexp_overtime <- function(InputGene,Norm_transcounts){
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  
  # linear regression for significant changes over time
  #for(isoform in unique(df$isoform)){
  #  cat("Processing",isoform,"\n")
  #  df1 <- df[df$isoform == isoform & df$group == "WT",]
  #  df1_WTmean <- df1 %>% group_by(time) %>% summarise(mean_exp = mean(value))
  #  df1_TG <- df[df$isoform == isoform & df$group == "TG",]
  #  df2 <- merge(df1_TG, df1_WTmean, by = "time") %>% mutate(diff = abs(value - mean_exp))
  #print(summary(lm(diff~0 + time,df2)))
  #}
  
  #df_WT <-  df[df$group == "WT",]
  #df_WTmean <- df %>% group_by(time, isoform) %>% summarise(mean_exp = mean(value),  .groups = 'drop')
  #df_TG <- df[df$group == "TG",]
  #df2 <- merge(df_TG, df_WTmean, by = c("time","isoform")) %>% mutate(diff = abs(value - mean_exp))
  
  #p1 <- ggplot(df2, aes(x = time, y = diff, colour = isoform)) + geom_point() +
  #  stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
  #  scale_y_continuous(trans = 'log10') + mytheme + labs(x = "Age (months)", y = "Fold Change of Isoform Expression \n (TG - WT)") + theme(legend.position = "right")
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point() + 
    facet_grid(~group,scales = "free", space = "free") +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = "Age (months)", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    theme(strip.background = element_blank(), legend.position = "right",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) 
  return(p)
}


