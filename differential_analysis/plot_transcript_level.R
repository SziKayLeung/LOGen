## ---------- Plot gene and transcript expression -----------------

# plot transcript expression of all samples as box plot by phenotype
# top 10 most abundantly expressed transcript 
# or all transcripts of all samples 
plot_trans_exp <- function(InputGene, Norm_transcounts, type, name){
  plot_title <- paste0(InputGene,"\n",name)
  df <-  Norm_transcounts %>% filter(associated_gene == InputGene) 
  
  groups = levels(Norm_transcounts$group)
  if("Control" %in% groups){
    df <- df %>% mutate(group = factor(group, levels = c("Control",groups[!groups %in% "Control"])))
  }else{
    df <- df %>% mutate(group = factor(group, levels = c("CONTROL",groups[!groups %in% "CONTROL"]))) 
  }
  
  if (type == "top10"){
    # show only the first highest expressed isoforms
    top10_isoforms = df %>% group_by(isoform) %>% tally(value) %>% arrange(-n) %>% as.data.frame() %>% .[1:10,1]
    df <- df %>% filter(isoform %in% top10_isoforms) 
  }
  
  p <- ggplot(df, aes(x = reorder(isoform,-value), y = value, colour = group)) + geom_boxplot() + 
    mytheme + labs(x = "", y = "Normalised Isoform Expression",title = plot_title) +
    scale_colour_manual(values = c(label_colour("Case"),label_colour("Control")), name = "Phenotype",
                        c(label_group("Control"),label_group("Case"))) +
    theme(strip.background = element_blank(), legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  if(nrow(df) > 1){
    p = p + facet_grid(~structural_category,  scales = "free", space = "free")
  }
  return(p)
}

# plot individual transcripts/isoforms
plot_trans_exp_individual <- function(transcript, Norm_transcounts){
  print(transcript)
  dat <- Norm_transcounts %>% filter(isoform == transcript)
  gene <- dat$associated_gene[1]
  group1 <- unique(dat$group)[1]
  group2 <- unique(dat$group)[2]
  
  p <- ggplot(dat, aes(x = group, y = value, fill = group)) + geom_boxplot() + 
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    labs(x = "", y = "Isoform Expression", 
         title = paste0(gene, ":", "\n",transcript,"")) + mytheme + 
    scale_fill_manual(values = c(label_colour(group1),label_colour(group2))) + 
    theme(legend.position = "none")
  
  return(p)
}

plot_trans_exp_individual_overtime <- function(transcript, Norm_transcounts){
  print(transcript)
  dat <- Norm_transcounts %>% filter(isoform == transcript)
  gene <- dat$associated_gene[1]
  group1 <- unique(dat$group)[1]
  group2 <- unique(dat$group)[2]
  
  p <- ggplot(dat, aes(x = time, y = value, colour = group)) + geom_point(size = 3) +
    labs(x = "Age (months)", y = "Isoform Expression", 
         title = paste0(gene, ":", "\n",transcript,"")) + mytheme + 
    scale_colour_manual(values = c(label_colour(group1),label_colour(group2))) + 
    theme(legend.position = "none")
  
  return(p)
}

plot_all_gene_exp <- function(){
  df1 = targetedtappas_isoexp$GeneExp %>% filter(associated_gene %in% TargetGene)
  df2 = targetedtappas_rnaexp$GeneExp %>% filter(associated_gene %in% TargetGene)
  
  p <- rbind(df1 %>% mutate(dataset = "IsoSeq"), df2 %>% mutate(dataset = "RNASeq"))%>% 
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("Ctrl", "AD"))) %>%
    ggplot(., aes(x = group, y = log10(Exp), fill = group)) + geom_boxplot() + facet_grid(dataset ~ associated_gene) + 
    labs(y = "Normalised gene expression (Log10)", x = "Phenotype") + 
    scale_fill_manual(values = c(label_colour("AD"),label_colour("Control")), labels = c("AD","Control"),
                      name = "Phenotype") + 
    theme(legend.position = "top")  
  
  return(p)
}


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

plot_transexp_overtime <- function(InputGene, Norm_transcounts, name, type, dataset, difftrans){
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(name," ", InputGene)
  
  if(type == "case_control"){
    # facet labels 
    df2 <- df
    levels(df2$group)[levels(df2$group) %in% c("Control","CONTROL")] <- label_group("Control")
    levels(df2$group)[levels(df2$group) %in% c("Case","CASE")]   <- label_group("Case")
    
    p <- ggplot(df2, aes(x = time, y = value, colour = Isoform)) + geom_point() + 
      facet_grid(~group, scales = "free", space = "free") +
      stat_summary(data=df2, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
      mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
      theme(strip.background = element_blank(), legend.position = "bottom",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) +
      guides(colour=guide_legend(ncol=3,bycol=TRUE))
    
  }else if(type == "time_series_sig"){
    levels(df$group)[levels(df$group) %in% c("Control","CONTROL")] <- label_group("Control")
    levels(df$group)[levels(df$group) %in% c("Case","CASE")]   <- label_group("Case")
    df <- df %>% mutate(group = factor(group, levels = c(label_group("Control"),label_group("Case"))))
    
    if(InputGene %in% c("Gfap","C4b","Ctsd","Gatm","H2-D1","Padi2","Cd34","Ubqln1")){
      if(dataset == "isoseq"){
        sigiso = c(difftrans[difftrans$associated_gene == InputGene, "isoform"])$isoform
        df = df %>% mutate(sig = ifelse(isoform %in% sigiso,Isoform,"NA")) %>% filter(sig != "NA")
      }else{
        sigisoiso = c(difftrans[difftrans$associated_gene == InputGene, "isoform"])$isoform
        df = df %>% mutate(Isoform = ifelse(isoform %in% sigisoiso,Isoform,"RNA-specific"))
      }
    }else{
      difftrans = difftrans %>% filter(associated_gene == InputGene) %>% arrange(`p-value`)
      # top 5 differentially-ranked isoforms 
      df <- df %>% filter(isoform %in% difftrans$isoform[1:3])
      
    }
    
    
    if(InputGene == "Gfap" & dataset == "rnaseq"){
      df <- df %>% mutate(Isoform = factor(Isoform, levels = c("PB.2972.13_NIC","PB.2972.16_FSM","PB.2972.36_FSM","RNA-specific")))
      colours = c("#F8766D","#00BFC4","#C77CFF",alpha("grey",0.6))
    } else if(InputGene == "Gfap"){
      colours = c("#F8766D","#7CAE00","#00BFC4","#C77CFF",alpha("grey",0.6))
    } else if(InputGene == "Ubqln1"){
      colours = c("#F8766D","#00BFC4",alpha("grey",0.6))
    }else{
      colours = c("#F8766D",alpha("grey",0.6))
    }
    
    if(nrow(df) != 0){
      p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point(size = 3) + 
        facet_grid(~group,scales = "free", space = "free") +
        stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted", size = 1.5) +
        mytheme + labs(x = "Age (months)", y = "Normalised counts",title = plot_title) +
        theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
              panel.spacing = unit(2, "lines"),
              legend.position = c(0.4,0.8)) 
      
      if(dataset == "isoseq"){
        p <- p + scale_colour_discrete(name = "Isoform")
      }else{
        p <- p + scale_colour_manual(values = colours) + theme(legend.position = "none")
      }
      
    }else{
      p <- ggplot() + theme_void()
    }
    
  }
  
  return(p)
}

twocate_plot_transexp_overtime <- function(InputGene,Norm_transcounts, name){
  print(InputGene)
  p <- plot_transexp_overtime(InputGene,Norm_transcounts,name) + 
    scale_x_discrete(limits = c(0,1), labels = c(label_group("0_timeanalysis"),label_group("1_timeanalysis"))) + 
    labs(x = "")
  
  return(p)
}


twocate_plot_transexp <- function(transcript, exp_matrix, phenotype){
  
  df <- exp_matrix %>% rownames_to_column(., var = "isoform") %>% 
    reshape2::melt(id = "isoform") %>% 
    filter(isoform == transcript) %>%
    right_join(., phenotype, by = c("variable" = "sample")) 
  
  p <- ggplot(df, aes(x = as.factor(time), y = value, fill = group)) + geom_boxplot() +
    scale_x_discrete(labels = c(label_group("0_timeanalysis"),label_group("1_timeanalysis"))) + 
    labs(x = "", y = "Normalised counts", title = transcript) +
    scale_fill_discrete(labels = c(label_group("Case_timeanalysis"),label_group("Control_timeanalysis"))) +
    mytheme
  
  return(p)
}


plot_transexp_overtime_filter <- function(InputGene,Norm_transcounts,original_norm, phenotype,name){
  
  exp_stats <- manual_filter_isoforms(original_norm, Norm_transcounts, InputGene,phenotype)
  df <- Norm_transcounts  %>% filter(associated_gene == InputGene) %>% filter(isoform %in% row.names(exp_stats$exp_filtered))
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  
  p <- ggplot(df, aes(x = group, y = value, colour = Isoform)) + geom_point() + 
    stat_summary(data=df, aes(x=group, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = " ", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), legend.position = "bottom",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines")) +
    guides(colour=guide_legend(ncol=3,bycol=TRUE)) 
  
  return(p)
}   

# abundance file before normalisation
plot_raw_expression <- function(FL_reads, transcript){
  p <- FL_reads[rownames(FL_reads) == transcript, ] %>% reshape2::melt() %>% 
    dplyr::rename(sample = variable) %>%
    full_join(., phenotype$ont, by = "sample") %>% 
    ggplot(., aes(x = group, y = value)) + geom_boxplot() + 
    mytheme + labs(x = " ", y = "FL Reads")
  
  return(p)
}


