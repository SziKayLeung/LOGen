## ---------- Plot gene and transcript expression -----------------


relabel_case_control <- function(df){
  
  levels(df$group)[levels(df$group) %in% c("Control","CONTROL")] <- label_group("Control")
  levels(df$group)[levels(df$group) %in% c("Case","CASE")]   <- label_group("Case")
  df <- df %>% mutate(group = factor(group, levels = c(label_group("Control"),label_group("Case"))))
  return(df)
  
}

time_case_boxplot <- function(normalised_counts, transcript){
  
  df <- normalised_counts %>% filter(isoform == transcript)
  df <- relabel_case_control(df)
  df$time <- as.factor(df$time)
  isoID <- word(transcript,c(3),sep = fixed("."))
  cate <- unique(df$structural_category)
  
  p <- ggplot(df, aes(x = group, y = normalised_counts)) + geom_boxplot(outlier.shape = NA) +
    geom_point(position="jitter",aes(color = time), size = 3) +
    mytheme +
    labs(x = "Genotype", y = "Normalized counts",
         title = paste0("LR.", unique(df$associated_gene),".",isoID),
         subtitle = paste0(df$associated_transcript," (",cate, ")")) +
    scale_colour_manual(name = "Age (months)",
                        values = c(wes_palette("Darjeeling2")[[5]], wes_palette("Zissou1")[[1]],
                                   wes_palette("Zissou1")[[3]],wes_palette("Zissou1")[[5]]))
  
  return(p)
}

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
    mytheme + labs(x = "", y = "Normalized Isoform Expression",title = plot_title) +
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
  group1 <- unique(Norm_transcounts$group)[1]
  group2 <- unique(Norm_transcounts$group)[2]
  
  p <- ggplot(dat, aes(x = group, y = value, fill = group)) + geom_boxplot() + 
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    labs(x = "", y = "Isoform Expression", 
         title = paste0(gene, ":",transcript,"")) + mytheme + 
    scale_fill_manual(values = c(label_colour(group1),label_colour(group2))) + 
    theme(legend.position = "none")
  
  return(p)
}

plot_trans_exp_individual_overtime <- function(transcript, Norm_transcounts,type=transcript){
  print(transcript)
  dat <- Norm_transcounts %>% filter(isoform == transcript)
  dat$time <- as.factor(dat$time)
  
  groups = levels(Norm_transcounts$group)
  if("Control" %in% groups){
    dat <- dat %>% mutate(group = factor(group, levels = c("Control",groups[!groups %in% "Control"])))
  }else{
    dat <- dat %>% mutate(group = factor(group, levels = c("CONTROL",groups[!groups %in% "CONTROL"]))) 
  }
  
  group1 <- levels(dat$group)[[1]]
  group2 <- levels(dat$group)[[2]]
  
  # title
  gene <- dat$associated_gene[1]
  associated_transcript <- dat$associated_transcript[1]
  if(type == "transcript"){
    title = paste0(associated_transcript, " (", transcript,")")
    subtitle = gene
  }else{
    title = gene
    subtitle = NULL
  }

  p <- ggplot(dat, aes(x = time, y = normalised_counts, colour = group)) + geom_point(size = 3) +
    stat_summary(data=dat, aes(x=time, y=normalised_counts, group=group), fun ="mean", geom="line", linetype = "dotted") +
    labs(x = "Age (months)", y = "Normalized counts", 
         title = title, subtitle = subtitle) + mytheme + 
    scale_colour_manual(values = c(label_colour(group1),label_colour(group2))) 
  
  return(p)
}

plot_all_gene_exp <- function(){
  df1 = targetedtappas_isoexp$GeneExp %>% filter(associated_gene %in% TargetGene)
  df2 = targetedtappas_rnaexp$GeneExp %>% filter(associated_gene %in% TargetGene)
  
  p <- rbind(df1 %>% mutate(dataset = "IsoSeq"), df2 %>% mutate(dataset = "RNASeq"))%>% 
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("Ctrl", "AD"))) %>%
    ggplot(., aes(x = group, y = log10(Exp), fill = group)) + geom_boxplot() + facet_grid(dataset ~ associated_gene) + 
    labs(y = "Normalized gene expression (Log10)", x = "Phenotype") + 
    scale_fill_manual(values = c(label_colour("AD"),label_colour("Control")), labels = c("AD","Control"),
                      name = "Phenotype") + 
    theme(legend.position = "top")  
  
  return(p)
}


subsetNormCounts <- function(inputGene, normCounts,design="time_series",show="all",rank=5, isoSpecific=NULL){
  
  df <-  normCounts  %>% filter(associated_gene == inputGene) 

  if(design=="time_series"){
    df$time <- as.factor(df$time)
  }else if(design != "multiple_case_control"){
    df <- relabel_case_control(df)
  }

  if(show == "toprank"){
    cat("Keeping only the top-ranked", rank, "most-expressed isoforms across entire dataset\n")
    topranked_expression <- df %>% group_by(isoform) %>% tally(normalised_counts) %>% arrange(-n)
    keepiso <- as.character(topranked_expression$isoform[1:rank])
    
    if(!is.null(isoSpecific)){
      cat("Keeping also isoforms", isoSpecific, "\n")
      keepiso <- c(keepiso, isoSpecific)
    }
    
    df <- df %>% filter(isoform %in% keepiso) 
  
  }else if(show == "specific"){
    df <- df %>% filter(isoform %in% isoSpecific) 
  }
  
  return(df)
}

replace_pbID <- function(PbID, gene){
  isoform <- word(PbID, c(3), sep = fixed("."))
  return(paste0("LR.",gene,".",isoform))
}


plot_transexp_overtime <- function(inputGene, normCounts,design="time_series",show="all",rank=5, isoSpecific=NULL, plotTitle=NULL,setorder=NULL,classfiles=NULL){
  
  df <- subsetNormCounts(inputGene,normCounts,design=design,show=show,rank=rank,isoSpecific=isoSpecific)
  df$LRID <- apply(df, 1, function(x) replace_pbID (x[["isoform"]], x[["associated_gene"]]))
  
  if(!is.null(classfiles)){
    df <- merge(df,classfiles[,c("isoform","structural_category")], by = "isoform")
    if(!"structural_category" %in% colnames(df)){
      df <- merge(df,classfiles[,c("isoform","structural_category")], by = "isoform")
    }
    df <- df %>% dplyr::mutate(LRID_struc = paste0(LRID," (", structural_category, ")"))
  }
  
  if(!is.null(setorder) & design != "multiple_case_control"){
    df$group <- factor(df$group, levels = setorder,
                       labels = c(label_group(setorder[1]),label_group(setorder[2])))
  }else{
    df$group <- factor(df$group, levels = setorder)
  }
  
  if(design != "time_series"){
  
    p <- ggplot(df, aes(x = group, y = normalised_counts, colour = LRID)) + geom_boxplot(outlier.shape = NA) + 
      geom_point(position = position_jitterdodge()) +
      mytheme + labs(x = " ", y = "Normalized counts", title = inputGene) +
      theme(strip.background = element_blank(), legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines"),
            text=element_text(size=16)) +
      guides(colour=guide_legend(ncol=3,bycol=TRUE))
    
    return(p)
    
  }else if(design == "time_series"){
    
    if(length(unique(df$isoform)) == 1){
      p <- ggplot(df, aes(x = time, y = normalised_counts, colour = group)) + geom_point(size = 3) +
        stat_summary(data=df, aes(x=time, y=normalised_counts, group=group), fun ="mean", geom="line", linetype = "dotted", size = 1.5) +
        mytheme + labs(x = "Age (months)", y = "Normalized counts", title = inputGene) +
        scale_colour_manual(values = c(label_colour(levels(df$group)[1]),label_colour(levels(df$group)[2]))) 
    }else{
      p <- ggplot(df, aes(x = time, y = normalised_counts, colour = LRID)) + geom_point(size = 3) + 
        facet_grid(~group,scales = "free", space = "free") +
        mytheme + labs(x = "Age (months)", y = "Normalized counts", title = inputGene)  +
        stat_summary(data=df, aes(x=time, y=normalised_counts, group=LRID), fun ="mean", geom="line", linetype = "dotted", size = 1.5)
    }
    p <- p + theme(strip.background = element_blank(), 
            text=element_text(size=16),
            plot.title = element_text(hjust = 0, size = 16),
            panel.spacing = unit(1, "lines"),
            legend.position = c(0.4,0.8))
  }
  
  if(show == "specific" & length(isoSpecific) == 1){
    if(!is.null(classfiles)){
      struc = unique(df$structural_category)
      isotitle = paste0("LR.", inputGene,".",word(isoSpecific,c(3),sep=fixed("."))," (",struc,")")
    }else{
      isotitle = paste0("LR.", inputGene,".",word(isoSpecific,c(3),sep=fixed("."))) 
    }
    p <- p + labs(title = isotitle) + theme(legend.position = "None")
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
    labs(x = "", y = "Normalized counts", title = transcript) +
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
    mytheme + labs(x = " ", y = "Normalized Counts",title = plot_title) +
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


