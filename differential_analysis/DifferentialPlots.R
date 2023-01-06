## ---------- Script -----------------
##
## Purpose: input variables for differential analysis of Iso-Seq targeted mouse transcriptome datasets 
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)





generate_diff_plots <- function(genelist,tappasinput,type,name,...){
  if(type == "Gene"){
    plist <- lapply(lapply(genelist, function(gene) plot_mergedexp(gene,"NA",tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],name)),ggplotGrob)
  }else if(type == "IsoGene_Targeted"){
    plist <- lapply(lapply(genelist, function(gene) plot_mergedexp(gene,"NA",tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],name) + 
                      labs(y = "Iso-Seq Normalised Counts", title = paste0(gene,"\n","Gene Expression","\n"))+
                        theme(axis.text.y = element_text(angle = 90, hjust = 0.5))),ggplotGrob)
  }else if(type == "ONTGene_Targeted"){
    plist <- lapply(lapply(genelist, function(gene) plot_mergedexp(gene,"NA",tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],name) + 
                      labs(y ="ONT Normalised Counts", title = paste0(gene,"\n","Gene Expression","\n")) +
                      theme(axis.text.y = element_text(angle = 90, hjust = 0.5))),ggplotGrob)
  }else if(type == "Transcript"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp(gene,tappasinput[["Norm_transcounts"]],"isoseq",name)),ggplotGrob)
  }else if(type == "Transcript_Rnaseq_Targeted"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp(gene,tappasinput[["Norm_transcounts"]],"rnaseq",name)),ggplotGrob)
  }else if(type == "Transcript_Isoseq Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],"isoseq",name)),ggplotGrob)
  }else if(type == "Transcript_Rnaseq Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],"rnaseq",name)),ggplotGrob)
  }else if(type == "Transcript Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_general_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],...,name)),ggplotGrob)
  }else if(type == "IF_Isoseq"){
    plist <- lapply(lapply(genelist, function(gene) IF_plot(gene,tappasinput[["gene_transcripts.tsv"]],tappasinput[["input_norm"]],"Iso_Targeted")[[3]]),ggplotGrob)
  }else if(type == "IF_ONT"){
    plist <- lapply(lapply(genelist, function(gene) IF_plot(gene,tappasinput[["gene_transcripts.tsv"]],tappasinput[["input_norm"]],"ONT_Targeted")[[3]]),ggplotGrob)
  }else{
    print("Type Required")
  }
  
  names(plist) = genelist
  return(plist)
}

### Differential Gene Expression #################################################################
tabulate_diffgene <- function(gene, GeneExp, Interaction){
  print(gene)
  df <- GeneExp %>% filter(associated_gene == gene) 
  
  meanexp = df %>% group_by(time, group) %>% summarise_at(vars(Exp), funs(mean(., na.rm=TRUE))) %>% as.data.frame() %>% 
    mutate(groupings = paste0(time,group)) %>% .[,c("Exp","groupings")] %>% spread(., groupings, Exp) %>% 
    mutate(log2fc = log2(`8CASE`/`2CASE`)) #%>% .[,c(5,2,4,1,3)]
  
  FDR = Interaction[Interaction$...1 == gene,c("p-value","R-squared")]
  
  if(nrow(FDR) == 0){FDR = data.frame(`p-value` = as.numeric("0"), `R-squared` = as.numeric("0"))}
  
  # 3 significant factor 
  meanexp = cbind(FDR, meanexp)
  meanexp = signif(meanexp, digits = 3)
  colnames(meanexp)[1] = "p-value"
  colnames(meanexp)[2] = "R-squared"
  
  return(meanexp)
}

all_tab_gene <- function(GeneExp, diffgene){
  # empty list
  meanexp_output = list()
  meanexp_output = lapply(TargetGene, function(x) tabulate_diffgene(x, GeneExp, diffgene))
  names(meanexp_output) = TargetGene
  meanexp_output <- do.call("rbind",meanexp_output) 
  return(meanexp_output)
}


#transExp = IsoExp$Norm_transcounts
#trans = "PB.1691.13"
simple_diff_stats <- function(trans, transExp, type){
  
  df <- transExp %>% filter(isoform == trans) %>% mutate(Exp = value)
  
  meanexp = df %>% group_by(time, group) %>% summarise_at(vars(Exp), funs(mean(., na.rm=TRUE))) %>% as.data.frame() %>% 
    mutate(groupings = paste0(time,group)) %>% .[,c("Exp","groupings")] %>% spread(., groupings, Exp) 
  
  if(nrow(meanexp != 0)){meanexp = meanexp %>% mutate(log2fc_pheno = log2(`8TG`/`8WT`),log2fc_age = log2(`8TG`/`2TG`))} 
  
  if(type == "isoseq"){
    meanexp = meanexp[,c(5,2,4,1,3)]
  }
  
  return(meanexp)
}

all_tab_iso <- function(gene, Norm_transcounts, difftrans, type){

  sig = difftrans[difftrans$associated_gene == gene,]
  sig_stat = do.call(rbind,lapply(sig$isoform, function(x) tabulate_difftrans(x,Norm_transcounts,difftrans,type))) 
  
  return(sig_stat)
}

final_targeted_difftrans <- function(Gene){
  print(Gene)
  isoSig = tTransStat$iso[tTransStat$iso$associated_gene == Gene,] %>% 
    mutate(dataset = "Iso-Seq", value = paste0(log2fc," (", `p-value`,",", `R-squared`,")")) %>% 
    dplyr::select(associated_gene, value, dataset)
  
  ontSig = tTransStat$ont[tTransStat$ont$associated_gene == Gene,] %>% 
    arrange(`p-value`) %>% 
    mutate(dataset = "ONT", value = paste0(log2fc," (", `p-value`,",", `R-squared`,")")) %>% 
    dplyr::select(associated_gene, value, dataset) 
  
  if(nrow(ontSig) >= 3){ontSig = ontSig[1:3,]}
  output = rbind(isoSig,ontSig) %>% rownames_to_column(., var = "isoform") %>% spread(., dataset, value) 
  
  if(nrow(isoSig) > 0){
    output = output %>% mutate(`Iso-Seq` = paste0(isoform," : ", `Iso-Seq`)) %>% mutate(`ONT` = paste0(isoform," : ", `ONT`))
  }else{
    output = output %>% mutate(`Iso-Seq` = "NA")  %>% mutate(`ONT` = paste0(isoform," : ", `ONT`))
  }
  
  return(output)
}



tabulate_difftrans <- function(trans, transExp, All_FDR, type){
  
  meanexp = simple_diff_stats(trans,transExp,type)
  FDR = All_FDR[All_FDR$isoform == trans,c("associated_gene","p-value","R-squared")]
  meanexp = cbind(FDR, meanexp)
  
  # 3 significant factor 
  meanexp = cbind(meanexp[,1],signif(meanexp[,2:ncol(meanexp)], digits = 3))
  colnames(meanexp)[1] = c("associated_gene")
  rownames(meanexp) = trans
  
  return(meanexp)
}

# plot gene expression or by isoformID 
plot_mergedexp <- function(InputGene,IsoformID,GeneExp,Norm_transcounts,type){
  if (InputGene != "NA"){
    df <- GeneExp %>% filter(associated_gene == InputGene) 
    df$group <- factor(df$group, levels = c("CONTROL","CASE"))
    plot_title <- paste0(InputGene,"\n",type)
    cat("Mean Expression across groups for", InputGene,"\n")
    
  }else if(IsoformID != "NA"){
    df <- Norm_transcounts  %>% filter(isoform == IsoformID) #%>% left_join(., phenotype$ont, by = c("variable" = "sample"))
    colnames(df)[6] <- "Exp"
    df$group <- factor(df$group, levels = c("WT","TG"))
    plot_title <- paste0(df$associated_gene,": ",df$associated_transcript)
  }else{
    print("2 variables required")
  }
  
  p <- ggplot(df, aes(x = time, y = Exp, colour = group)) + geom_point(size = 3) + 
    stat_summary(data=df, aes(x=time, y=Exp, group=group), fun="mean", geom="line", linetype = "dotted") + 
    labs(y = "Normalised Counts", x = "Age (months)", title = paste0(plot_title,"\n\n")) + mytheme +
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
# Input: tappassig with the all the sheets
# Plots: 
# P1: venn diagram of genes that are differentially expressed between RNA+RNA(Isabel),Iso+RNA,Iso+Iso
tappas_genesig <- function(){
  
  ## Genes that are already filtered by significance (p <0.05 and R > 0.5)
  # Genes associated with interaction effects (results from using RNA-Seq or Iso-Seq as expression)
  WholeRNA_Interaction = c(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1)
  
  WholeIso_Interaction = c(gene_sigs_WholeIso_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 5 Interaction`$...1,
                           gene_sigs_WholeIso_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 7 Interaction`$...1)
  
  # Genes associated with genotype effects
  WholeRNA_Genotype = gene_sigs_WholeRNA_lst$models$`Model 1 Genotype`$...1
  WholeIso_Genotype = gene_sigs_WholeIso_lst$models$`Model 1 Genotype`$...1
  
  
  cat("Nummber of DEG from Isabel's supp, Genotype Effect:", nrow(Isabel_gene_Tg4510AgeGenotypeDEG),"\n")
  cat("Nummber of DEG from Isabel's supp, Genotype & Age Effect:", nrow(Isabel_gene_Tg4510GenotypeDEG),"\n")
  
  # first column of the input table = gene list 
  p1 <- venn.diagram(x = list(Isabel_gene_Tg4510AgeGenotypeDEG$Gene, WholeRNA_Interaction, WholeIso_Interaction), 
                     category.names = c("Reference genome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n Iso-Seq Expression"), 
                     filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC","#CBD5E8"), main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                     cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                     print.mode = "raw")
  
  p2 <- venn.diagram(x = list(Isabel_gene_Tg4510GenotypeDEG$Gene, WholeRNA_Genotype, WholeIso_Genotype), 
                     category.names = c("Reference genome, \n RNA-Seq Expression",
                                        "Iso-Seq transcriptome, \n RNA-Seq Expression","Iso-Seq transcriptome, \n Iso-Seq Expression"), 
                     filename = NULL, output=TRUE, lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC","#CBD5E8"), 
                     main = "\n", cex = 1,fontface = "bold",fontfamily = "ArialMT",
                     cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-30, 0, 30), cat.dist = c(0.055, 0.055, 0.065),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n"
                     print.mode = "raw")
  
  output <- list(p1,p2)
  names(output) <- c("p1","p2")
  return(output)
  
}

RNASeq_IsoSEQ_Geneexp <- function(){
  DEA_uIiso = setdiff(WholeIso_Interaction,WholeRNA_Interaction)
  DEA_uRNA = setdiff(WholeRNA_Interaction,WholeIso_Interaction)
  DEA_commonRNAIso = intersect(WholeRNA_Interaction,WholeIso_Interaction)
  
  plot_mergedexp("Cd34","NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Cd34","NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Thra","NA",wholetappas_isoexp$GeneExp,wholetappas_isoexp$Norm_transcounts,"Iso-Seq Expression")
  plot_mergedexp("Thra","NA",wholetappas_rnaexp$GeneExp,wholetappas_rnaexp$Norm_transcounts,"Iso-Seq Expression")
  
  tappasiso$gene_matrix_mean <- apply(tappasiso$gene_matrix.tsv,1,mean) %>% reshape2::melt() %>% `colnames<-`(c("meanIsoGeneExp"))
  wholetappas_isoexp$GeneExp <- merge(wholetappas_isoexp$GeneExp,tappasiso$gene_matrix_mean %>% rownames_to_column(var = "associated_gene"), by = "associated_gene")
  
  tappasrna$gene_matrix_mean <- apply(tappasrna$gene_matrix.tsv,1,mean) %>% reshape2::melt() %>% `colnames<-`(c("meanRNAGeneExp"))
  wholetappas_rnaexp$GeneExp <- merge(wholetappas_rnaexp$GeneExp,tappasrna$gene_matrix_mean %>% rownames_to_column(var = "associated_gene"), by = "associated_gene")
  
  meanGeneExp = merge(distinct(wholetappas_rnaexp$GeneExp[,c("associated_gene","meanRNAGeneExp")]),
                      distinct(wholetappas_isoexp$GeneExp[,c("associated_gene","meanIsoGeneExp")]), by = "associated_gene")
  for(i in 1:nrow(meanGeneExp)){
    meanGeneExp$type[[i]] = if(meanGeneExp$associated_gene[[i]] %in% DEA_uIiso){"IsoSeq_Only"
    }else if(meanGeneExp$associated_gene[[i]] %in% DEA_uRNA){"RNASeq Only"
    }else if(meanGeneExp$associated_gene[[i]] %in% DEA_commonRNAIso){"Common DEA"}else{"Not DEA"}}
  
  meanGeneExp %>% filter(type != "Not DEA") %>% 
    ggplot(., aes(x = meanRNAGeneExp, y = meanIsoGeneExp, colour = type)) + geom_point() + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") 
  
  rbind(tappassiggene$WholeRNA_Genexp[tappassiggene$WholeRNA_Genexp$...1 %in% DEA_uRNA,c("R-squared")] %>% mutate(type = "RNASeq"),
        tappassiggene$WholeRNA_Genexp[tappassiggene$WholeRNA_Genexp$...1 %in% DEA_commonRNAIso,c("R-squared")] %>% mutate(type = "Both")) %>% 
    ggplot(., aes(`R-squared`, fill = type)) + geom_density(alpha = 0.2)
  
  rbind(tappassiggene$WholeIso_Genexp[tappassiggene$WholeIso_Genexp$...1 %in% DEA_uIiso,c("R-squared")] %>% mutate(type = "IsoSeq"),
        tappassiggene$WholeIso_Genexp[tappassiggene$WholeIso_Genexp$...1 %in% DEA_commonRNAIso,c("R-squared")] %>% mutate(type = "Both")) %>% 
    ggplot(., aes(`R-squared`, fill = type)) + geom_density(alpha = 0.2)
  
}

### Differential Transcript Expression #################################################################
tabulate_all_diffgenes <- function(GeneExp){
  Allsig_output = list()
  # Interaction genes only (365 genes)
  Allsig_output = lapply(Interaction$...1, function(x) tabulate_diffgene(x, GeneExp, Interaction))
  names(Allsig_output) = Interaction$...1
  Allsig_output <- do.call("rbind",Allsig_output) %>% mutate(direction = ifelse(log2fc > 0, "Up", "Down")) 
  
  # binomial test
  threshold <- 0.05
  ntotal <- nrow(Allsig_output) # number of sig genes I identified
  nupregulated <- nrow(Allsig_output[Allsig_output$direction == "Up",])
  test <- binom.test(nupregulated, ntotal, 0.5)
  pvalue <- test[3]
  pvalue <- signif(as.numeric(pvalue), 3)
  
  ## Volcano plot
  #Allexp_output = list()
  #Allexp_output = lapply(tappassiggene$WholeIsoAll_Genexp$...1, function(x) tabulate_diffgene(x, IsoExp$GeneExp, tappassiggene$WholeIsoAll_Genexp))
  #names(Allexp_output) = tappassiggene$WholeIsoAll_Genexp$...1
  #Allexp_output <- do.call("rbind",Allexp_output) %>% mutate(direction = ifelse(log2fc > 0, "Up", "Down")) 
  #Allexp_output = Allexp_output %>% mutate(sig = ifelse(`p-value` < 0.05 & `R-squared` > 0.5, "sig","not-sig"))
  #Allexp_output %>% filter(log2fc < 2.5) %>% 
  #  ggplot(., aes(x=log2fc, y=-log10(`p-value`), colour = sig)) + geom_point() + mytheme
}

tabulate_all_difftrans <- function(){
  Allsig_output = list()
  # Interaction genes only (365 genes)
  Allsig_output = lapply(Interaction_Difftrans$isoform, function(x) tabulate_difftrans(x, IsoExp$Norm_transcounts))
  names(Allsig_output) = Interaction_Difftrans$isoform
  Allsig_output <- do.call("rbind",Allsig_output) %>% mutate(direction = ifelse(log2fc > 0, "Up", "Down")) 
  
  Allsig_output %>% group_by(direction) %>% tally()
  
  # binomial test
  threshold <- 0.05
  ntotal <- nrow(Allsig_output) # number of sig genes I identified
  nupregulated <- nrow(Allsig_output[Allsig_output$direction == "Up",])
  test <- binom.test(nupregulated, ntotal, 0.5)
  pvalue <- test[3]
  pvalue <- signif(as.numeric(pvalue), 3)
  
  ## Volcano plot
  #Allexp_output = list()
  #Allexp_output = lapply(tappassiggene$WholeIsoAll_Genexp$...1, function(x) tabulate_diffgene(x, IsoExp$GeneExp, tappassiggene$WholeIsoAll_Genexp))
  #names(Allexp_output) = tappassiggene$WholeIsoAll_Genexp$...1
  #Allexp_output <- do.call("rbind",Allexp_output) %>% mutate(direction = ifelse(log2fc > 0, "Up", "Down")) 
  #Allexp_output = Allexp_output %>% mutate(sig = ifelse(`p-value` < 0.05 & `R-squared` > 0.5, "sig","not-sig"))
  #Allexp_output %>% filter(log2fc < 2.5) %>% 
  #  ggplot(., aes(x=log2fc, y=-log10(`p-value`), colour = sig)) + geom_point() + mytheme
}

#gene = "Trpa1"
#Norm_transcounts = IsoExp$Norm_transcounts
#cf = class.files$iso
#type = "targeted"
draw_heatmap_gene <- function(gene, cf, Norm_transcounts, type){
  print(gene)
  # Draw the heatmap for the Isoform expression of the gene 
  
  if(type == "targeted"){
    # retained isoforms from merged annotations
    retained= c(Merged_noISM[Merged_noISM$associated_gene == gene,"IsoSeq_isoform"],
                Merged_noISM[Merged_noISM$associated_gene == gene,"ONT_isoform"])
    
    Norm_transcounts = Norm_transcounts %>% filter(isoform %in% retained)
  }
  
  # Subset the normalised expression count to gene, and datawrangle for plot
  dat = Norm_transcounts[Norm_transcounts$associated_gene == gene,c("isoform","value","variable")] %>%
    mutate(value = log2(value)) %>%
    spread(., isoform, value) %>% tibble::column_to_rownames(var = "variable") 
  
  # remove isoforms that have been removed by tappAS due to very low count 
  # replace infinity value from log value with 0 
  # rotate the dataframe for visualisation ease
  dat <- dat[,colSums(is.na(dat))<nrow(dat)]
  dat[dat == "-Inf"] <- 0
  dat.t <- t(dat)
  
  # set the order for the column (Age, Genotype)
  coldata = Norm_transcounts %>% 
    dplyr::select(sample, time, group) %>% distinct(.keep_all = TRUE) %>% column_to_rownames(var = "sample") %>% 
    mutate(time = as.factor(time))
  colnames(coldata) = c("Age (months)","Genotype")

  
  # set the order for the row (isoform structural category)
  rowdata = cf[cf$isoform %in% colnames(dat),c("isoform","structural_category")] %>% dplyr::select(-isoform)
  colnames(rowdata) = c("Category")
  
  # set annotation colours
  if(type == "targeted"){
    annotation_colors = list(
      Genotype = c(WT=wes_palette("Royal1")[1], TG=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white","4"="#CFCFCF","6"="#777777","8"="black"), 
      Category = c(FSM = alpha("#00BFC4",0.8),ISM = alpha("#00BFC4",0.3),NIC = alpha("#F8766D",0.8),NNC = alpha("#F8766D",0.3)))
  }else{
    annotation_colors = list(
      Genotype = c(WT=wes_palette("Royal1")[1], TG=wes_palette("Royal1")[2]),
      `Age (months)` = c("2"="white", "8"="black"), 
      Category = c(FSM = alpha("#00BFC4",0.8),ISM = alpha("#00BFC4",0.3),NIC = alpha("#F8766D",0.8),NNC = alpha("#F8766D",0.3)))
  }

  
  # draw heatmap
  if(nrow(dat.t) > 1){
    p = pheatmap(dat.t, annotation_col=coldata, annotation_row = rowdata, annotation_legend = FALSE,
                 show_colnames = FALSE,show_rownames = FALSE, color = viridis(10),annotation_colors = annotation_colors,
                 fontsize_col = 20)
    }else{
    p = ggplot()
  }
  
  return(p)
}

#InputGene = "Abca7"
#Norm_transcounts = IsoExp$Norm_transcounts
plot_transexp <- function(InputGene,Norm_transcounts,type, name){
  plot_title <- paste0(InputGene,"\n","Transcript Expression","\n")
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene)
  df$grouping1 <- paste0(df$group, "_", df$time)
  df <- df[order(match(df$grouping1, c("WT_2","WT_8","TG_2","TG_8"))),]
  #df <- df %>% arrange(desc(group)) %>% arrange(desc(time)) 
  df$Age <- as.factor(df$time)
  # df$group <- factor(df$group, level = c("WT","TG"))
  #df$time <- factor(paste(df$time,"months"), levels = c("2 months","8 months"))
  df$groupings <- paste0(df$isoform,"_", df$group, "_", df$time)
  

  # retained isoforms from merged annotations
  retained= c(Merged_noISM[Merged_noISM$associated_gene == InputGene,"IsoSeq_isoform"],
              Merged_noISM[Merged_noISM$associated_gene == InputGene,"ONT_isoform"])
  
  df = df %>% filter(isoform %in% retained)
  top_mean_isoforms = aggregate(df$value, list(df$isoform), sum) %>% arrange(desc(x)) %>% .[c(1:3),"Group.1"]
  
  p <- ggplot(df[df$isoform %in% top_mean_isoforms,],aes(x = Age, y = value, colour = group)) + 
    #ggplot(df, aes(x = Age, y = value, colour = group)) + 
    geom_point(size = 3) +
    #geom_line(aes(group = groupings, colour = group), position = position_dodge(width = 0.3)) + 
    #geom_point(size = 4, aes(group = groupings, shape = Age, colour = group), position = position_dodge(width = 0.3)) +
    mytheme + labs(x = "Age (months)", y = "Iso-Seq Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines")) +
    legend_theme +
    stat_summary(data=df[df$isoform %in% top_mean_isoforms,], aes(x=Age, y=value, group=group), fun ="mean", geom="line", linetype = "dotted") +
    scale_y_continuous(trans='log10') + facet_grid(~ isoform) + 
    theme(legend.position = "none") +
    scale_fill_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT")))
  #+ facet_grid(~ Age, labeller=as_labeller(c(`2` = "2 mos", `4` = "4 mos", `6` = "6 mos",`8` = "8 mos"))) +
    #guides(shape = guide_legend(order = 2),col = guide_legend(order = 1)) + 
  
  if(type == "isoseq"){
    p <- p + #scale_shape_manual(name = "Age (months)", values=c(1, 17), label = c("2 months", "8 months")) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  } else {
    p <- p + #scale_shape_manual(name = "Age (months)", values=c(1, 16, 2, 17)) + 
      scale_colour_manual(name = "Genotype", values = c(label_colour("TG"),label_colour("WT"))) 
  }
  
  return(p)
}



plot_transexp_overtime <- function(InputGene,Norm_transcounts,type, name){
  print(InputGene)
  
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(InputGene,"\n",name,"\n\n")

  if(InputGene %in% c("Gfap","C4b","Ctsd","Gatm","H2-D1","Padi2","Cd34","Ubqln1")){
    if(type == "isoseq"){
      sigiso = c(tappassigtrans$WholeIso_Transexp[tappassigtrans$WholeIso_Transexp$associated_gene == InputGene, "isoform"])$isoform
      df = df %>% mutate(sig = ifelse(isoform %in% sigiso,Isoform,"NA")) %>% filter(sig != "NA")
    }else{
      sigisoiso = c(tappassigtrans$WholeIso_Transexp[tappassigtrans$WholeIso_Transexp$associated_gene == InputGene, "isoform"])$isoform
      df = df %>% mutate(Isoform = ifelse(isoform %in% sigisoiso,Isoform,"RNA-specific"))
    }
  }
  
  if(InputGene == "Gfap"){
    colours = c("#F8766D","#7CAE00","#00BFC4","#C77CFF",alpha("grey",0.6))
  } else if(InputGene == "Ubqln1"){
    colours = c("#F8766D","#00BFC4",alpha("grey",0.6))
  }else{
    colours = c("#F8766D",alpha("grey",0.6))
  }
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point(size = 3) + 
    facet_grid(~group,scales = "free", space = "free") +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted", size = 1.5) +
    mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),
          panel.spacing = unit(2, "lines"),
          legend.position = c(0.4,0.8)) 
  
  if(type == "isoseq"){
    p <- p + scale_colour_discrete(name = "Isoform")
  }

  if(type != "isoseq"){
    p <- p + theme(legend.position = "none") + scale_colour_manual(values = colours)
  }
  

  return(p)
}

#InputGene = "Trem2"
#Norm_transcounts = OntExp$Norm_transcounts
#difftrans = tappassigtrans$ont$TargetedOnt_Transexp
#name = "ONT_Expression"
#plot_general_transexp_overtime("App",IsoExp$Norm_transcounts,tappassigtrans$iso$TargetedIso_Transexp,"IsoSeq_Expression")
plot_general_transexp_overtime <- function(InputGene,Norm_transcounts,difftrans, name){
  
  difftrans = difftrans %>% filter(associated_gene == InputGene) %>% arrange(`p-value`)
  # top 5 differentially-ranked isoforms 
  df <- Norm_transcounts  %>% filter(associated_gene == InputGene) %>% filter(isoform %in% difftrans$isoform[1:3])
  plot_title <- paste0(InputGene,"\n",name,"\n\n")
  if(nrow(df) != 0){
    p <- ggplot(df, aes(x = time, y = value, colour = isoform)) + geom_point(size = 3) + 
      facet_grid(~group,scales = "free", space = "free") +
      stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
      mytheme + labs(x = "Age (months)", y = "Normalised Counts",title = plot_title) +
      theme(strip.background = element_blank(), legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines"),
            axis.text.y = element_text(angle = 90, hjust = 0.5)) 
  }else{
    p <- ggplot() + theme_void()
  }
  
  if(name == "Iso-Seq Expression"){
    p = p + scale_colour_manual(values=plot_target_diff_col(InputGene,"IsoSeq")) + 
      labs(y = "Iso-Seq Normalised Counts",title = paste0(InputGene,"\n","Transcript Expression","\n"))
    }else{
    p = p + scale_colour_manual(values=plot_target_diff_col(InputGene,"ONT")) + 
      labs(y = "ONT Normalised Counts", title = paste0(InputGene,"\n","Transcript Expression","\n"))
  }
  
  return(p)
}

DIU_analysis_output_venn <- function(){
  
  # whole trancriptome: Iso-Seq scaffold + Iso-Seq expression 
  # venn diagram of the genes identified as differentially isoform usage 
  # methods for removing lowly expressed isoforms:
  # proportion: threshold at 0.2 (20%)
  # fold change: threshold at 2.5 (relative fold change)
  v1 = twovenndiagrams(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_isoseq_fc$gene,"DIU Genes: Iso-Seq Expression \n filtered by proportion","DIU Genes: Iso-Seq Expression \n filtered by fold change")
  v2 = twovenndiagrams(tappasDIU$DIU_rnaseq_prop$gene,tappasDIU$DIU_rnaseq_fc$gene,"DIU Genes: RNA-Seq Expression \n filtered by proportion","DIU Genes: RNA-Seq Expression \n filtered by fold change")
  v3 = twovenndiagrams(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_rnaseq_prop$gene,"DIU Genes: Iso-Seq Expression \n filtered by proportion","DIU Genes: RNA-Seq Expression \n filtered by proportion")
  v4 = twovenndiagrams(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene,"DIU Genes: Iso-Seq Expression \n filtered by fold change","DIU Genes: RNA-Seq Expression \n filtered by fold change")
  
  cat("Number of Isoforms with DIU - Isoseq Fold Change:", length(tappasDIU$DIU_isoseq_fc$gene),"\n")
  cat("Number supported by RNA-Seq:", length(intersect(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene)),"\n")
  
  output <- list(v1,v2,v3,v4)
  names(output) <- c("v1","v2","v3","v4")
  return(output)
}



### Differential Isoform Usage #################################################################
Gene_overall_exp <- function(gene_transcripts_input, normexp_input, type){
  genetrans = gene_transcripts_input
  normexp = normexp_input
  
  if(type == "Iso_Targeted"){  
    phenotype$iso = phenotype$iso %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                             col = paste0(sample,"_",pheno_group,"_",time))
    names(normexp) <- phenotype$iso$col[match(names(normexp), phenotype$iso$sample)]
    
  }else{
    phenotype$ont = phenotype$ont %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                             col = paste0(sample,"_",pheno_group,"_",time))
    #names(normexp) <- phenotype$ont$col[match(names(normexp), phenotype$ont$sample)]
  }
  
  # data wrangle for the gene name and isoform type
  dat = as.data.frame(normexp) %>% rownames_to_column(., var = "isoform") %>% 
    merge(., genetrans[,c("transcript","geneName")], by.x = "isoform", by.y = "transcript") %>% 
    mutate(isoform_type = ifelse(grepl("ENS", isoform), "Known","Novel"))
  
  # function to parse through each associated gene to find the mean across the WT and TG (merged dataset) for IF
  parse_genotype_IF <- function(gene){
    d = dat[dat$geneName == gene,] %>% select(-c("isoform","geneName","isoform_type")) %>% apply(.,1,mean) %>%
      as.data.frame() %>% 
      apply(.,2, function(x) x/sum(x) * 100) %>%
      cbind(., dat[dat$geneName == gene,c("isoform","geneName","isoform_type")]) %>%
      `colnames<-`(c("perc", "isoform", "associated_gene","isoform_type"))
    
    return(d)
  }
  
  gIF_lst = lapply(TargetGene, function(x) parse_genotype_IF(x))
  gIF_lst = do.call(rbind, gIF_lst)
  
  p = gIF_lst %>% mutate(isoform_col = ifelse(perc < 5,"Other",isoform)) %>%
    ggplot(., aes(x = isoform_col, y = perc,fill = isoform_type)) + geom_bar(stat = "identity") + 
    facet_grid(~associated_gene, scales='free') + mytheme + 
    theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom",
          strip.background = element_blank(), strip.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20), legend.text=element_text(size=20)) + 
    labs(y = "Isoform Fraction (%)", x = "Isoform") + 
    scale_fill_manual(name = "Isoform Classification", values = c(label_colour("known"),label_colour("novel"))) 

  #ggplot(dat, aes(x = associated_gene, y = perc, fill = isoform_type)) + geom_bar(stat = "identity")
  return(p)
}

#gene = "Trem2"
#genetrans = tappasont$gene_transcripts.tsv
#normexp = tappasont$input_norm
#type = "ONT_Targeted"
IF_plot <- function(gene, gene_transcripts_input, normexp_input, type){
  
  print(gene)
  genetrans = gene_transcripts_input
  normexp = normexp_input
  
  # subset expression by all the isoforms associated with the gene 
  iso = subset(genetrans, geneName == gene) 
  isoexp = normexp[which(rownames(normexp) %in% iso$transcript),]
  
  if(type == "Iso_Targeted"){  
    phenotype$iso = phenotype$iso %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                             col = paste0(sample,"_",pheno_group,"_",time))
    names(isoexp) <- phenotype$iso$col[match(names(isoexp), phenotype$iso$sample)]
    
    # perform prefiltering by fold change 
    #normexp = DIU_time_analysis_filteronly(tappas_input_dir$iso,"3")
    
    
  }else if(type == "ONT_Targeted"){
    phenotype$ont = phenotype$ont %>% mutate(pheno_group = ifelse(group == "CONTROL","WT","TG"),
                                             col = paste0(sample,"_",pheno_group,"_",time))
    names(isoexp) <- phenotype$ont$col[match(names(isoexp), phenotype$ont$sample)]
    
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

# plot expression of DIU genes
# DIU genes identified from using Iso-Seq expression as abundance, and commonly identified from both filtering methods (proportion and fold change)
# Gene expression = sum of normalised expression or FL read counts from associated isoforms that were prefiltered 
# prefiltered isoforms = filtered lowly-expressed isoforms as part of normalisation (refer to prefiltering prior to TMM)
# mean gene expression/median gene expression = mean across all the samples 
# 1. Determined normalised gene expression/read counts across each sample by summing normalised isoform expression/read counts 
# 2. Determined mean or median of normalised gene expression across all samples (n = 12)
DIU_genes_exp <- function(){
  # common genes identified by both proporition and fold change
  # commonDIU_isoiso: whole transcriptome: Iso-Seq scaffold + RNA_Seq expresion 
  commonDIUisoiso = intersect(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_isoseq_fc$gene)
  
  #commonDIUisorna_prop = intersect(tappasDIU$DIU_isoseq_prop$gene,tappasDIU$DIU_rnaseq_prop$gene)
  #commonDIUisorna_fc = intersect(tappasDIU$DIU_isoseq_fc$gene,tappasDIU$DIU_rnaseq_fc$gene)
  #commonDIUisorna = intersect(commonDIUisorna_fc,commonDIUisorna_prop)
  
  ###1. Using classification files (whole transcriptome), obtain the total gene FL counts associated with the isoforms that were retained from pre-filtering 
  # normalised_matrix.tsv = isoforms that were retained after prefiltering
  # subset the classification files on retained isoforms, and obtain total gene FL counts (across all samples)
  # genes_iso = number of prefiltered isoforms associated per gene
  retained_transcripts = row.names(tappasiso$input_norm)
  FL_genes = class.files %>% filter(isoform %in% retained_transcripts) %>% group_by(associated_gene) %>% tally(FL) 
  genes_iso = class.files %>% filter(isoform %in% retained_transcripts) %>% group_by(associated_gene) %>% tally()
  
  ###2. Obtain the median and mean of gene FL counts across the samples 
  # sum_dat = sum of the FL counts per gene per sample (FL counts of retained isoforms)
  # med_sumdat = median of the gene across all samples 
  # mean_sumdat mean of the gene across all samples
  dat = class.files %>% filter(isoform %in% retained_transcripts) %>% dplyr::select(associated_gene,isoform, contains("FL.")) 
  rownames(dat) <- paste(dat$isoform,"_",dat$associated_gene)
  sum_dat = dat %>% dplyr::select(-isoform) %>% group_by(associated_gene) %>% summarise(across(everything(), ~ sum(., is.na(.), 0))) %>% column_to_rownames(., var = "associated_gene")
  med_sumdat = data.frame(apply(sum_dat,1,median)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "median_FL_expression"))
  mean_sumdat = data.frame(apply(sum_dat,1,mean)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "mean_FL_expression"))
  
  # Merge all the counts for DIU genes 
  ###3. Use only the commonly identified genes from lowly-expressed filtering methods: proportion and fold change 
  merged_counts = merge(merge(merge(tappasDIU$DIU_isoseq_prop, FL_genes, by.x = "gene", by.y = "associated_gene"),mean_sumdat, by = "gene"),med_sumdat, by = "gene") %>% filter(gene %in% commonDIUisoiso)
  
  ## check normalised gene expression make sense using example of 1110032F04Rik
  ## iso = isoforms associated with 1110032F04Rik, note can not just use pb.id as wrong assignment 
  ## iso_exp = obtain the normalised expression of those isoforms 
  ## sum the iso_exp per column across the samples --> gene expression --> determine median and mean across samples
  ## median and mean == merged_counts$mean_expression, and merged_counts$median_expression
  #iso = tappasiso$result_gene_trans.tsv[tappasiso$result_gene_trans.tsv$gene == "1110032F04Rik","isoform"]
  #iso_exp = tappasiso$input_norm %>% rownames_to_column(var = "isoform") %>% filter(isoform %in% iso) 
  #median(apply(iso_exp[,-1],2,sum))
  #mean(apply(iso_exp[,-1],2,sum))
  #head(merged_counts)
  
  # plots of gene expression of commonly identified DTU genes, with the normalised gene expression, and FL reads
  p1 = density_plot(merged_counts,"sum_expression","n","Normalised Gene Counts","Gene FL Reads","") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
  p2 = density_plot(merged_counts,"median_expression","median_FL_expression","Median Normalised Gene Counts","Median Gene FL Reads","") + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10")
  p3 = density_plot(merged_counts,"mean_expression","mean_FL_expression","Mean Normalised Gene Counts","Mean Gene FL Reads","") + scale_y_continuous(trans = "log10", labels = scales::number_format(accuracy = 1)) + scale_x_continuous(trans = "log10",labels = scales::number_format(accuracy = 1))
  
  output = list(p1,p2,p3)
  names(output) = c("p1","p2","p3")
  return(output)
}

# Using Fold Change from RNASeq differential isoform usage
# Podium Change = TRUE or FALSE for whether differential isoform usage
DIU_RNASEQ_results <- function(){
  cat("Number of Genes with DIU:", nrow(tappasDIU$rnaseq %>% filter(adjPValue < 0.05)),"\n")
  cat("Number of Genes with DIU with no major switching:", nrow(tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "NO")),"\n")
  cat("Number of Genes with DIU with no switching:", nrow(tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "YES")),"\n")
  
  # significant DEG with interaction effects using Iso-Seq as expression
  genesigs_interaction_wholeiso = c(gene_sigs_WholeIso_lst$models$`Model 4 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 5 Interaction`$...1,
                                    gene_sigs_WholeIso_lst$models$`Model 6 Interaction`$...1,gene_sigs_WholeIso_lst$models$`Model 7 Interaction`$...1)
  
  # significant DEG with interaction effects using RNA-Seq as expression
  genesigs_interaction_wholerna = c(gene_sigs_WholeRNA_lst$models$`Model 4 - 7 Interaction`$...1)
  
  # Cases of Differentially Expressed Genes with Differential Isoform Usage
  DEG_DIU = list(
    #case1 = DIU genes with major switching isoform, and DEG
    #case2 = DIU genes but no major switching, and DEG
    #case3 = DIU genes with no major switching, but not DEG
    #case4 = DIU genes with major switching, but not DEG
    tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "YES" & gene %in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "NO" & gene %in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "NO" & gene %!in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue)),
    tappasDIU$rnaseq %>% filter(adjPValue < 0.05 & podiumChange == "YES" & gene %!in% c(genesigs_interaction_wholerna)) %>% arrange(-desc(adjPValue))
  )
  names(DEG_DIU) = c("DIU_DEG_maj","DIU_DEG_nomaj","DIU_notDEG_nomaj","DIU_notDEG_maj")
  
  cat("Number of Genes with DIU with major switching, but genes are differentially expressed:", nrow(DEG_DIU$DIU_DEG_maj),"\n")
  cat("Number of Genes with DIU with no major switching, but genes are differentially expressed:", nrow(DEG_DIU$DIU_DEG_nomaj),"\n")
  cat("Number of Genes with DIU with no major switching, but genes are not differentially expressed:", nrow(DEG_DIU$DIU_notDEG_nomaj),"\n")
  cat("Number of Genes with DIU with major switching, but genes are not differentially expressed:", nrow(DEG_DIU$DIU_notDEG_maj),"\n")
  
  return(DEG_DIU)
}


diff_across_rnaseq <- function(gene){
  colours = c(wes_palette("Darjeeling2")[[2]],wes_palette("Darjeeling2")[[3]],
              alpha(wes_palette("Cavalcanti1")[[2]],0.5),
              wes_palette("Darjeeling2")[[5]],
              wes_palette("Darjeeling2")[[4]])
  p <- list(plot_mergedexp(gene,"NA",RNAExp$GeneExp,RNAExp$Norm_transcounts,"RNA-Seq Gene Expression"),
            plot_transexp_overtime(gene,RNAExp$Norm_transcounts,"isoseq","RNA-Seq Isoform Expression") + 
              theme(legend.position = "left", legend.direction = "vertical") +
              scale_colour_manual(values = colours),
            IF_plot(gene,tappasrna$gene_transcripts.tsv, tappasrna$input_norm, "rnaseq")[[1]] +
              theme(legend.position = "left"),
            IF_plot(gene,tappasrna$gene_transcripts.tsv, tappasrna$input_norm, "rnaseq")[[3]] +
              theme(legend.position = "none") + scale_fill_manual(values = colours))
  return(p)
}

IR_ORF <- function(){
  dat = lapply(group.class.files, function(x) x %>% filter(subcategory == "intron_retention") %>% group_by(associated_gene) %>% tally())
  IR = bind_rows(dat, .id = "Sample") %>% spread(., Sample, n) %>% mutate(TG_WT_8mos = `TG_8mos` - `WT_8mos`,TG_WT_2mos = `TG_2mos` - `WT_2mos`)
  
  total_genes = sapply(group.class.files, function(x) length(unique(x$associated_gene))) %>% reshape2::melt(value.name = "total_genes") %>% rownames_to_column(var = "Sample")
  total_iso = sapply(group.class.files, function(x) nrow(x)) %>% reshape2::melt(value.name = "total_isoforms") %>% rownames_to_column(var = "Sample")
  IR_genes = bind_rows(dat, .id = "Sample") %>% group_by(Sample) %>% tally() %>% full_join(.,total_genes, by = "Sample") %>% mutate(perc = n/total_genes * 100)
  
  
  plot_Genes <- function(dataset,grp, age_name, y_name){
    x.var <- rlang::sym(quo_name(enquo(grp)))
    cols <- enquo(grp) 
    
    
    # difference in number of isoforms tally 
    dat2 = dataset %>% group_by_at(vars(!!cols)) %>% tally()
    dat2$sign_freq = dat2$n * ifelse(dat2[[grp]] > 0 , 1, -1)
    dat2[[grp]]  = reorder(dat2[[grp]] , dat2$sign_freq, sum)
    
    if(grp == "TG_WT_8mos"){
      dat2_agg = aggregate(sign_freq ~ TG_WT_8mos, FUN = sum, data = dat2)
    } else{
      dat2_agg = aggregate(sign_freq ~ TG_WT_2mos, FUN = sum, data = dat2) 
    }
    dat2_agg[[grp]] = as.numeric(as.character(dat2_agg[[grp]]))
    dat2_agg$col = ifelse(dat2_agg[[grp]] == "0", "neutral", ifelse(dat2_agg$sign_freq < 0, "negative","positive"))
    
    print(dat2_agg)
    p <- ggplot(dat2_agg, aes(x = !! x.var , y = abs(sign_freq), fill = factor(col))) +
      geom_bar(position = 'identity', stat = "identity") +
      guides(fill = FALSE) + mytheme +
      labs(y = y_name, x = paste0("Difference in number of IR isoforms per gene between WT and TG \n at ", age_name,"(TG - WT)"))
    
    return(p)
  }
  
  # number of ORFs
  ORF = lapply(group.class.files, function(x) x %>% filter(predicted_NMD == "TRUE") %>% group_by(associated_gene) %>% tally())
  ORF_genes = bind_rows(ORF, .id = "Sample") %>% group_by(Sample) %>% tally() %>% full_join(.,total_genes, by = "Sample") %>% mutate(perc = n/total_genes * 100)
  ORF_transcripts = bind_rows(ORF, .id = "Sample") %>% group_by(Sample) %>% tally(n) %>% full_join(.,total_iso, by = "Sample") %>% mutate(perc = n/total_isoforms * 100)
  
  ORF_counts = bind_rows(ORF, .id = "Sample") %>% spread(., Sample, n) %>% mutate(TG_WT_8mos = `TG_8mos` - `WT_8mos`,TG_WT_2mos = `TG_2mos` - `WT_2mos`)
  
  p1 = plot_Genes(IR,"TG_WT_8mos","8 months","Number of Genes with intron-retained isoforms")
  p2 = plot_Genes(IR,"TG_WT_2mos","2 months","Number of Genes with intron-retained isoforms")
  p3 = plot_Genes(ORF_counts,"TG_WT_8mos","8 months","Number of Genes with isoforms predicted for ORF")
  p4 = plot_Genes(ORF_counts,"TG_WT_2mos","2 months","Number of Genes with isoforms predicted for ORF")
  
  return(list(p1,p2,p3,p4))
}

### Functional Diversity Plots #################################################################
FDA_plot <- function(FDA_input, title){
  
  TranscriptFeatures = c("NMD","repeat","3UTR Motif","5UTR Motif","PAS","uORF","miRNA Binding","3UTR Length","5UTR Length","CDS","PolyA Site","3UTR Motif","5UTR Motif")
  
  # data wrangle for plot
  dat = reshape2::melt(FDA_input, id = "Gene", variable = "Category") %>% group_by(Category,value) %>% tally() %>% filter(value != "") %>% 
    mutate(type = factor(ifelse(Category %in% TranscriptFeatures, "Transcript","Protein"), levels = c("Transcript","Protein"))) %>%
    mutate(perc = n/nrow(FDA_FP) * 100) 
  
  tab = dat %>% arrange(desc(perc))
  
  p = ggplot(dat, aes(x = reorder(Category,n), y = n, fill = type, label = n)) + geom_bar(aes(alpha = value),stat = "identity") + coord_flip() + 
    scale_fill_manual(values = alpha(c("red", "blue"), .3)) + 
    facet_wrap(~type,nrow = 2, scales = "free_y") + mytheme + labs(y = "", x = "", title = title) + 
    guides(fill="none",alpha=guide_legend(title="Status")) + theme(legend.position = c(0.8,0.1),strip.background = element_blank()) 
  
  return(list(tab,p))
}

# Plot venn diagram of the genes varying by 5'UTR, 3'UTR, CDS, PolyA
FDA_lengths_plot <- function(){
  
  x <- list(
    `3'UTR` = FDA_lengths$`3UTR`[FDA_lengths$`3UTR`$Status == "YES","Gene"], 
    `5'UTR` = FDA_lengths$`5UTR`[FDA_lengths$`5UTR`$Status == "YES","Gene"] , 
    CDS = FDA_lengths$`CDS`[FDA_lengths$`CDS`$Status == "YES","Gene"],
    PolyA = FDA_lengths$PolyASite[FDA_lengths$PolyASite$Status == "YES","Gene"]
  )
  
  # Example of Genes with varying 3'UTR, 5'UTR, CDS, PolyA site
  intersect(intersect(intersect(intersect(FDA_lengths$`3UTR`[FDA_lengths$`3UTR`$Status == "YES","Gene"],
                                          FDA_lengths$`5UTR`[FDA_lengths$`5UTR`$Status == "YES","Gene"]),
                                FDA_lengths$`CDS`[FDA_lengths$`CDS`$Status == "YES","Gene"]),
                      FDA_lengths$PolyASite[FDA_lengths$PolyASite$Status == "YES","Gene"]),
            FDA_lengths$mir3225[FDA_lengths$mir3225$Status == "YES","Gene"]) 
  
  p = ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 0.5, set_name_size = 4)
  return(p)
}

DFI_generateplots <- function(){
  cat("Number of DFI:", DFI %>% filter(DFI == "DFI") %>% nrow(),"\n")
  cat("Number of DFI:", DFI %>% filter(DFI == "DFI") %>% select(Gene) %>% unique() %>% nrow(),"\n")
  
  lstTrans = c("3UTRmotif","5UTRmotif","RNA_binding","uORF","miRNA","miRNA_Binding","PAS","repeat")
  
  dfi_results = data.frame()
  for(feature in unique(DFI$Feature)){
    newRow = data.frame(Feature=feature,
                        TotalCount=nrow(DFI[DFI$Feature==feature,]),
                        DFICount=nrow(DFI[DFI$Feature==feature & DFI$DFI == "DFI",]),
                        TotalAnnot=sum(DFI_counts[DFI_counts$Feature==feature,]$Total),
                        Level=ifelse(feature %in% lstTrans, "Transcript","Protein"), 
                        Time = "0")
    dfi_results = rbind(dfi_results,newRow)
  }
  
  dfi_results["Freq1"] = dfi_results[,"TotalCount"]/sum(dfi_results[,"TotalCount"])*100
  dfi_results["Freq2"] = dfi_results[,"TotalAnnot"]/sum(dfi_results[,"TotalAnnot"])*100
  dfi_results$Feature = as.character(dfi_results$Feature)
  dfi_results = dfi_results[order(dfi_results$Feature, method = ),]
  
  # original plot
  p1 = ggplot(dfi_results) +
    geom_bar(aes(x=Feature,y=Freq2,fill=as.factor(Time)),width =0.8, size=0.5,
             color="black", stat = "identity", position = position_dodge()) +
    geom_line(aes(x=Feature,y=Freq1), stat="identity", group = 1) +
    geom_point(aes(x=Feature,y=Freq1, shape="Features"), stat="identity", group = 1) +
    ylab("% Features") +
    xlab("Category") +
    labs(fill= "", shape=NULL) +
    scale_fill_manual(labels = c("DFI Features"), values = c(myPalette[c(2)])) +
    theme_classic() +
    theme(axis.title.x = element_text(size=17, margin=margin(5,0,0,0)),
          axis.text.x  = element_text(margin=margin(7,0,0,0), size=17, hjust = 1, angle = 45),
          axis.title.y = element_text(size=17,  margin=margin(0,15,0,0)),
          axis.text.y  = element_text(vjust=0.5, size=17))  +
    facet_grid(~ Level, scales = "free",  space="free") +  theme(strip.text.x = element_text(size = 15, face="bold"), strip.background = element_rect(colour="black", fill=c("white")))
  
  # original plot reformatted
  p1b = ggplot(dfi_results) +
    geom_bar(aes(x=Feature,y=Freq2,fill=as.factor(Time)),width =0.8, size=0.5,
             color="black", stat = "identity", position = position_dodge()) +
    geom_line(aes(x=Feature,y=Freq1), stat="identity", group = 1) +
    geom_point(aes(x=Feature,y=Freq1, shape="Features"), stat="identity", group = 1) +
    ylab("Features (%)") +
    xlab("Category") +
    labs(fill= "", shape=NULL) +
    scale_fill_manual(labels = c("DFI Features"), values = c(myPalette[c(2)])) +
    mytheme +  facet_grid(~ Level, scales = "free",  space="free") +  
    theme(strip.text.x = element_text(size = 15, face="bold"), strip.background = element_rect(colour="white", fill=c("white"))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "bottom") 
  
  
  # Plot 2
  dfFavored = DFI[DFI$DFI == "DFI",c("Feature","Favored")]
  total_favored = data.frame()
  
  for(feature in unique(dfFavored$Feature)){
    group = unique(dfFavored$Favored)
    for(name in group){
      newRow = data.frame(
        "Feature" = feature,
        "Favored" = name,
        "Count" = nrow(dfFavored[dfFavored$Feature==feature & dfFavored$Favored==name,]),
        "Level" = ifelse(feature %in% lstTrans, "Transcript","Protein")
      )
      total_favored = rbind(total_favored,newRow)
    }
    total_favored[total_favored$Feature==feature,]$Count <- (total_favored[total_favored$Feature==feature,]$Count/
                                                               sum(total_favored[total_favored$Feature==feature,]$Count))*100
  }
  
  if(length(which(total_favored$Favored=="N/A"))>0)total_favored = total_favored[-which(total_favored$Favored=="N/A"),]
  total_favored = total_favored[order(total_favored$Level, method = ),] 
  total_favored$Level = factor(total_favored$Level, levels = c("Transcript","Protein"))
  
  p2 = ggplot(total_favored, aes(x=Feature, y = Count, fill=Favored)) +
    geom_bar(stat = "identity", position = "fill") +
    mytheme + labs(x = "Category", y = "DFI Features") +
    facet_grid(~ Level, scales = "free",  space="free") +
    theme(strip.text.x = element_text(size = 15, face="bold"), strip.background = element_rect(colour="white", fill=c("white"))) +
    geom_hline(yintercept=0.5, linetype="dashed",
               color = myPalette[1], size=1)+
    scale_fill_manual(values = c(label_colour("TG"),label_colour("WT")), labels = c("WT","TG")) +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "bottom") 
  
  # Co inclusion 
  simplify_features <- function(feature){
    simplify = if(grepl("miR", feature) == "TRUE"){"miRNA"}else if(grepl("U00", feature) == "TRUE"){"uORF"}else(as.character(feature))
    return(simplify)
  }
  
  #DFI_co$FeatureID1_simplified = lapply(DFI_co$FeatureID1, function(x) simplify_features(x))
  #DFI_co$FeatureID2_simplified = lapply(DFI_co$FeatureID2, function(x) simplify_features(x))
  #DFI_co$Pair = paste0(DFI_co$FeatureID1_simplified,"_",DFI_co$FeatureID2_simplified)
  #DFI_co = DFI_co %>% filter(GeneswithBothDFI > 10)
  
  #merge(DFI_co %>% group_by(Pair) %>% tally(MutualExclusive) %>% dplyr::rename("MutualExclusive" = "n"),
  #      DFI_co %>% group_by(Pair) %>% tally(Coinclusion) %>% dplyr::rename("Coinclusion" = "n")) %>% reshape2::melt() %>% 
  #  ggplot(., aes(x = Pair, y = value, fill = variable)) + geom_bar(stat = "identity") +
  #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "bottom")
  
  
  return(list(p1,p1b,p2))
}


DEI_genes_unique <- function(x){
  # differentially expressed isoforms across interactions detected in Iso-Seq dataset
  DEI_IsoSeq = unique(c(trans_sigs_WholeIso_lst$models$`Model 4 Interaction`$isoform,
                        trans_sigs_WholeIso_lst$models$`Model 5 Interaction`$isoform,
                        trans_sigs_WholeIso_lst$models$`Model 6 Interaction`$isoform,
                        trans_sigs_WholeIso_lst$models$`Model 7 Interaction`$isoform))
  
  # differentially expressed isoforms unique in Iso-Seq dataset
  DEI_IsoSeq_only = setdiff(DEI_IsoSeq,unique(trans_sig_WholeRNA_lst$models$`Model 4 - 7 Interaction`$isoform))
  DEI_IsoSeq_RNAseq = intersect(DEI_IsoSeq,unique(trans_sig_WholeRNA_lst$models$`Model 4 - 7 Interaction`$isoform))
  cat("Number of Differentially expressed isoforms (Iso-Seq) associated with genotype and interaction effects:", 
      length(DEI_IsoSeq),"\n")
  cat("Number of Differentially expressed isoforms (Iso-Seq) not supported by RNA-Seq expression input:",
      length(DEI_IsoSeq_only),"\n")
  
  DEI_Isoseq_only_FL = class.files[class.files$isoform %in% DEI_IsoSeq_only,] %>% dplyr::select(starts_with("FL.")) %>% apply(.,1,mean) %>% reshape2::melt() %>% merge(.,class.files[,c("isoform","associated_gene")],by.x = 0, by.y = "isoform", all.x = T) 
  
  DEI_IsoSeq_RNAseq_FL = class.files[class.files$isoform %in% DEI_IsoSeq_RNAseq,] %>% dplyr::select(starts_with("FL.")) %>% apply(.,1,mean) %>% reshape2::melt() %>% merge(.,class.files[,c("isoform","associated_gene")],by.x = 0, by.y = "isoform", all.x = T) 
  cat("Number of DIE (IsoSeq) with median FL reads < 24):", nrow(DEI_Isoseq_only_FL %>% filter(value < 24)))
}

