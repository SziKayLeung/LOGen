# plot gene expression of all samples as box plot by phenotype
plot_gene_exp <- function(InputGene, GeneExp, Norm_transcounts, type, name){
  
  df <- GeneExp %>% filter(associated_gene == InputGene) 
  plot_title <- paste0(name," ",InputGene)
  
  groups = levels(Norm_transcounts$group)
  if("Control" %in% groups){
    df <- df %>% mutate(group = factor(group, levels = c("Control",groups[!groups %in% "Control"])))
  }else{
    df <- df %>% mutate(group = factor(group, levels = c("CONTROL",groups[!groups %in% "CONTROL"]))) 
  }
  
  if(type == "case_control"){
    p <- ggplot(df, aes(x = group, y = Exp, colour = group)) + geom_boxplot() + 
      labs(y = "Normalised counts", x = "", title = paste0(plot_title)) 
  }else{
    # time series
    p <- ggplot(df, aes(x = time, y = Exp, colour = group)) + geom_point(size = 3) + 
      stat_summary(data=df, aes(x=time, y=Exp, group=group), fun="mean", geom="line", linetype = "dotted") + 
      labs(y = "Normalised counts", x = "Age (months)", title = paste0(plot_title,"\n\n")) 
  }
  
  p <- p + mytheme +
    # colours in the opposite direction of the levels 
    scale_colour_manual(values = c(label_colour("Case"),label_colour("Control")), name = "Phenotype",
                        c(label_group("Control"),label_group("Case"))) + 
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  
  
  return(p)
}


plot_raw_expression_bygene <- function(FL_reads, gene){
  iso <- class.files[class.files$associated_gene == gene,"isoform"]
  
  p <- FL_reads[rownames(FL_reads) %in% iso,] %>% rownames_to_column(var = "isoform") %>%
    reshape2::melt() %>% 
    dplyr::rename(sample = variable) %>%
    full_join(., phenotype$ont, by = "sample") %>% 
    left_join(., class.files[,c("isoform","structural_category")], by = "isoform") %>%
    ggplot(., aes(x = isoform, y = value)) + geom_boxplot(aes(colour = group)) + 
    mytheme + labs(x = " ", y = "FL Reads") + facet_grid(~structural_category,scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          legend.position = "bottom") + 
    labs(x = " ", y = "FL reads", title = gene) + 
    scale_colour_discrete(name = "", label = c(label_group("Case"),label_group("Control")))
  
  return(p)
}