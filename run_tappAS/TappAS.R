# TAPPAS Results for TargetedMouse

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("wesanderson"))
suppressMessages(library("stringr"))
suppressMessages(library("tibble"))
suppressMessages(library("edgeR"))
suppressMessages(library("cowplot"))

library(tidyr)

library(extrafont)
#font_install('fontcm')
loadfonts()

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=18,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 18,family="CM Roman"),
                 axis.text.x= element_text(size=16,family="CM Roman"),
                 axis.text.y= element_text(size=16,family="CM Roman"))

label_colour <- function(genotype){
  if(genotype == "WT"){colour = wes_palette("Royal1")[2]}else{
    if(genotype == "TG"){colour = wes_palette("Royal1")[1]}else{
      if(genotype == "mouse"){colour = wes_palette("Royal1")[4]}else{
        if(genotype == "novel"){colour = wes_palette("Darjeeling1")[4]}else{
          if(genotype == "known"){colour = wes_palette("Darjeeling1")[5]}else{
          }}}}}
  return(colour)
}

output_plot_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/Targeted_Transcriptome/"
input_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/Results"
phenotype <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% mutate(col_names = paste0(group,".",sample))

targeted.class.files <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt", sep = "\t", as.is = T, header = T)

# read in files generated from TAPPAS
files <- list.files(path = input_dir, pattern = ".tsv", full.names = T)
files <- lapply(files, function(x) read.table(x, sep = "\t", header = T))
names(files) <- list.files(path = input_dir, pattern = ".tsv")

# number of transcripts that are filtered for statistical purposes
files$tappAS_Transcripts_InputExpressionMatrix.tsv <- 
  merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], files$tappAS_Transcripts_InputExpressionMatrix.tsv, by.x = "isoform", by.y = "Id")

# tally of the number of transcripts filtered per gene 
files$tappAS_Transcripts_InputExpressionMatrix.tsv %>% group_by(associated_gene, structural_category, Filtered) %>% tally() %>% ggplot(., aes(x = associated_gene, y = n, fill = Filtered)) + geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms") + mytheme +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(labels = c("Retained","Removed due to low coverage"), values = c(wes_palette("Rushmore1")[3],wes_palette("Rushmore1")[5])) + theme(legend.position = c(0.85,0.8))

files$tappAS_Transcripts_InputExpressionMatrix.tsv %>% group_by(associated_gene, structural_category, Filtered) %>% tally() %>% filter(Filtered != "NO") %>% ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + labs(x = "Target Gene", y = "Number of Isoforms Removed") + mytheme +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = c(0.8,0.8))


# Normalised Gene Expression Counts (Already filtered for low expression counts) 
Norm_transcounts <- merge(targeted.class.files[,c("isoform","associated_gene","associated_transcript","structural_category")], files$input_normalized_matrix.tsv, by.x = "isoform", by.y = 0) %>% reshape2::melt() %>% 
  left_join(., phenotype, by = c("variable" = "sample")) %>% 
  mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG")),
         structural_category=recode(structural_category, 
                                    `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",
                                    `novel_in_catalog`="NIC"),
         Isoform = paste0(isoform,"_",structural_category))

GeneExp <- Norm_transcounts %>% group_by(associated_gene,variable) %>% dplyr::summarise(Exp = sum(value)) %>%
  left_join(., phenotype, by = c("variable" = "sample"))


plot_mergedexp <- function(InputGene,IsoformID){
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

plot_transexp <- function(InputGene,IsoformID){
  df <-  Norm_transcounts  %>% filter(associated_gene == InputGene) %>% left_join(., phenotype, by = c("variable" = "sample")) %>% 
    mutate(group = factor(group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG")),
           structural_category=recode(structural_category, 
                                      `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM",
                                      `novel_in_catalog`="NIC"),
           Isoform = paste0(isoform,"_",structural_category))
  #df$group <- factor(df$group, levels = c("CONTROL", "CASE"),labels = c("WT", "TG"))
  #df <- df %>% mutate(structural_category=recode(structural_category, `full-splice_match`="FSM",`novel_not_in_catalog`="NNC",`incomplete-splice_match`="ISM"))
  
  p <- ggplot(df, aes(x = time, y = value, colour = Isoform)) + geom_point() + 
    facet_grid(~group) +
    stat_summary(data=df, aes(x=time, y=value, group=Isoform), fun ="mean", geom="line", linetype = "dotted") +
    mytheme + labs(x = "Age (months)", y = "Normalised Isoform Expression",title = paste0(InputGene,"\n\n")) +
    theme(strip.background = element_blank(), legend.position = "right",plot.title = element_text(hjust = 0.5, size = 16,face = "italic"),panel.spacing = unit(2, "lines"))
  return(p)
}

IFtrans_plot <- function(InputGene){
  dat <- Norm_transcounts  %>% filter(associated_gene == InputGene)
 
   # group by the samples
  dat <- dat %>% group_by(variable) %>% mutate(IF = value/sum(value) * 100)
  # group by isoform, age and group 
  dat2 <- aggregate(.~Isoform+time+group, dat[,c("Isoform","time","group","IF")], mean) %>% mutate(IF = replace_na(IF, 0)) %>%
    mutate(time = as.factor(time)) 
  p <- ggplot(dat2, aes(x = time, y = IF, fill = Isoform)) + geom_bar(stat = "identity") + facet_grid(~group) + 
    labs(x = "Age (months)", y = "Isoform Fraction (%)", title = paste0(InputGene,"\n\n")) + mytheme + 
    guides(fill=guide_legend(ncol=3,bycol=TRUE)) + 
    theme(legend.position="bottom", strip.background = element_blank(),plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  
  p1 <- ggplot(dat2, aes(x = time, y = IF, colour = Isoform, group = Isoform)) + 
    geom_line() + geom_point() + facet_grid(~group) + 
    labs(x = "Age (months)", y = "Isoform Fraction (%)", title = paste0(InputGene,"\n\n")) + mytheme + 
    guides(fill=guide_legend(ncol=3,bycol=TRUE)) + 
    theme(legend.position="bottom", strip.background = element_blank(),plot.title = element_text(hjust = 0.5, size = 16,face = "italic"))
  return(p1)
}


# Isoform Fraction Plots
IF_plots <- lapply(unique(GeneExp$associated_gene), function(gene) IFtrans_plot(gene))
IF_plotsgrobs <- lapply(IF_plots, ggplotGrob)
pdf (paste0(output_plot_dir,"/DifferentialAnalysis_IF.pdf"), width = 10, height = 15)
for(i in 1:20){print(plot_grid(IF_plotsgrobs[[i]]))}
dev.off()

IF_plots <- lapply(unique(GeneExp$associated_gene), function(gene) IFtrans_plot(gene))
IF_plotsgrobs <- lapply(IF_plots, ggplotGrob)
pdf (paste0(output_plot_dir,"/DifferentialAnalysis_IF2.pdf"), width = 10, height = 15)
for(i in 1:20){print(plot_grid(IF_plotsgrobs[[i]]))}
dev.off()


# Gene Expression Plots 
geneexp_plots <- lapply(unique(GeneExp$associated_gene), function(gene) plot_mergedexp(gene,"NA"))
geneexp_plotsgrobs <- lapply(geneexp_plots , ggplotGrob)
names(geneexp_plotsgrobs) <- unique(GeneExp$associated_gene)

# Target Expression Plots
transexp_plots <- lapply(unique(GeneExp$associated_gene), function(gene) plot_transexp(gene) + theme(legend.position = "none"))
transexp_plotsgrobs <- lapply(transexp_plots, ggplotGrob)
names(transexp_plotsgrobs) <- unique(GeneExp$associated_gene)

Norm_transcounts
# Transcript DEA 
transDEA_plots <- lapply(IsoformDEA, function(i) plot_mergedexp("NA",i))
transDEA_plotsgrobs <- lapply(transDEA_plots, ggplotGrob)

pdf (paste0(output_plot_dir,"/DifferentialAnalysis_Zoom.pdf"), width = 10, height = 15)
#plot_transexp("App") + guides(colour=guide_legend(ncol=3,bycol=TRUE)) + theme(legend.position="bottom")
for(gene in targeted.class.files[targeted.class.files$isoform %in% IsoformDEA, "associated_gene"]){
  print(plot_transexp(gene) + guides(colour=guide_legend(ncol=3,bycol=TRUE)) + theme(legend.position="bottom"))}
plot_grid(transDEA_plotsgrobs[[1]],transDEA_plotsgrobs[[2]],transDEA_plotsgrobs[[3]],transDEA_plotsgrobs[[4]], labels = c("a","b","c","d"), ncol = 2, scale = 0.9)
dev.off()

pdf (paste0(output_plot_dir,"/TargetedDifferentialAnalysis.pdf"), width = 10, height = 15)
### Gene expression Analysis
# EOAD, Regulator of AD
plot_grid(geneexp_plotsgrobs$App, geneexp_plotsgrobs$Mapt, geneexp_plotsgrobs$Fyn, geneexp_plotsgrobs$Trpa1, geneexp_plotsgrobs$Vgf, NULL, 
          labels = c("a","b","c","d","e"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
# GWAS
plot_grid(geneexp_plotsgrobs$Apoe, geneexp_plotsgrobs$Abca1, geneexp_plotsgrobs$Abca7, geneexp_plotsgrobs$Bin1, geneexp_plotsgrobs$Cd33, geneexp_plotsgrobs$Picalm, 
          labels = c("a","b","c","d","e","f"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
plot_grid(geneexp_plotsgrobs$Ptk2b, geneexp_plotsgrobs$Sorl1, geneexp_plotsgrobs$Trem2, NULL,NULL,NULL, 
          labels = c("g","h","i"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
# FTD/ALS
plot_grid(geneexp_plotsgrobs$Snca, geneexp_plotsgrobs$Fus, geneexp_plotsgrobs$Tardbp, NULL,NULL,NULL, 
          labels = c("a","b","c"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
# EWAS
plot_grid(geneexp_plotsgrobs$Ank1, geneexp_plotsgrobs$Rhbdf2,NULL,NULL,NULL,NULL, 
          labels = c("a","b"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
dev.off()
## Transcript
#plot_grid(transexp_plotsgrobs$App, transexp_plotsgrobs$Mapt, transexp_plotsgrobs$Apoe, transexp_plotsgrobs$Bin1, transexp_plotsgrobs$Snca, transexp_plotsgrobs$Trem2, 
          #labels = c("a","b","c","d","e","f"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
#plot_grid(transexp_plotsgrobs$Abca1, transexp_plotsgrobs$Abca7, transexp_plotsgrobs$Cd33, transexp_plotsgrobs$Clu, transexp_plotsgrobs$Fus, transexp_plotsgrobs$Fyn, 
          #labels = c("g","h","i","j","k","l"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
#plot_grid(transexp_plotsgrobs$Sorl1, transexp_plotsgrobs$Tardbp, transexp_plotsgrobs$Trpa1, transexp_plotsgrobs$Vgf, transexp_plotsgrobs$Ptk2b, transexp_plotsgrobs$Picalm, 
#          labels = c("m","n","o","p","q","r"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
#plot_grid(transexp_plotsgrobs$Ank1, transexp_plotsgrobs$Rhbdf2,NULL,NULL,NULL,NULL, 
#          labels = c("s","t"), label_size = 30, label_fontfamily = "CM Roman", ncol = 2, scale = 0.9)
#dev.off()

