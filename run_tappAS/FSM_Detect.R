suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("cowplot"))

mytheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=14,  family="CM Roman"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black",family="CM Roman"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.80, 0.80),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6),
                 legend.text = element_text(size = 12,family="CM Roman"),
                 axis.text.x= element_text(size=12,family="CM Roman"),
                 axis.text.y= element_text(size=12,family="CM Roman"),
                 strip.background = element_rect(colour="white",fill="white"))



# SQANTI filtered file
sqanti.filtered.names.files = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite_classification.txt"
sqanti_filtered_inputfile  <- read.table(sqanti.filtered.names.files, header = T)

# Target Genes and Groups for plots
TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")
ADReg_Genes = c("App","Mapt","Fyn","Trpa1","Vgf")
GWAS_Genes = c("Apoe","Abca1","Abca7","Bin1","Cd33","Picalm","Ptk2b","Sorl1","Trem2")
FTDEPI_Genes = c("Snca","Fus","Tardbp","Ank1","Rhbdf2")

class.files <- sqanti_filtered_inputfile %>% filter(toupper(associated_gene) %in% TargetGene) %>% arrange(length)
class.files$isoform <- factor(class.files$isoform, levels = class.files$isoform)

boxplot_expression <- function(input_gene){
  dat = class.files %>% filter(toupper(associated_gene) == input_gene, associated_transcript != "novel") %>%
    select(isoform, associated_transcript, structural_category, starts_with("FL.")) %>% reshape2::melt()

  p <- ggplot(dat, aes(x = isoform, y = value, fill = structural_category)) + geom_boxplot(aes(colour = structural_category, alpha = 0.8)) + 
    facet_grid(~associated_transcript,scales = "free_x") + 
    labs(x = "Isoform", title = input_gene, y = "FL Read Counts") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + mytheme + 
    guides(alpha="none") + theme(legend.position = "bottom") + scale_y_continuous(trans = "log10")
  return(p)
}

fsm_boxplot_expression <- function(input_gene){
  ncounts_FSM <-  class.files %>% filter(structural_category == "full-splice_match") %>% group_by(associated_transcript) %>% tally()
  
  dat = class.files %>% filter(structural_category == "full-splice_match") %>% 
    filter(associated_gene %in% input_gene) %>%
    select(isoform, associated_gene, associated_transcript, structural_category, starts_with("FL.")) %>% reshape2::melt() %>% 
    left_join(., ncounts_FSM) %>% 
    mutate(unique = ifelse(n == 1, "Yes","No")) 
  
  p <- ggplot(dat, aes(x = isoform, y = value)) + 
    geom_boxplot(aes(colour = unique), alpha = 0.8) + 
    facet_grid(~associated_gene,scales = "free_x") + scale_y_continuous(trans = "log10") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Isoform", y = "FL Read Counts") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + mytheme + theme(legend.position = "bottom")
  
  return(p)
}


TG_plots <- lapply(lapply(TargetGene, function(x) boxplot_expression(x)),ggplotGrob)
names(TG_plots) <- TargetGene

pdf(paste0(output_plot_dir,"/Expression.pdf"), width = 11, height = 8.5)
for(gene in TargetGene){print(plot_grid(TG_plots[[gene]]))}
fsm_boxplot_expression(ADReg_Genes)
fsm_boxplot_expression(GWAS_Genes)
fsm_boxplot_expression(FTDEPI_Genes)
dev.off()


# Create threshold to use most highly expressed FSM as annotation for RNA-Seq
calc = class.files %>% select(starts_with("FL.")) %>% mutate(n_detected = rowSums(. != 0), total = rowSums(.)) 
class.files <- cbind(class.files, calc[,c("n_detected","total")])

FSM = class.files %>% filter(structural_category == "full-splice_match")

# only remove if more than one FSM transcript for that gene 
length(which(FSM$associated_gene == "Snca"))

