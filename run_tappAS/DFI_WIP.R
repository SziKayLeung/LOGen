suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))

DFI = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/IsoSeq_Expression/Data/DFI/dfi_result.0675479350.tsv", header = T)
DFI <- cbind(DFI,data.frame(do.call('rbind', strsplit(as.character(DFI$gene),';',fixed=TRUE))) %>% `colnames<-`(c("associated_gene", "f1", "f2", "f3","id")))
DFI <- DFI %>% mutate(type = word(id,c(1),sep = fixed("_")), pbid = word(id,c(3),sep = fixed("_")))
DFI$associated_gene = str_replace_all(DFI$associated_gene, fixed(" "), "")

# adPvalue < 0.05
DFI_sig <- DFI %>% filter(adjPValue < 0.05)

# Number of transcripts/proteins with/without podium change 
DFI_sig %>% group_by(podiumChange) %>% tally()

# With podium change 
DFI_sig_pod = DFI_sig %>% filter(podiumChange == "TRUE")
DFI_sig_pod_gene = DFI_sig_pod %>% group_by(associated_gene) %>% tally()
length(unique(DFI_sig_pod$associated_gene))

DFI_sig_pod[DFI_sig_pod$associated_gene == "Nfatc2",]
common = intersect(DFI_sig_pod$associated_gene,tappasDIU$DIU_rnaseq_fc$gene)
length(unique(common))
DFI_sig_pod = DFI_sig %>% filter(associated_gene %in% common)

for(i in unique(DFI_sig_pod$f1)){
  cat("##################",i,"\n")
  #print(unique(DFI_sig_pod[DFI_sig_pod$f1 == i,"f2"]))
  print(unique(DFI_sig_pod[DFI_sig_pod$f1 == i,"type"]))
}


### 
DFI_sig_pod_plots = DFI_sig_pod  %>%  mutate(f1 = recode(f1,SIGNALP_EUK = "SignalP", "UTRsite" = "UTR Sites","TranscriptAttributes" = "Transcript Attributes"),
                                            f2 = recode(f2,DOMAIN = "Pfam Domains", "repeat" = "Repeat Regions", "SIGNAL" = "Signal Peptides", 
                                                        "3UTR_Length" = "3'UTR Length","5UTR_Length" = "5'UTR Length",
                                                        "3UTRmotif" = "3'UTR Motif","5UTRmotif" = "5'UTR Motif",
                                                        "COILED" = "Coiled-coil Regions","TRANSMEM" = "Transmembrane Regions",
                                                        "miRNA_Binding" = "miRNA Binding","polyA_site" = "PolyA-site", "PTM" = "Post-Translational Modifications",
                                                        "COMPBIAS" = "Compositional Bias", "ACT_SITE" = "Active Site","MOTIF" = "Motifs","BINDING" = "Binding Site")) 


DFI_sig_pod_plots %>% filter(type == "TRANS") %>% group_by(f1,f2) %>% tally() %>% ggplot(., aes(x = reorder(f2,n), y = n, fill = f1)) + geom_bar(stat = "identity") + facet_grid(rows = vars(f1), scales = "free") + mytheme + theme(legend.position = c(0.8,0.2), strip.text.y = element_blank()) + labs(x = "Features", y = "Number of Transcript Features with Differential Inclusion") + coord_flip() + scale_fill_discrete(name = "")

DFI_sig_pod_plots %>% filter(type == "PROTEIN") %>% group_by(f1,f2) %>% tally() %>% ggplot(., aes(x = reorder(f2, n), y = n, fill = f1)) + geom_bar(stat = "identity") + facet_grid(rows = vars(f1), scales = "free") + mytheme + theme(legend.position = c(0.8,0.8), strip.text.y = element_blank()) + labs(x = "Features", y = "Number of Protein Features with Differential Inclusion") + coord_flip() + scale_fill_discrete(name = "")

unique(DFI_sig_pod[DFI_sig_pod$f2 == "miRNA_Binding","associated_gene"])
DFI_sig_pod[DFI_sig_pod$associated_gene == "Hspb6",]

                       