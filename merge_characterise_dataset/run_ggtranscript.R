suppressMessages(library(ggtranscript))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

replace_pbID <- function(PbID, gene){
  isoform <- word(PbID, c(3), sep = fixed("."))
  return(paste0("LR.",gene,".",isoform))
}

ggTranPlots <- function(inputgtf,classfiles,isoList,colours,lines,selfDf=NULL,gene=NULL){
  
  if(any(duplicated(isoList))){
    print("Duplicated isoforms")
    print(isoList[duplicated(isoList)])
  }
  
  GTF <- as.data.frame(subset(inputgtf, transcript_id %in% isoList & type == "exon"))
  relCols <- c("isoform","structural_category","associated_transcript")
  xaxislevelsF1 <- c("Reference","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  GTF <- merge(GTF,classfiles[classfiles[["isoform"]] %in% isoList, relCols], by.x = "transcript_id", by.y = "isoform", all = T) %>% 
    mutate(structural_category = ifelse(grepl("ENSM",transcript_id),"Reference",as.character(structural_category))) %>%
    mutate(transcript_id = factor(transcript_id, levels = isoList))
  
  if(!is.null(selfDf)){
    structural_colours <- rbind(
      data.frame(
        structural_category = c("Reference","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron"),
        structural_col = c("#0C0C78","#00BFC4",alpha("#00BFC4",0.3),"#F8766D",alpha("#F8766D",0.3),"grey1","grey2","grey3","grey4","grey5")
      ),
      data.frame(structural_category = selfDf[selfDf$Category == "DTE","colour"], 
                 structural_col = selfDf[selfDf$Category == "DTE","colour"])
    )
    
    GTF <- merge(GTF,selfDf, by.x = "transcript_id", by.y = "Isoform", all.x = TRUE) 
    GTF <- merge(GTF, structural_colours, by = "structural_category",all.x=T) %>%
      mutate(structural_col = ifelse(Category != "DTE",as.character(structural_col),as.character(colour))) %>%
      mutate(structural_category = ifelse(Category != "DTE",as.character(structural_category),as.character(colour)))
    GTF$structural_category = factor(GTF$structural_category, levels = c(xaxislevelsF1,selfDf[selfDf$Category == "DTE","colour"]))
    GTF$transcript_id <- factor(GTF$transcript_id, levels = isoList, ordered = TRUE)
    GTF$Category <- factor(GTF$Category, levels = unique(selfDf$Category))
  }

  if(!is.null(gene)){
    GTF <- GTF %>% mutate(PBtranscript = ifelse(grepl("PB.",transcript_id),word(transcript_id,c(3),sep=fixed(".")),"NA"))
    GTF <- GTF %>% mutate(transcript_id = ifelse(grepl("PB.",transcript_id), paste0("LR.",gene,".",PBtranscript),as.character(transcript_id)))
  }
  
  if(!is.null(selfDf)){
    
    p <- as.data.frame(GTF)  %>%
      ggplot(aes(xstart = start,xend = end, y = transcript_id)) +
      geom_range(aes(fill = structural_category)) +
      geom_intron(data = to_intron(GTF, "transcript_id"),aes(strand = strand, colour = structural_category)) +
      labs(y ="") + 
      facet_grid(rows = vars(Category), scales = "free_y", space='free') +
      scale_fill_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"])) +
      scale_colour_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"]))
    
  }else{
    GTF$transcript_id <- factor(GTF$transcript_id, levels = unlist(lapply(isoList, function(x) replace_pbID(x,gene))), ordered = TRUE)
    p <-  as.data.frame(GTF)  %>%
      ggplot(aes(xstart = start,xend = end, y = transcript_id)) +
      geom_range(aes(fill = transcript_id)) +
      geom_intron(data = to_intron(GTF, "transcript_id"),aes(strand = strand)) +
      scale_fill_manual(values = colours) +
      labs(y ="") + 
      facet_grid(rows = vars(structural_category), space = 'free_y', scales = "free_y")
  }
  
  
  if(!missing(lines)){
    p <- p + geom_intron(data = to_intron(GTF, "transcript_id"),aes(strand = strand, colour = transcript_id)) +
      scale_colour_manual(values = lines)
  }
  
  p <- p + theme_classic() +
    theme(legend.position = "None", 
          axis.line.x = element_line(colour = "grey80"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.border = element_rect(fill = NA, color = "grey50", linetype = "dotted"),
          axis.text.y= element_text(size=12),
          strip.text.y = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill = "white", colour = "grey50"))
  
  return(p)
}
