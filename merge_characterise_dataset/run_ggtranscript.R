suppressMessages(library(ggtranscript))
suppressMessages(library(ggplot2))
suppressMessages(library(rtracklayer))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(wesanderson))

replace_pbID <- function(PbID, gene){
  if(grepl("PB",PbID)){
    isoform <- word(PbID, c(3), sep = fixed("."))
    return(paste0("LR.",gene,".",isoform))
  }else{
    return(PbID)
  }
}

# inputDir = directory path of rnaseq aligned SJ.out.bed files
# subGTF = GTF to be subsetted with seqnames, start, and end
ggTranRNASeqPlots <- function(inputDir, inputgtf, transcriptID, filterCount=0){
  #rnaseq = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/C_RNASeq/2_aligned/All/mergedAll.SJ.out.bed")
  #rnaseq = rnaseq %>% filter(V1 == pGTF$seqnames[1],  V2 >= min(pGTF$start), V3 <= max(pGTF$end))
  #rnaseqCov = rnaseq %>% group_by(V1,V2,V3) %>% tally(V7)
  #meanGeneJuncFiles <- rnaseqCov 
  
  GTF <- as.data.frame(subset(inputgtf, transcript_id == transcriptID & type == "exon")) %>% mutate(structural_category = "RNA-Seq")
  CDS <- as.data.frame(subset(inputgtf, transcript_id == transcriptID & type == "CDS"))
  
  # read in SJ.out.bed files from RNA-Seq alignment
  message("Reading files")
  juncFilesNames = list.files(path = inputDir, pattern = "SJ.out.bed", full.names = T)
  print(juncFilesNames)
  juncFiles <- lapply(juncFilesNames, function(x) data.table::fread(x, header = "auto", sep = "\t"))
  
  # filter coordinates by the start and end of the transcriptID
  geneJuncFiles <- lapply(juncFiles, function(x) x %>% filter(V1 == GTF$seqnames[1],  V2 >= min(GTF$start), V3 <= max(GTF$end)))
  names(geneJuncFiles) <- word(list.files(path = inputDir, pattern = "SJ.out.bed"), c(1), sep = fixed("."))
  geneJuncFiles <- bind_rows(geneJuncFiles, .id = "sample") 
  
  # determine mean counts across samples
  meanGeneJuncFiles <- aggregate(geneJuncFiles$V7, by=list(geneJuncFiles$V1,geneJuncFiles$V2,geneJuncFiles$V3), FUN=mean)
  
  # format for plot
  meanGeneJuncFiles$strand = GTF$strand[1]
  colnames(meanGeneJuncFiles) <- c("seqnames","start","end","count","strand")
  meanGeneJuncFiles <- meanGeneJuncFiles %>% dplyr::mutate(transcript_id = transcriptID)
  
  # filter reads for plotting
  message("Filtering reads less than ", filterCount)
  meanGeneJuncFiles <- meanGeneJuncFiles %>% filter(count >= filterCount)
  
  
  p <- GTF %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(fill = "white", height = 0.25) +
    geom_range(data = CDS) + 
    geom_intron(data = to_intron(GTF, "transcript_id")) + 
    geom_junction(data = meanGeneJuncFiles, aes(size = count), junction.y.max = 0.5, colour = "purple") +
    geom_junction_label_repel(data = meanGeneJuncFiles, aes(label = round(count, 0)), segment.color = 'transparent',junction.y.max = 0.5) +
    scale_size_continuous(range = c(0.1, 1), guide = "none") + 
    theme_classic() + labs(y ="") +
    theme(legend.position = "None", 
          axis.line.x = element_line(colour = "grey80"),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.border = element_rect(fill = NA, color = "grey50", linetype = "dotted"),
          axis.text.y= element_text(size=12, colour = "white"),
          strip.text.y = element_text(size = 12, color = "black"),
          strip.background = element_rect(fill = "white", colour = "grey50"),
          axis.ticks.y=element_blank()) +
    facet_grid(rows = vars(structural_category))
  
  return(p)
}


# inputPfam: dataframe with 3 columns: <domain> <gene> <coordinates>
  # i.e <Apolipoprotein> <Apoe> <chr7:19696464-19697099> 
ggTranPlots <- function(inputgtf,classfiles,isoList,colours,lines,selfDf=NULL,gene=NULL,inputPfam=NULL,rnaseqDir=NULL,rnaseqTransID=NULL,rnaseqFilterCount=10,simple=FALSE,
                        inputCpat=NULL,cpatSpecies=NULL){
  
  if(any(duplicated(isoList))){
    print("Duplicated isoforms")
    print(isoList[duplicated(isoList)])
  }
  
  GTF <- as.data.frame(subset(inputgtf, transcript_id %in% isoList & type == "exon"))
  relCols <- c("isoform","structural_category","associated_transcript")
  xaxislevelsF1 <- c("Reference","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  GTF <- merge(GTF,classfiles[classfiles$isoform %in% isoList, relCols], by.x = "transcript_id", by.y = "isoform", all = T) %>% 
    mutate(structural_category = ifelse(grepl("ENSM",transcript_id),"Reference",as.character(structural_category))) %>%
    mutate(transcript_id = factor(transcript_id, levels = isoList))
  
  CDS <- as.data.frame(subset(inputgtf, transcript_id %in% isoList & type == "CDS"))
  
  if(!is.null(inputPfam)){
    if(is.null(gene)){
      stop("Error: Require gene name input")
    }
    Pfam <- inputPfam %>% filter(associated_gene == gene) %>% dplyr::rename("transcript_id" = "domain") %>% 
      mutate(seqnames = word(coordinates,c(1),sep = ":"),
             start = as.numeric(word(word(coordinates,c(2),sep = ":"),c(1),sep="-")),
             end = as.numeric(word(word(coordinates,c(2),sep = ":"),c(2),sep="-")),
             strand = GTF$strand[1],
             structural_category = "Pfam",
             Category = "Pfam"
             )%>% 
      select(-coordinates) 
    Pfam$transcript_id = factor(Pfam$transcript_id, levels = unique(Pfam$transcript_id))
  }
  
  if(!is.null(selfDf)){
    
    structural_colours <- rbind(
      data.frame(
        structural_category = c("Reference","FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron","coding","non-coding","noORF"),
        structural_col = c("#0C0C78","#00BFC4",alpha("#00BFC4",0.3),"#F8766D",alpha("#F8766D",0.3),"grey1","grey2","grey3","grey4","grey5",wes_palette("Darjeeling1")[2],wes_palette("Royal1")[2],"#0C0C78")
      ),
      data.frame(structural_category = selfDf[selfDf$Category == "DTE","colour"], 
                 structural_col = selfDf[selfDf$Category == "DTE","colour"])
    )
    
    if(!is.null(inputPfam)){
      structural_colours <- structural_colours %>% add_row(structural_category = "Pfam", structural_col = "#1400FA")
    }
    
    GTF <- merge(GTF,selfDf, by.x = "transcript_id", by.y = "Isoform", all.x = TRUE) 
    if(!is.null(inputPfam)){GTF <- plyr::rbind.fill(GTF, Pfam)}
    GTF <- merge(GTF, structural_colours, by = "structural_category",all.x=T) %>%
      mutate(structural_col = ifelse(Category != "DTE",as.character(structural_col),as.character(colour))) %>%
      mutate(structural_category = ifelse(Category != "DTE",as.character(structural_category),as.character(colour)))

    if(!is.null(inputPfam)){
      GTF$structural_category = factor(GTF$structural_category, levels = c(xaxislevelsF1,selfDf[selfDf$Category == "DTE","colour"],"Pfam"))
      GTF$transcript_id <- factor(GTF$transcript_id, levels = unique(c(isoList,as.character(GTF[GTF$Category == "Pfam","transcript_id"]))), ordered = TRUE)
      GTF$Category <- factor(GTF$Category, levels = c(as.character(unique(selfDf$Category)),"Pfam"))
    }else{
      GTF$structural_category = factor(GTF$structural_category, levels = c(xaxislevelsF1,selfDf[selfDf$Category == "DTE","colour"]))
      GTF$transcript_id <- factor(GTF$transcript_id, levels = isoList, ordered = TRUE)
      GTF$Category <- factor(GTF$Category, levels = unique(selfDf$Category))
    }
    
  }
  
  coding_prob_colours <- function(num,cpatSpecies){
    
    if(cpatSpecies == "human"){
      threshold = 0.364
    }else{
      threshold = 0.44
    }
    
    if(is.na(num)){
      return("noORF")
    }else if(num < threshold){
      return("non-coding")
    }else if(num >= threshold){
      return("coding")
    }else{
      return("noORF")
    }
  }

  if(!is.null(inputCpat)){
    if(!cpatSpecies %in% c("human","mouse") ){
      stop("cpatSpecies argument required <human/mouse>")
    }else{
      GTF <- merge(GTF, inputCpat[,c("ID","Coding_prob")], by.x = "transcript_id", by.y = "ID", all.x = T)
      GTF <- GTF %>% mutate(PBtranscript = ifelse(grepl("PB.",transcript_id),word(transcript_id,c(3),sep=fixed(".")),"NA"))
      GTF <- GTF %>% mutate(transcript_id = ifelse(grepl("PB.",transcript_id), paste0("LR.",gene,".",PBtranscript),as.character(transcript_id)))
      GTF$structural_category <- unlist(lapply(GTF$Coding_prob, function(x) coding_prob_colours(x,cpatSpecies)))
    }

  }
  
  if(!is.null(gene)){
    GTF <- GTF %>% mutate(PBtranscript = ifelse(grepl("PB.",transcript_id),word(transcript_id,c(3),sep=fixed(".")),"NA"))
    GTF <- GTF %>% mutate(transcript_id = ifelse(grepl("PB.",transcript_id), paste0("LR.",gene,".",PBtranscript),as.character(transcript_id)))
  }
  
  gexons <- GTF %>% dplyr::filter(type == "exon")
  gintrons <- gexons %>% to_intron(group_var = "transcript_id")
  grescaled <- shorten_gaps(gexons, gintrons, group_var = "transcript_id")
  
  if(!is.null(selfDf)){
    
    p <- grescaled %>%
      dplyr::filter(type == "exon") %>%
      ggplot(aes(xstart = start,xend = end,y = transcript_id)) +
      geom_range(aes(fill = structural_category)) +
      labs(y ="") + 
      geom_intron(data = grescaled %>% dplyr::filter(type == "intron"),arrow.min.intron.length = 300) +
      facet_grid(rows = vars(Category), scales = "free_y", space='free') +
      scale_fill_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"])) +
      scale_colour_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"]))
    
    
    #as.data.frame(GTF)  %>%
    #  ggplot(aes(xstart = start,xend = end, y = transcript_id)) +
    #  geom_range(aes(fill = structural_category)) +
    #  geom_intron(data = to_intron(GTF, "transcript_id"),aes(strand = strand, colour = structural_category)) +
    #  labs(y ="") + 
    #  facet_grid(rows = vars(Category), scales = "free_y", space='free') +
    #  scale_fill_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"])) +
    #  scale_colour_manual(values = as.character(structural_colours[structural_colours$structural_category %in% unique(GTF$structural_category),"structural_col"]))
    
  }else{
    GTF$transcript_id <- factor(GTF$transcript_id, levels = unlist(lapply(isoList, function(x) replace_pbID(x,gene))), ordered = TRUE)
    

    if(isFALSE(simple)){
      p <-  grescaled %>%
        ggplot(aes(xstart = start,xend = end, y = transcript_id)) +
        geom_range(data = CDS) +
        geom_intron(data = grescaled %>% dplyr::filter(type == "intron"),arrow.min.intron.length = 300) +
        scale_fill_manual(values = colours) +
        labs(y ="") + 
        facet_grid(rows = vars(structural_category), space = 'free_y', scales = "free_y")
      
    }else{
           
      p <- grescaled %>% filter(type == "exon") %>%
        ggplot(aes(xstart = start,xend = end, y = transcript_id)) +
        geom_range(aes(fill = transcript_id)) +
        geom_intron(data = grescaled %>% dplyr::filter(type == "intron"),arrow.min.intron.length = 300)  +
        scale_fill_manual(values = colours) +
        labs(y ="") +
        theme_minimal() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.x=element_blank(), 
              axis.text.x=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position = "None") 
    }

  }
  
  if(isFALSE(simple)){
    p <- p + theme_classic() +
      theme(legend.position = "None", 
            axis.line.x = element_line(colour = "grey80"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            panel.border = element_rect(fill = NA, color = "grey50", linetype = "dotted"),
            axis.text.y= element_text(size=12),
            strip.text.y = element_text(size = 12, color = "black"),
            strip.background = element_rect(fill = "white", colour = "grey50"))
  }

  
  if(is.null(rnaseqDir)){
    return(p)
  }else{
    p2 <- ggTranRNASeqPlots(inputDir=rnaseqDir,inputgtf=inputgtf,transcriptID=rnaseqTransID,filterCount=rnaseqFilterCount)
    return(list(p,p2))
  }
}
