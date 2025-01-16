descriptives_summary <- function(class_files_all, class_files_annoGenes) {

  annoGenesStats <- list(
  novelTrans = class_files_annoGenes[class_files_annoGenes$associated_transcript == "novel",],
  annoTrans = class_files_annoGenes[class_files_annoGenes$associated_transcript != "novel",],
  NIC = class_files_annoGenes[class_files_annoGenes$structural_category == "NIC",],
  NNC = class_files_annoGenes[class_files_annoGenes$structural_category == "NNC",]
  )

  # Initialize an empty list to store results
  results <- list()
  
  # Number of total transcripts and genes
  results[["Total Genes"]] <- length(unique(class_files_all$associated_gene))
  results[["Total Transcripts"]] <- nrow(class_files_all)

  # Annotated genes and transcripts to known genes
  results[["Annotated Genes"]] <- length(unique(class_files_annoGenes$associated_gene))
  results[["Total Transcripts to Annotated Genes"]] <- nrow(class_files_annoGenes)
  
  # Length summary for known genes
  results[["Mean Transcript Length"]] <- paste0(round(mean(class_files_annoGenes$length), 2), 
                                                " (sd = ", round(sd(class_files_annoGenes$length), 2), ")")
  results[["Transcript Length Range"]] <- paste0(round(min(class_files_annoGenes$length), 2), 
                                                 " - ", round(max(class_files_annoGenes$length), 2))
  
  # Number of transcripts summary per gene
  numIsogeneTally <- class_files_annoGenes %>%
    group_by(associated_gene) %>%
    tally()
  results[["Mean Number of Isoforms per Gene"]] <- paste0(round(mean(numIsogeneTally$n), 2), 
                                                         " (sd = ", round(sd(numIsogeneTally$n), 2), ")")
  results[["Genes with 10 or More Isoforms"]] <- paste0(nrow(numIsogeneTally[numIsogeneTally$n >= 10,]), 
                                                       " (", round(nrow(numIsogeneTally[numIsogeneTally$n >= 10,]) / 
                                                                    length(unique(numIsogeneTally$associated_gene)) * 100, 2), "%)")
  
  # Exon summary
  results[["Mean Number of Exons"]] <- paste0(round(mean(class_files_annoGenes$exons), 2), 
                                              " (sd = ", round(sd(class_files_annoGenes$exons), 2), ")")
  meanExonGene <- aggregate(class_files_annoGenes[,"exons"], 
                            list(class_files_annoGenes$associated_gene), mean)
  results[["Mean Exons per Gene"]] <- paste0(round(mean(meanExonGene$x), 2), 
                                             " (sd = ", round(sd(meanExonGene$x), 2), ")")
  
  # Novel and known transcript summaries
  results[["Total Transcripts to Annotated Genes"]] <- nrow(class_files_annoGenes)
  results[["Novel Transcripts"]] <- paste0(nrow(annoGenesStats$novelTrans), " (", 
                                           round(nrow(annoGenesStats$novelTrans) / 
                                                 nrow(class_files_annoGenes) * 100, 2), "%)")
  results[["Annotated Genes with Novel Transcripts"]] <- paste0(length(unique(annoGenesStats$novelTrans$associated_gene)), " (", 
                                                                round(length(unique(annoGenesStats$novelTrans$associated_gene)) / 
                                                                      length(unique(class_files_annoGenes$associated_gene)) * 100, 2), "%)")
  
  results[["Mean Novel Transcript Length"]] <- paste0(round(mean(annoGenesStats$novelTrans$length), 2), 
                                                      " (sd = ", round(sd(annoGenesStats$novelTrans$length), 2), ")")
  results[["Mean Novel Exons"]] <- paste0(round(mean(annoGenesStats$novelTrans$exons), 2), 
                                          " (sd = ", round(sd(annoGenesStats$novelTrans$exons), 2), ")")
  
  results[["Known Transcripts"]] <- paste0(nrow(annoGenesStats$annoTrans), " (", 
                                           round(nrow(annoGenesStats$annoTrans) / 
                                                 nrow(class_files_annoGenes) * 100, 2), "%)")
  
  results[["NIC Transcripts"]] <- paste0(nrow(annoGenesStats$NIC), " (", 
                                         round(nrow(annoGenesStats$NIC) / 
                                               nrow(annoGenesStats$novelTrans) * 100, 2), "%)")
  
  results[["NNC Transcripts"]] <- paste0(nrow(annoGenesStats$NNC), " (", 
                                         round(nrow(annoGenesStats$NNC) / 
                                               nrow(annoGenesStats$novelTrans) * 100, 2), "%)")
  
  # Convert results list to a data frame
  results_df <- as.data.frame(t(as.data.frame(results)))
  colnames(results_df) <- "Value"
  
  return(results_df)
} # <- Corrected final closing parenthesis for the function

monoExonic_summary <- function(preFilteredClassFile, demuxExpression){
  
  # subset to mono exonic isoforms in original SQANTI classification file
  monoExonic <- preFilteredClassFile[preFilteredClassFile$exons == 1,]
  
  # expresion data (totaFL of all reads across sample)
  demuxTotal <- demuxExpression %>% mutate(TotalFL = rowSums(select(., -id)))
  monoExonicdemux <- demuxTotal %>% filter(id %in% monoExonic$isoform)
  if(isFALSE(nrow(monoExonicdemux) == nrow(monoExonic))){
    message("WARNING: number of mono-exonic isoforms in classification file not the same as the expression data")
  }
  
  # initiate output list summary
  monoExonicStats <- list()
  
  monoExonicStats[["Total number of mono-exonic transcripts"]] <- nrow(monoExonic)
  
  monoExonicStats[["Percentage of mono-exonic transcripts out of all transcripts (before SQANTI filtering) "]] <- 
    round(nrow(monoExonic)/nrow(input$preFilteredClassFile) * 100,2)
  
  monoExonicStats[["Total number of FL reads of mono-exonic transcripts"]] <- sum(monoExonicdemux$TotalFL)
  
  monoExonicStats[["Percentage of FL reads of mono-exonic transcripts"]] <- 
    round(sum(monoExonicdemux$TotalFL)/sum(demuxTotal$TotalFL) * 100,2)
  
  
  monoExonicSqantiKept <- monoExonic[monoExonic$filter_result == "Isoform",]
  monoExonicStats[["Number of mono-exonic transcripts kept after SQANTI filtering"]] <- nrow(monoExonicSqantiKept)
  
  monoExonicStats[["Percentage of mono-exonic transcripts kept out of all mono-exonic transcripts"]] <- 
    round(nrow(monoExonicSqantiKept)/nrow(monoExonic) * 100,2)
  
  monoExonicMeaningful <- monoExonicSqantiKept[!monoExonicSqantiKept$structural_category %in% c("genic_intron","intergenic","genic"),]
  
  monoExonicStats[["Number of mono-exonic transcripts kept after further filtering"]] <- nrow(monoExonicMeaningful)
  
  monoExonicStats[["Percentage of mono-exonic transcripts kept after further filtering"]] <- 
    round(nrow(monoExonicMeaningful)/nrow(monoExonic) * 100,2)
  
  results_df <- as.data.frame(t(as.data.frame(monoExonicStats)))
  colnames(results_df) <- "Value"
  
  return(results_df)
}
