get_gene_information <- function(type){
  
  # library("biomaRt")
  # ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # attributes <- c("ensembl_gene_id","external_gene_name","chromosome_name","strand","end_position","start_position","gene_biotype")
  # gene.location <- biomaRt::getBM(attributes = attributes,filters = "chromosome_name",values = c(seq_len(19), "X", "Y"),mart = ensembl)
  # write.csv(gene.location,"gene.location.csv")
  
  gene.location <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/MethReg/gene.location.csv", row.names = 1)
  gene.location$strand[gene.location$strand == 1] <- "+"
  gene.location$strand[gene.location$strand == -1] <- "-"
  gene.location$chromosome_name <- paste0("chr", gene.location$chromosome_name)
  
  if (missing(type)) {
    gene.location <- gene.location %>%
      makeGRangesFromDataFrame(
        seqnames.field = "chromosome_name",
        start.field = "start_position",
        end.field = "end_position", keep.extra.columns = TRUE
      )
  }
  
  return(gene.location)
}


#' Create a Granges object from a genmic region string
#' @description Given a region name such as chr22:18267969-18268249, we will create a Granges
#' object
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @param names A region name as "chr22:18267969-18268249" or a vector of region names.
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' @export
#' @return A GRanges
make_granges_from_names <- function(names){
  names %>%
    data.frame %>%
    separate(col = ".",into = c("chr","start","end")) %>%
    makeGRangesFromDataFrame()
}

#' Create region name from Granges
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @importFrom stringr str_c
#' @importFrom dplyr %>%
#' @importFrom GenomicRanges start end seqnames
#' @param region A GenomicRanges object
#' @examples
#' regions.names <- c("chr22:18267969-18268249","chr23:18267969-18268249")
#' regions.gr <- make_granges_from_names(regions.names)
#' make_names_from_granges(regions.gr)
#' @export
#' @return A string 
make_names_from_granges <- function(region){
  str_c(
    region %>% seqnames %>% as.character,":",
    region %>% start %>% as.character,"-",
    region %>% end %>% as.character)
}

#' @title Get promoter genes using biomart
#' @description Subset a granges object to those overlapping promoter regions
#' (default +- 2kb away from TSS)
#' @noRd
#' @importFrom GenomicRanges promoters strand strand<-
get_promoter_regions <- function(
  genome,
  upstream = 2000,
  downstream = 2000
){
  
  genes <- get_gene_information()
  promoters.gr <- promoters(genes, upstream = upstream, downstream = downstream)
  strand(promoters.gr) <- "*"
  return(promoters.gr %>% unique)
}



get_region_target_gene_by_promoter_overlap <- function(
  regions.gr,
  genome,
  upstream=target.promoter.upstream.dist.tss,
  downstream=target.promoter.downstream.dist.tss
){
  
  # Get gene information
  gene.promoters <- get_promoter_regions(
    genome = genome,
    upstream = upstream,
    downstream = downstream
  )
  
  hits <- findOverlaps(
    query = regions.gr,
    subject = gene.promoters,
    ignore.strand = TRUE,
    select = "all"
  )
  
  # overlap region and promoter
  neargenes <- gene.promoters[subjectHits(hits)] %>% as.data.frame(row.names = NULL)
  
  regions.gr <- regions.gr[queryHits(hits)]
  neargenes <- cbind(
    neargenes[,c(
      "seqnames",
      "start",
      "end",
      "external_gene_name",
      "ensembl_gene_id"
    )]
  )
  
  colnames(neargenes)[1:3] <- c(
    "target_gene_chrom",
    "target_gene_start",
    "target_gene_end"
  )
  colnames(neargenes)[4] <- "target_gene_name"
  colnames(neargenes)[5] <- "target"
  
  regionID <- regions.gr %>% data.frame %>% dplyr::select(1:3)
  regionID <- paste0(regionID[[1]],":",regionID[[2]],"-",regionID[[3]])
  out <- dplyr::bind_cols(
    data.frame("regionID" = regionID, stringsAsFactors = FALSE),
    neargenes[,4:5]
  ) %>% tibble::as_tibble()
  out
  
  out <- get_distance_region_target(out, genome = "mm10")
  out$target_tss_pos_in_relation_to_region <- NULL
  out$region_pos_in_relation_to_gene_tss <- NULL
  
  #out <- out %>% dplyr::rename(target_symbol = .data$target_gene_name)
  return(out)
  
}

get_distance_region_target <- function(
  region.target,
  genome = c("hg38","hg19","mm10")
){
  genome <- match.arg(genome)
  
  if ( !all( c("target","regionID") %in%  colnames(region.target))) {
    stop("Input requires columns regionID (chrx:start:end) and target (ENSG)")
  }
  region.target.only <- region.target %>%
    dplyr::filter(!is.na(.data$regionID)) %>%
    dplyr::filter(!is.na(.data$target)) %>%
    dplyr::select(c("target","regionID")) %>%
    unique()
  
  
  # resize is used to keep only the TSS to calculate the distance to TSS
  genes.gr <- get_gene_information() %>% resize(1)
  
  # We only need to calculate the distance of genes in the input
  genes.gr <- genes.gr[genes.gr$ensembl_gene_id %in% region.target.only$target]
  
  # Adding new information
  region.target.only <- region.target.only %>%
    cbind(
      data.frame(
        "distance_region_target_tss" = NA,
        "target_tss_pos_in_relation_to_region" = NA,
        "region_pos_in_relation_to_gene_tss" = NA
      )
    )
  
  # If the gene has no information the distance is NA
  region.target.no.info <- region.target.only %>%
    dplyr::filter(!.data$target %in% genes.gr$ensembl_gene_id)
  
  # If the gene has information the distance will be calculated
  region.target.info <- region.target.only %>%
    dplyr::filter(.data$target %in% genes.gr$ensembl_gene_id)
  
  regions.gr <- make_granges_from_names(
    names = region.target.info$regionID
  )
  
  
  idx <- match(region.target.info$target,genes.gr$ensembl_gene_id)
  dist <- distance(
    regions.gr,
    genes.gr[idx]
  )
  
  region.target.info$distance_region_target_tss <- dist
  region.target.info$region_pos_in_relation_to_gene_tss <- ifelse(
    as.logical(strand(genes.gr[idx]) != "-"),
    ifelse(start(regions.gr) < start(genes.gr[idx]), "upstream", "downstream"),
    ifelse(start(regions.gr) < start(genes.gr[idx]), "downstream", "upstream")
  )
  
  region.target.info$target_tss_pos_in_relation_to_region <- ifelse(
    start(regions.gr) < start(genes.gr[idx]), "right", "left"
  )
  
  
  region.target.info$distance_region_target_tss <-
    region.target.info$distance_region_target_tss * ifelse(region.target.info$region_pos_in_relation_to_gene_tss == "downstream",1, -1)
  
  # output both results together
  region.target.only <- plyr::rbind.fill(region.target.info, region.target.no.info)
  
  # using target and region keys, map the distance to the original input
  region.target$distance_region_target_tss <-
    region.target.only$distance_region_target_tss[
      match(
        paste0(
          region.target$regionID,
          region.target$target
        ),
        paste0(
          region.target.only$regionID,
          region.target.only$target
        )
      )
      ]
  
  region.target$target_tss_pos_in_relation_to_region <-
    region.target.only$target_tss_pos_in_relation_to_region[
      match(
        paste0(
          region.target$regionID,
          region.target$target
        ),
        paste0(
          region.target.only$regionID,
          region.target.only$target
        )
      )
      ]
  
  region.target$region_pos_in_relation_to_gene_tss <-
    region.target.only$region_pos_in_relation_to_gene_tss[
      match(
        paste0(
          region.target$regionID,
          region.target$target
        ),
        paste0(
          region.target.only$regionID,
          region.target.only$target
        )
      )
      ]
  
  return(region.target)
}

register_cores <- function(cores){
  
  parallel <- FALSE
  if (cores > 1) {
    if (cores > parallel::detectCores()) cores <- parallel::detectCores()
    doParallel::registerDoParallel(cores)
    parallel = TRUE
  }
  return(parallel)
}


get_tf_in_region <- function(
  region,
  window.size = 0,
  genome = c("hg19","hg38","mm39","mm10"),
  p.cutoff = 1e-8,
  cores = 1,
  TF.peaks.gr = NULL,
  verbose = FALSE
) {
  
  parallel <- register_cores(cores)
  
  if (is(region,"character") | is(region,"factor")) {
    region.gr <- make_granges_from_names(region)
    region.names <- region
  } else if (is(region,"GenomicRanges")) {
    region.gr <- region
    region.names <- make_names_from_granges(region)
  }
  
  region.gr <- region.gr + (window.size/2)
  # region <- region.gr %>% resize(width(.) + window.size, fix = "center")
  
  genome <- match.arg(genome)
  
  if (is.null(TF.peaks.gr)) {
    
    if (min(IRanges::width(region.gr)) < 8) {
      stop("Minimun region size is 8, please set window.size argument")
    }
    
    opts <- list()
    
    if(genome %in% c("mm10","mm39")){
      opts[["species"]] <- 10090 # Mus musculus
    } else {
      opts[["species"]] <- 9606 # homo sapiens
    }
    # opts[["all_versions"]] <- TRUE
    PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
    motifs.names <- lapply(PFMatrixList, function(x)(TFBSTools::name(x)))
    names(PFMatrixList) <- motifs.names
    PFMatrixList <- PFMatrixList[grep("::|var",motifs.names,invert = TRUE)]
    
    if(verbose)  message("Evaluating ", length(PFMatrixList), " JASPAR ",ifelse(genome %in% c("mm10","mm39"),"Mouse","Human")," TF motifs")
    if(verbose)  message("This may take a while...")
    
    suppressWarnings({
      motif.matrix <- motifmatchr::matchMotifs(
        pwms = PFMatrixList,
        subject = region.gr,
        genome = genome,
        p.cutoff = p.cutoff
      ) %>% SummarizedExperiment::assay()
    })
    rownames(motif.matrix) <- region.names
    
    # remove motifs not found in any regions
    motif.matrix <- motif.matrix[,DelayedArray::colSums(motif.matrix) > 0, drop = FALSE]
    
    if (ncol(motif.matrix) == 0) {
      message("No motifs found")
      return(NULL)
    }
    
    if (is(motif.matrix, "lgCMatrix")) {
      motif.matrix <- motif.matrix[!duplicated(rownames(motif.matrix)),, drop = FALSE]
      motif.matrix <- motif.matrix %>% as.matrix() %>% as.data.frame()
    }
    
    if(verbose)  message("Preparing output")
    
    motifs.probes.df <- plyr::alply(
      colnames(motif.matrix),
      .margins = 1,
      function(colum.name){
        colum <- motif.matrix[,colum.name, drop = FALSE]
        regions <- rownames(colum)[which(colum %>% pull > 0)];
        tfs <- colum.name
        expand.grid(regions,tfs,stringsAsFactors = FALSE)
      }, .progress = "time",.parallel = parallel)
    motifs.probes.df <- dplyr::bind_rows(motifs.probes.df)
    colnames(motifs.probes.df) <- c("regionID","TF_symbol")
    
    motifs.probes.df$TF <- map_symbol_to_ensg(
      motifs.probes.df$TF_symbol,
      genome = genome
    )
    
    motifs.probes.df <- motifs.probes.df %>% na.omit
    
  } else {
    
    hits <- findOverlaps(TF.peaks.gr, region.gr, ignore.strand = TRUE)
    motifs.probes.df <- data.frame(
      "regionID" = region.names[subjectHits(hits)],
      "TF_symbol" = TF.peaks.gr$id[queryHits(hits)]
    )
    motifs.probes.df$TF <- map_symbol_to_ensg(motifs.probes.df$TF_symbol, genome=genome)
  }
  return(motifs.probes.df %>% unique)
}


map_symbol_to_ensg <- function(
  gene.symbol, 
  genome
){
  gene.location <- get_gene_information()
  ensembl_gene_id <- ifelse(gene.symbol %in% gene.location$external_gene_name, 
                            as.character(gene.location$ensembl_gene_id[match(gene.symbol, gene.location$external_gene_name)]), "NA")
  return(ensembl_gene_id)
}

map_ensg_to_symbol <- function(
  ensembl.gene.id,
  genome 
){
  gene.location <- get_gene_information()
  symbols <- gene.location[match(ensembl.gene.id,gene.location$ensembl_gene_id),]$external_gene_name
  return(symbols)
}

#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param ensembl.gene.id Gene ensembl ID. A character vectors
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @examples
#' data(gene.exp.chr21.log2)
#' gene.region <- map_ensg_to_region(rownames(gene.exp.chr21.log2))
#' @noRd
get_target_tss_to_region_distance <- function(
  regionID,
  ensembl.gene.id,
  genome
){
  region.gr <- make_granges_from_names(regionID)
  gene.tss.location <- get_gene_information()  %>% resize(1)
  gene.tss.location <- gene.tss.location[match(ensembl.gene.id,gene.tss.location$ensembl_gene_id),]
  distance <- paste0(format(distance(region.gr,gene.tss.location)/1000,scientific = F), " kbp")
  return(distance)
}

#' @param genome Human genome of reference. Options: hg38, hg19.
#' @param ensembl.gene.id Gene ensembl ID. A character vectors
#' @description Given a GRanges returns region name such as chr22:18267969-18268249
#' @examples
#' data(gene.exp.chr21.log2)
#' gene.region <- map_ensg_to_region(rownames(gene.exp.chr21.log2))
#' @noRd
map_ensg_to_region <- function(
  ensembl.gene.id,
  genome
){
  gene.location <- get_gene_information(type="csv")
  gene.location <- gene.location[match(ensembl.gene.id,gene.location$ensembl_gene_id),]
  region <- paste0(gene.location$chromosome_name,":",gene.location$start_position,"-",gene.location$end_position)
  return(region)
}

#' @examples
#' results <- data.frame(
#'   "regionID" = c("chr3:203727581-203728580","chr4:203727581-203728580"),
#'   "TF" = c("ENSG00000232886","ENSG00000232887"),
#'   "target" = c("ENSG00000232889","ENSG00000242887"),
#'   "quant_pval_metGrp" = c(0.96,0.96),
#'   "quant_pval_rna.tf" = c(0.2, 0.2),
#'   "quant_pval_metGrp:rna.tf" = c(0.2, 0.2)
#' )
#' results <- calculate_stage_wise_adjustment(results)
#' @noRd
calculate_stage_wise_adjustment <- function(results){
  
  interactiol.col <- grep("RLM_metGrp.rna.tf_pvalue",colnames(results),value = TRUE)
  results <- stage_wise_adjustment(results, interactiol.col)
  
  dnam.col <- grep("metGrp_pvalue",colnames(results),value = TRUE)
  results <- stage_wise_adjustment(results, dnam.col)
  
  tf.col <- grep("tf_pvalue$",colnames(results),value = TRUE)
  tf.col <- grep("metGrp",tf.col,invert = TRUE,value = TRUE)
  results <- stage_wise_adjustment(results, tf.col)
  
  return(results)
}

