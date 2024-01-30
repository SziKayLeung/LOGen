# Differential Isoform Usage Analysis R script
#
# Differential Isoform Usage analysis will be performed by using two different methods: 
#   edgeR spliceVariants function and
#   DEXSeq package
# tappAS functions(): https://github.com/ConesaLab/tappAS/blob/master/scripts/DIU.R

minorFoldfilterTappas <- function(data, gen, minorfilter, minorMethod=c("PROP","FOLD")){
  print ("Removing low expressed minor isoforms")
  moreOne <- names(which(table(gen) > 1))
  head(moreOne)
  iso.sel <- NULL
  iso.sel.prop <- NULL
  
  gene.sel <- NULL
  if(minorMethod=="FOLD"){
    for ( i in moreOne) {
      which(gen==i)
      gene.data <- data[which(gen==i),]
      row.names(gene.data) <- row.names(data)[which(gen==i)]
      # sum the number of reads for each isoform across all the samples
      isoSUM <- apply(gene.data, 1, sum)
      names(isoSUM) <- row.names(data)[which(gen==i)]
      # identify major isoform as the isoform with the highest abundance
      major <- names(which(isoSUM == max(isoSUM)))[1]
      minors <- names(which(isoSUM != max(isoSUM)))
      # for each minor isoform, divide major isoform FL read/minor isoform FL read
      # i.e. major isoform = 45, minor isoform = 11, therefore fold change = 4.09
      #div <- as.numeric(matrix(rep(gene.data[major,], length(minors)), ncol = ncol(data), length(minors), byrow = T)) / as.matrix(gene.data[minors,])
      
      major_index <- which(row.names(gene.data) == major)
      major_row <- gene.data[major_index, ]
      minors_rows <- gene.data[which(row.names(gene.data) %in% minors), ]
      div <- as.matrix(data.table::rbindlist(rep(list(major_row), length(minors))))/as.matrix(minors_rows)
      row.names(div) <- minors
      
      # keep if fold change < 2
      # note if the fold change would be smaller, if the minor isoform reads are higher (as number/big number = small number)
      is <- names(which(apply(div, 1, min, na.rm = T) < minorfilter))
      iso.sel <- c(iso.sel, major, is)
      
    }
  }else{
    for ( i in moreOne) {
      which(gen==i)
      gene.data <- data[which(gen==i),]
      
      # by proportion
      geneSUM <- apply(gene.data, 2, sum)
      proportion = t(t(gene.data)/geneSUM)
      is.prop = rownames(proportion[apply(proportion, 1, function(x) any(x>minorfilter)),,drop=F])
      iso.sel <- c(iso.sel, is.prop) 
      
    }}
  
  print(length(iso.sel))
  print ("Done")
  return(iso.sel)
}

### EdgeR spliceVariants

spliceVariant.DS <- function(raw.counts, feature_association, factors) {
  library(edgeR)
  tmm_factors = calcNormFactors(raw.counts)
  DGEList_object_counts = DGEList(counts = round(raw.counts),
                                  lib.size = colSums(raw.counts),
                                  group = as.factor(factors[,"Replicate"]), genes = NULL,
                                  norm.factors = tmm_factors, 
                                  remove.zeros = FALSE) 
  isoform_names <- row.names(feature_association)
  feature_association <- as.data.frame(feature_association)
  row.names(feature_association) <- isoform_names
  
  results = spliceVariants(DGEList_object_counts, geneID=as.character(feature_association[rownames(feature_association) %in% rownames(raw.counts),"id"]), 
                           dispersion=NULL, estimate.genewise.disp = TRUE)$table[,"PValue",drop=F]
  isoPerGene = table(feature_association["id"])
  genesMoreOneIso = isoPerGene[isoPerGene>1] 
  results_multi = results[names(genesMoreOneIso),,drop=F]
  return(results_multi)
}

### DEXseq

DEXSeq.DS <- function(raw.counts, feature_association, factors) {
  library(DEXSeq)
  sample_table <- data.frame(row.names = colnames(raw.counts),
                             condition = as.factor(factors[,1]))
  dxd <- DEXSeqDataSet(round(raw.counts), sampleData = sample_table, design = ~ sample + exon + condition:exon,
                       featureID = rownames(raw.counts), 
                       groupID = as.character(feature_association[rownames(raw.counts),"id"]))
  tmm_librarySize_factors = tmm_factors_function(raw.counts, factors)
  sizeFactors(dxd) = c(tmm_librarySize_factors, tmm_librarySize_factors)
  ## pass number of workers? - 2 minimum: , BPPARAM = MulticoreParam(workers = 4))
  #dxd <- estimateDispersions(dxd, BPPARAM = MulticoreParam(workers = 4))
  if (.Platform$OS.type == "unix") {
    bp_param <- MulticoreParam(workers=4);
  } else if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers=4);
  }
  dxd <- estimateDispersions(dxd, BPPARAM = bp_param);
  ## Run the test and get results for exon bins 
  dxd <- testForDEU(dxd)
  ## Summarizing results on gene level
  res <- DEXSeqResults(dxd)
  print(dim(res))
  print(head(res))
  pgq <- perGeneQValue(res, p = "pvalue")
  results <- data.frame(row.names = names(pgq), q.value = pgq) 
  return(results)
}


#### Calculation of TMM normalization factors 

tmm_factors_function <- function(counts, myfactors) {
  library(edgeR)
  myedgeR = DGEList(counts = counts, group = myfactors[,1])
  myedgeR = calcNormFactors(myedgeR, method = "TMM", refColumn = 1)  
  norm_factors_edgeR =  myedgeR$samples$norm.factors
  total = colSums(as.matrix(counts))
  norm_factors_edgeR = norm_factors_edgeR * (total/mean(total))
  return(norm_factors_edgeR)
}

#### Shift in the position of isoforms by level of expression

podiumChange <- function(norm.counts, feature_association, myfactors) {
  
  isoform_names <- row.names(norm.counts)
  norm.counts <- as.data.frame(norm.counts)
  row.names(norm.counts) <- isoform_names
  
  isoform_names <- row.names(feature_association)
  feature_association <- as.data.frame(feature_association)
  row.names(feature_association) <- isoform_names
  
  as.character(feature_association[rownames(norm.counts),"id"])

  factors = split(myfactors, myfactors[,"Replicate"])
  norm.counts_mean = data.frame(cond1=apply(norm.counts[,rownames(factors[[1]])], 1, function(x) mean(x)), 
                                cond2=apply(norm.counts[,rownames(factors[[2]])], 1, function(x) mean(x)), 
                                gene=feature_association[rownames(norm.counts),"id"])
  changeOfPosition = as.list(by(norm.counts_mean[,c(1,2)], norm.counts_mean[,"gene"], function(x) {
    prop = t(t(x)/colSums(x))*100
    prop[is.nan(prop)] <- 0
    prop = as.data.frame(prop)
    position = data.frame(cond1_pos = order(prop$cond1, decreasing = TRUE), cond2_pos = order(prop$cond2, decreasing = TRUE), delta=abs(prop$cond1-prop$cond2))
    position <<- position
    if (position[1,"cond1_pos"] != position[1,"cond2_pos"]) { r = TRUE }
    else {r = FALSE}
    return (list(podiumChange = r, totalChange = sum(position$delta)/2))
  }))
  A = unlist(sapply(changeOfPosition, function(x) x$podiumChange))
  B = unlist(sapply(changeOfPosition, function(x) x$totalChange))
  return(list(podiumChange=A, totalChange=B))
}
