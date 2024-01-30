# Szi Kay Leung: sl693@exeter.ac.uk
# Differential Isoform Usage Analysis R script from tappAS

minorFoldfilterTappas <- function(data, gen, minorfilter, minorMethod=c("PROP","FOLD")){
  print ("Removing low expressed minor isoforms")
  # list the genes with more than one isoform 
  # table(gen) = number of times gene repeated
  moreOne <- names(which(table(gen) > 1))
  
  iso.sel <- NULL
  iso.sel.prop <- NULL
  gene.sel <- NULL
  
  if(minorMethod=="FOLD"){
    for ( i in moreOne) {
      # subset from the normalised expression matrix, the isoforms associated with that gene
      which(gen==i)
      gene.data <- data[which(gen==i),]
      
      # by fold change
      # for each isoform, sum the expression across the sample
      isoSUM <- apply(gene.data, 1, sum)
      
      # determine the major isoform and mino isoform from the summed expression across all the samples
      major <- names(which(isoSUM == max(isoSUM)))[1]
      minors <- names(which(isoSUM != max(isoSUM)))
      
      # for each minor isoform across each sample, divide the expression by the expression of the major isoform
      # major/minor i.e. if small number then minor is a large proportion
      div <- as.numeric(matrix(rep(gene.data[major,], length(minors)), ncol = ncol(data), length(minors), byrow = T)) / as.matrix(gene.data[minors,])
      
      # determine the minimum expression for the isoform of all the samples, ignoring infinite 
      # note the isoform if the minimum expression < threshold 
      is <- names(which(apply(div, 1, min, na.rm = T) < minorfilter))
      iso.sel <- c(iso.sel, major, is)
      
    }
  }else{
    for ( i in moreOne) {
      # subset from the normalised expression matrix, the isoforms associated with that gene
      which(gen==i)
      gene.data <- data[which(gen==i),]
      
      # by proportion
      # for each sample, sum the expression of the isoforms and then determine the proportion 
      # keep only the isoforms that have a proportion higher than the user-defined threshold any of the samples
      geneSUM <- apply(gene.data, 2, sum)
      proportion = t(t(gene.data)/geneSUM)
      geneSUM <<- geneSUM
      proportion <<- proportion
      # for each isoform i.e. row; TRUE if any of the samples have an isoform proportion > threshold, otherwise drop
      is.prop = rownames(proportion[apply(proportion, 1, function(x) any(x>minorfilter)),,drop=F])
      # add the passed isoforms to the growing list of iso.sel
      iso.sel <- c(iso.sel, is.prop) 
      
    }}
  
  print(length(iso.sel))
  print ("Done")
  return(iso.sel)
}

# DIU_prefilter 
# mff = minimum foldchange
# filteringType = filtering by proportion or foldchange i.e. PROP or FOLD
DIU_prefilter <- function(dir, mff, filteringType){
  
  transMatrix <- NULL
  genes <- NULL

  # read transcript expression matrix normalized
  cat("\nReading normalized transcript matrix file data...")
  transMatrix = read.table(file.path(dir, "/Data/transcript_matrix.tsv"), row.names=1, sep="\t", header=TRUE)
  cat("\nRead ", nrow(transMatrix), " normalized transcripts expression data rows")
  
  # read gene transcripts map - this file also contains the transcript lengths, (geneName/transcript/length),
  cat("\nReading gene transcripts map...")
  geneTrans <- read.table(file.path(dir, "Data/gene_transcripts.tsv"), sep="\t", header=TRUE)
  genes <- data.frame(row.names = geneTrans$transcript, "id" = as.character(geneTrans$geneName))
  genes <<- genes

  ### filter matrix
  cat("\nRead ", nrow(transMatrix), " expression data rows")
  
  #Filter transMatrix - need genes
  # note same file as geneTrans except without the length
  infoGenes = read.table(file.path(dir, "Data/result_gene_trans.tsv"), sep="\t", quote=NULL, header=FALSE,  stringsAsFactors=FALSE)
  genes_trans = c()
  index = c()
  
  transMatrix = transMatrix[order(rownames(transMatrix)),]
  cat("\nIntersecting DIU information with transcript matrix...")
  # for each isoform with expression
  for(i in (rownames(transMatrix))){
    # find the corresponding gene name and index of where it is on that file
    if(i %in% infoGenes[,2]){
      genes_trans = c(genes_trans, infoGenes[which(infoGenes[,2]==i),1])
      index = c(index, which(rownames(transMatrix)==i))
    }
  }
  
  transMatrix = transMatrix[index,]
  
  # filter the new transcript matrix with the filtering method 
  cat(paste0("\nFiltering new transcript matrix by ",filteringType,"..."))
  trans = minorFoldfilterTappas(transMatrix, genes_trans, mff, minorMethod=filteringType)
  transMatrix=transMatrix[trans,]
  
  return(transMatrix)
}

#---------------------------------------------------------------------------------------------------
# Auxiliar internal functions: Model, REP, Formula0, Formula1  
#---------------------------------------------------------------------------------------------------

modelIso <- function(Formula, design, family, data1, Genes1, epsilon, Q, gen1){
  pval<-NULL
  g <- length(Genes1)
  dis <- as.data.frame(design$dis)
  mycolnames <- colnames(dis)
  for(i in 1:g)
  {
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) 
      print(paste(c("fitting gene", i, "out of", g), collapse = " "))
    
    zz <-data1[gen1==Genes1[i],]
    nt <- nrow(zz) 
    
    dis.gen <- REP(dis,nt)
    y <-  c(t(as.matrix(zz)))
    transcript<- factor(rep(c(1:nt),each = ncol(zz)))
    ydis <- cbind(y, dis.gen, transcript)
    
    model0 <- glm(Formula(mycolnames,design$edesign), data=ydis,  family = family, epsilon = epsilon )
    model1 <- glm(Formula1(mycolnames), data=ydis,  family = family, epsilon = epsilon )
    
    if(family$family=="gaussian")  {
      pvali <- anova(model0, model1, test="F")[2,6] }
    else {
      pvali <- anova(model0,model1, test = "Chisq")[2,5] }
    names(pvali) <- Genes1[i]
    pval <- c(pval, pvali)
  }
  p.adjusted <- p.adjust(pval, method="fdr")
  num.genes <- sum(p.adjusted<Q, na.rm=TRUE)
  selected.genes <-names(sort(p.adjusted)[1:num.genes])
  
  results = list(p.adjusted, selected.genes)
  names(results) = c("p.adjusted", "selected.genes")
  
  return(results)
  
}

REP <- function(D,k)
{
  r<-nrow(D)
  c<-ncol(D)
  DD<-NULL
  for(i in 1:c)
  {
    DDi <- rep(D[,i],k)
    DD <- cbind(DD,DDi)
  }
  colnames(DD)<-colnames(D)
  as.data.frame(DD)
}

#---------------------------------------------------------------------------

Formula0 <- function(names,edesign=NULL)
{
  formula <- "y~"
  
  if(length(names)==1){ formula=paste(formula,names[1],"+ transcript") }
  else if(length(names)>1)
  {
    
    for (i in 1:(length(names)))
    {
      formula <- paste(formula,names[i],"+")
    }
    formula <- paste(formula, "transcript")
  }
  formula <- as.formula(formula)
  formula
}

#---------------------------------------------------------------------------

Formula1 <- function(names)
{
  formula <- "y~"
  
  if(length(names)==1){ formula=paste(formula,names[1],"* transcript") }
  else if(length(names)>1)
  {
    formula <- paste(formula,"(" )
    for (i in 1:(length(names)-1))
    {
      formula <- paste(formula,names[i],"+")
    }
    formula <- paste(formula,names[length(names)])
    formula <- paste(formula, ") * transcript")
  }
  formula <- as.formula(formula)
  formula
}

Formula000 <- function(names,edesign)
{
  formula <- "y~"
  if(length(names)==1){ formula=paste(formula,names[1],"+ transcript") }
  else if(length(names)>1)
  {
    name.time = colnames(edesign)[1]
    names.inter = paste(name.time,c("x","2x","3x","4x","5x"), sep="")
    COL.inter=NULL
    i=1
    col.i = grep(names.inter[i],names)
    COL.inter=col.i
    while(length(col.i)!=0)
    {
      i=i+1
      col.i = grep(names.inter[i],names)
      COL.inter=c(COL.inter,col.i)
    }
    names1=names[-COL.inter]
    names2=names[COL.inter]
    
    formula <- paste(formula,"(" )
    for (i in 1:(length(names1)-1))
    {
      formula <- paste(formula,names1[i],"+")
    }
    formula <- paste(formula,names1[length(names1)])
    formula <- paste(formula, ") * transcript")
  }
  
  formula <- paste(formula, "+")
  
  if(length(names2)>1){
    for (i in 1:(length(names2)-1))
    {formula <- paste(formula, names2[i],"+") }
  }
  formula <- paste(formula, names2[length(names2)])
  formula <- as.formula(formula)
  formula
}

#---------------------------------------------------------------


Formula00 <- function(names,edesign)
{
  formula <- "y~"
  if(length(names)==1){ formula=paste(formula,names[1],"+ transcript") }
  else if(length(names)>1)
  {
    name.time = colnames(edesign)[1]
    names.conds = colnames(edesign)[3:ncol(edesign)]
    names.inter = paste(name.time,c("x","2x","3x","4x","5x"), sep="")
    COL.group = unique(unlist(sapply(names.conds, function(m) grep(m, names))))
    COL.inter = unique(unlist(sapply(names.inter, function(m) grep(m, names))))
    
    COL.out = unique(c(COL.group,COL.inter))
    
    names1=names[-COL.out]
    names2=names[COL.out]
    
    formula <- paste(formula,"(" )
    for (i in 1:(length(names1)-1))
    {
      formula <- paste(formula,names1[i],"+")
    }
    formula <- paste(formula,names1[length(names1)])
    formula <- paste(formula, ") * transcript")
  }
  
  formula <- paste(formula, "+")
  
  if(length(names2)>1){
    for (i in 1:(length(names2)-1))
    {formula <- paste(formula, names2[i],"+") }
  }
  formula <- paste(formula, names2[length(names2)])
  formula <- as.formula(formula)
  formula
}

# function for tappas IsoModel

IsoModel_tappas <- function(data, gen, design = NULL, degree = 2, Q = 0.05, min.obs = 6, minorFoldfilter = NULL, 
                            theta = 10, epsilon = 1e-05, triple = FALSE) 
{
  #---------------------------------------------------------------------------------------------------
  # data is a matrix containing isoform expression. Isoforms must be in rows and experimental conditions in columns
  # gen is a vector with the name of the gene each isoform belongs to
  #---------------------------------------------------------------------------------------------------
  
  Genes <- unique(gen)
  g <- length(Genes)
  
  # assuming binomial negative distribution
  family = negative.binomial(theta)

  #---------------------------------------------------------------------------------------------------
  # STEP -1: Remove cases with low expressed isoforms:
  #---------------------------------------------------------------------------------------------------
  
  print (paste(nrow(data), "transcripts"))
  print (paste(length(unique(gen)), "genes"))
  
  totaldata = data  #### NEW!
  totalgen = gen   #### NEW!
  
  #---------------------------------------------------------------------------------------------------
  #  STEP 0: Remove cases with 1 transcript:
  #---------------------------------------------------------------------------------------------------
  
  NT <- tapply(rownames(data),gen,length)
  Genes1 <- names(which(NT!=1))
  data1 <- data[gen%in%Genes1,]
  gen1 <- gen[gen%in%Genes1]
  Genes1 <- unique(gen1)
  print (paste(nrow(data1), "remaining transcripts to analyse DS"))   # changed
  print (paste(length(unique(Genes1)), "remaining genes to analyse DS")) # changed
  
  #---------------------------------------------------------------------------------------------------
  # STEP 1: Gene models comparison. 
  #---------------------------------------------------------------------------------------------------
  
  # make.design.matrix 
  results = NULL
  
  # epsilon = 1e-05
  # Q = 0.05
  # gen1 = genes with more than one transcripts 
  
  if (ncol(design$edesign)==3){
    results <- modelIso(Formula=Formula0, design, family, data1, Genes1, epsilon, Q, gen1)
  } 
  else{
    # for (k in 3:ncol(design$edesign)){
    #   singleG_dis = design$edesign[which(design$edesign[,k]==1), c(1,2,k)]
    #   dis = make.design.matrix(singleG_dis, degree = degree)
    #   results[[colnames(design$edesign)[k]]] <- function(Formula=Formula0, dis, family, data1, Genes1, epsilon, Q)
    # }
    if(triple){Formula=Formula000} else {Formula=Formula00}
    results <- modelIso(Formula=Formula, design, family, data1, Genes1, epsilon, Q, gen1)
  }
  
  
  #---------------------------------------------------------------------------------------------------
  # STEP 2: p.vector and T.fit for DE to the transcripts that belong to selected.genes
  #---------------------------------------------------------------------------------------------------
  # data2 <- data[gen%in%selected.genes,]
  # gen2 <- gen[gen%in%selected.genes]
  # pvector2 <- p.vector(data2, design, counts=counts, item="isoform")
  # Tfit2 <- T.fit(pvector2, item="isoform")
  # 
  #---------------------------------------------------------------------------------------------------
  # Output  
  #---------------------------------------------------------------------------------------------------
  
  #ISO.SOL <-list(data, gen, design, selected.genes, pvector2, Tfit2)
  #names(ISO.SOL)<- c("data","gen", "design", "DSG","pvector","Tfit")
  #ISO.SOL
  
  
  ISO.SOL <-list(totaldata, totalgen, design, results$selected.genes, results$p.adjusted)
  names(ISO.SOL)<- c("data","gen", "design", "DSG", "adjpvalue")
  ISO.SOL  
  
}


#-------------------------------------------------------
# PodiumChange_tappas
#-------------------------------------------------------

PodiumChange_tappas <- function(iso, only.sig.genes = FALSE, comparison=c("any","group","specific"), group.name="Ctr", time.points=0){
  gen <- iso$gen
  data <- iso$data
  edesign <- iso$design$edesign
  repvect = edesign[,2]
  
  if(only.sig.genes){
    data.clust<-as.matrix(data[gen%in%iso$DSG,])  
    sig.iso2<-rownames(data.clust)
    gen.sig.iso2 <- as.character(gen[rownames(data)%in%sig.iso2])
    # Here, there is not any mDSG because in this analysis it is considered only (>1 iso)
  }else{
    # remove rows that are na 
    gen2 <- names(which(tapply(rownames(data), gen, length) >1))
    data.clust<-as.matrix(data[gen%in%gen2,])  
    sig.iso2<-rownames(data.clust)
    gen.sig.iso2 <- as.character(gen[rownames(data)%in%sig.iso2])
  }
  
  #-------------------------------------------------------
  # Major Isoform identification
  #-------------------------------------------------------
  
  time.M <- tapply(edesign[,1],repvect,mean)
  groups.M <- apply(edesign[,3:ncol(edesign),drop=F],2,function(x){tapply(x,repvect,mean)})        ##### Added -> drop=F!!!!!
  
  unic <- unique(gen.sig.iso2)
  Mayor=NULL
  LIST = NULL
  TIMELIST = list()
  #if(comparison=="any"){for(group.name in colnames(groups.M)){TIMELIST[[group.name]] = c()}}
  for(i in 1:length(unic)){
    # zz = subset the isoform expression relating to that gene
    # M = find the major isoform with the summed higest expressed across all the samples
    # zzM subset the row of isoform expression relating to the major isoform
    zz<-data.clust[gen.sig.iso2==unic[i],]
    M <- MayorIso(zz)
    zzM<-zz[M==1,]
    if(length(zzM)>length(repvect)){
      #print(zzM)
      zzM <- zzM[1,]
      MzzM <- tapply(zzM,repvect,mean)
      zzm <- zz[M!=max(M),]
      if(length(zzm)==0){
        zzm <- zz[nrow(zz),]
      }
    }else{
      # with the major isoform, for each group (in this case, 4 groups with WT vs TG in 2 age points), determine the mean expression
      # MzzM = group mean expression of major isoform
      MzzM <- tapply(zzM,repvect,mean)
      zzm <- zz[M!=max(M),]
    }
    
    # also determine mean expression across the groups for the minor isoforms 
    # Mzzm = group mean expression of other minor isoforms
    if(is.null(nrow(zzm))) ni=1 else ni=nrow(zzm)
    
    if(ni==1) Mzzm=tapply(zzm,repvect,mean) else Mzzm <- t(apply(zzm, 1, function(x){tapply(x, repvect, mean)}))
    
    # dif = Major isoform (mean) - Minor isoform (mean)
    if(ni==1) dif=MzzM - Mzzm else dif <- t(apply(Mzzm, 1, function(x){MzzM-x}))
    
    # Comparison = "a" ----------------------------------------------    
    if(comparison=="any"){  
      if( any(dif<0) ){ 
        LIST <- c( LIST, unic[i])
        if(ni==1) t.points=names(dif[dif<0]) else t.points=colnames(dif[,apply(dif,2,function(x) any(x<0)), drop=F])
        r = data.frame(groups.M[t.points,,drop=F], t.point=time.M[t.points])
        for (group.name in colnames(groups.M)){
          if(length(r[which(r[,group.name]==1),'t.point'])>0) {TIMELIST[[group.name]][unic[i]] = paste(r[which(r[,group.name]==1),'t.point'], collapse = ",")}
        }
      }
    }
    
    # Comparison = "specific" ----------------------------------------------
    else if(comparison=="specific"){
      col <- groups.M[,colnames(groups.M)==group.name]
      if(ni==1) change <- all(dif[col==1 & time.M==time.points]<0)
      else change <- apply(dif[,col==1 & time.M==time.points,drop=F],1,function(x){all(x<0)})       ##### Added -> drop=F!!!!!
      if(any(change)) LIST <- c( LIST, unic[i]) } 
    # Comparison = group ----------------------------------------------
    else if(comparison=="group"){
      # WT vs TG (independent of age)
      mayors = NULL
      for (k in 3:ncol(edesign))
      {
        mayors = cbind(mayors, MayorIso(zz[,edesign[,k]==1]))
      } 
      #When all columns match, substraction with any of them will be 0:
      if (any(mayors-mayors[,1]!=0) )  LIST <- c( LIST, unic[i]) }                 #### any instead of all  !!!!
  }
  
  # lists of genes and isoforms:
  gen.L <- gen.sig.iso2 [gen.sig.iso2 %in% LIST]
  data.L <- data.clust[gen.sig.iso2 %in% LIST,]
  
  output <- list(LIST, TIMELIST, data.L,  gen.L, edesign)
  names(output) <- c("L","T", "data.L", "gen.L","edesign")
  output
}


#---------------------------------------------------------------------------------------------------
# Auxiliar internal functions: f3, MayorIso  
#---------------------------------------------------------------------------------------------------

f3 <- function(x)
{
  x<-x[!is.na(x)]
  y<-x[1]
  for (i in 2:length(x))
  {
    y<-paste(y,x[i],sep="&")
  }
  y
}
# ejemplo:
# x<-c(1,2,4,NA,NA)
# f3(x)


#-------------------------------------------------------

MayorIso<-function(zz){
  if( is.null(nrow(zz))){sol=1} else{
    M <-apply(zz,1,sum)  # sum expression across all the samples for each isoform
    sol=as.numeric(M==max(M)) # major isoform with  highest expression 
  }
  sol
}

#
# END maSigPro custom changes - should be incorporated into R package later
#

# run DIU Analysis based on selected method (fold filter already applied)
run_DIU_analysis <- function(degree, triple, transMatrix, myfactors){
  
  #degree = "1"
  #triple = FALSE
  design <- make.design.matrix(myfactors, degree)
  minobs = (as.numeric(degree) + 1) * groups
  result <- IsoModel_tappas(data=transMatrix, gen=genes[rownames(transMatrix),"id"], design=design, degree=degree,
                            min.obs = minobs, triple=triple) 
  
  # get the podium change information
  cmpType = "group"
  pcList = PodiumChange_tappas(iso = result, only.sig.genes = FALSE, comparison = "group", time.points=timepoints)
  
  
  # need to merge the two sets of data into a single dataframe
  result_diu <- NULL
  
  names(result)
  cat("\nlenGen: ", length(result$gen), ", lenDSG: ", length(result$DSG), ", lenadjPV: ", length(result$adjpvalue), "\nDSG:\n")
  # add tested genes, default to podium False
  result_sig <- data.frame("gene"=names(result$adjpvalue), "adjPValue"=result$adjpvalue, "podiumChange"=rep(FALSE, times=length(result$adjpvalue)), "podiumTimes"=rep('.', times=length(result$adjpvalue)), stringsAsFactors=FALSE)
  
  # add all remaining genes that have more than one transcript, default to podium False
  allGenes <- tapply(rownames(transMatrix), genes[rownames(transMatrix),"id"], length)
  multIsoGenes <- names(which(allGenes!=1))
  print(paste("Total multiple isoform genes = ", length(multIsoGenes)))
  #if there are genes not DS
  if(length(multIsoGenes[!(multIsoGenes %in% result_sig$gene)]!=0)){
    result_remaining <- data.frame("gene"=as.character(multIsoGenes[!(multIsoGenes %in% result_sig$gene)]), "adjPValue"=1.0, "podiumChange"=FALSE, "podiumTimes"='.', stringsAsFactors=FALSE)
    result_diu <- rbind(result_sig, result_remaining)
  }else{
    result_diu <- rbind(result_sig)
  }
  
  # set podium change flags
  result_diu[result_diu$gene %in% pcList$L, "podiumChange"] <- TRUE
  # set podium change time points if single series
  if(groups == 1) {
    # convert podium information from nested list to data frame
    podium_info <- data.frame("gene"=as.character(names(pcList$T[[1]])), "podiumTimes"=matrix(unlist(pcList$T[[1]]), nrow=length(pcList$T[[1]]), byrow=T), stringsAsFactors=FALSE)
    print(head(podium_info))
    # get index of rows to update in result_diu table - order in match call is important
    # we should never get an NA vector value in idx since all genes are in result_diu
    idx <- match(podium_info$gene, result_diu$gene)
    # if you need to remove at a later point use: idx <- idx[!is.na(idx)]
    print(paste("Number of podium gene matches:", length(idx)))
    # update podium times
    result_diu[idx, "podiumTimes"] <- podium_info["podiumTimes"]
    
    
    return(result_diu)
  }
  
  return(result_diu)
}



### Own functions
merge_exp <- function(input_result, gene_matrix){
  meangeneexp = data.frame(apply(gene_matrix,1,mean)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "mean_expression"))
  medgeneexp = data.frame(apply(gene_matrix,1,median)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "median_expression"))
  sumgeneexp = data.frame(apply(gene_matrix,1,sum)) %>% rownames_to_column(var = "gene") %>% `colnames<-`(c("gene", "sum_expression"))
  geneexp = merge(merge(meangeneexp,medgeneexp, by = "gene"), sumgeneexp, by = "gene")
  
  dat = merge(input_result, geneexp, by = "gene")
  return(dat)
}

DIU_time_analysis <- function(input_dir,degreevalue){
  # Whole Transcriptome: Iso-Seq scaffold + Iso-Seq reads
  gene_matrix=read.table(file.path(input_dir,"Data/gene_matrix.tsv")) 
  myfactors=read.table(file.path(input_dir, "Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
  groups = ncol(myfactors) - 2
  times = length(unique(myfactors[,1]))
  if(groups==1){timepoints = times}
  
  # proportion
  IsoIso_transMatrix_prop = DIU_prefilter(input_dir, mff=0.2, filteringType = "PROP")
  DIU_isoseq_prop = run_DIU_analysis(degree = degreevalue, triple = FALSE, IsoIso_transMatrix_prop,myfactors)
  DIU_isoseq_prop = merge_exp(DIU_isoseq_prop,gene_matrix)
  
  # fold change
  IsoIso_transMatrix_fc = DIU_prefilter(input_dir, mff=2.5, filteringType = "FOLD")
  DIU_isoseq_fc = run_DIU_analysis(degree = degreevalue, triple = FALSE, IsoIso_transMatrix_fc,myfactors)
  DIU_isoseq_fc = merge_exp(DIU_isoseq_fc,gene_matrix)
  
  output = list(DIU_isoseq_prop,DIU_isoseq_fc)
  names(output) = c("DIU_prop","DIU_fc")
  return(output)
}

DIU_time_analysis_filteronly <- function(input_dir, degreevalue){
  # Whole Transcriptome: Iso-Seq scaffold + Iso-Seq reads
  gene_matrix=read.table(file.path(input_dir,"Data/gene_matrix.tsv")) 
  myfactors=read.table(file.path(input_dir, "Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
  groups = ncol(myfactors) - 2
  times = length(unique(myfactors[,1]))
  if(groups==1){timepoints = times}
  
  # fold change
  IsoIso_transMatrix_fc = DIU_prefilter(input_dir, mff=2.5, filteringType = "FOLD")
}
