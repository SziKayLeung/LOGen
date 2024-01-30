# Szi Kay Leung: sl693@exeter.ac.uk
# Differntial Expression Analaysis script from tappAS

# run time course DEA analysis
DEA_time_analysis <- function(dtname, degree, siglevel, r2cutoff, knum, usemclust, groups, indir) {
  # read expression factors
  cat("\nReading factors file data...")
  factors=read.table(file.path(indir, "/Data/time_factors.txt"), row.names=1, sep="\t", header=TRUE)
  groups = ncol(factors) - 2
  times = length(unique(factors[,1]))
  
  # read expression matrix
  cat("\nReading normalized ", dtname, " matrix file data...")
  expMatrix = read.table(file.path(indir, paste0("/Data/",dtname, "_matrix.tsv")), row.names=1, sep="\t", header=TRUE)
  cat("\nRead ", nrow(expMatrix), " normalized ", dtname, " expression data rows")
  
  # create a regression matrix for regression model 
  design <- make.design.matrix(factors, degree)
  minobs = (degree + 1) * groups
  #cat("\ngroups: ", groups, ", degree: ", degree, ", minobs: ", minobs)
  #cat("\nsiglevel: ", siglevel, ", r2cutoff: ", r2cutoff, ", knum: ", knum, ", useMclust: ", usemclust)
  #print(design$edesign)
  #cat("\nRunning p.vector()...\n")
  
  ### Finding significant genes
  # regression fit for each gene for each condition 
  # p value for each regression fit associated with F statistic after correction with BH FDR
  # counts = TRUE uses a GLM and applies NB distribution
  # A gene with different profile (regression fit) between 2 conditions will have statistically significant coefficient
  fit <- p.vector(expMatrix, design, Q=siglevel, MT.adjust="BH", min.obs=minobs, counts=TRUE)
  cat("\nFound ", fit$i, " out of ", fit$g, " items to be significant\n")
  #cat("\ndesign:\n")
  print(head(fit$SELEC))
  #fit$p.vector
  
  ### Finding significant differences
  # stepwise regression - iterative regression approach that selects from a pool of potential variables the ‘best’ ones (according to a specified criterion) to fit the available data.
  # Apply a threshold of Rsqaured coefficient cutoff = goodness of fit for meaningful differences
  
  cat("\nRunning T.fit() using step.method='backward' and alfa = ", siglevel, "...\n")
  tstep <- T.fit(fit, step.method="backward", alfa=siglevel)
  tstep$sol
  
  cat("\nRunning get.siggenes() using rsq=", r2cutoff, " and vars='groups'...\n")
  sigs <- get.siggenes(tstep, rsq=r2cutoff, vars="groups")
  print("YESSSSSSS")

  for(i in 1:length(names(sigs$sig.genes))) {
    gname <- names(sigs$sig.genes)[i]
    cat("\n", gname, " significant ", dtname, ": ")
    print(sigs$sig.genes[[gname]]$g)
    df <- data.frame(sigs$sig.genes[[gname]]$sig.pvalues[, 1:2])
    df <- df[order(df$p.value), ]
    #df <<- df
    print(nrow(df))
    print(head(df))
    cat("\nSaving ", dtname, " DEA results to file...")
  }
  return(tstep$sol)
}
