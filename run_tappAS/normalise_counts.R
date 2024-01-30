
# handle command line arguments - don't want to use external package so users don't have to install it
args = commandArgs(TRUE)
cat("Expression Matrix Filtering and Normalization script arguments: ", args, "\n")
srcfile = args[1]
expfile = args[2]
classfile= args[3]
dstfile = args[4]

suppressMessages(library("MASS"))
suppressMessages(library("NOISeq"))

tappas_filterNormalize <- function(data, length, factors, filterCPM, filterCOV) {
  cat("\nProcessing ", nrow(data), " transcript expression data rows")
  factors[,1] = as.factor(factors[,1])
  dataFiltered = filtered.data(data, factor = factors[,1], norm = FALSE, depth = NULL, method = 1, cv.cutoff = filterCOV, cpm = filterCPM, p.adj = "fdr")
  cat("\n", nrow(dataFiltered), " transcript expression data rows left after filtering")
  data_object = readData(data = dataFiltered, length = length, factors = factors)
  nlengths = as.vector(as.matrix(data_object@featureData@data))
  dataNormalized = tmm(assayData(data_object)$exprs, long = nlengths, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, k = 0, lc = 0)
  cat("\nReturning ", nrow(dataNormalized), " transcript expression data rows after filtering and normalization.")
  return(dataNormalized)
}

# read data files
#srcfile = "/lustre/home/sl693/AllMouseTargeted_sqantisubset.expression.txt"
#srcfile = "/lustre/home/sl693/input_matrix.tsv"
#expfile = "/lustre/home/sl693/exp_factors.txt"
classfile = read.table(classfile, header = T, row.names = 1, sep = "\t", as.is = T)


cat("\nReading input matrix file data...")
inputMatrix=read.table(srcfile, row.names=1, sep="\t", header=TRUE)
cat("\nRead ", nrow(inputMatrix), " transcripts expression data rows")
cat("\nReading factors file data...")
myfactors=read.table(expfile, row.names=1, sep="\t", header=TRUE)
cat("\nReading transcript length file data...")
#mylengths=read.table(translenfile, row.names=1, sep="\t", header=FALSE)
mylengths = subset(classfile, rownames(classfile) %in% row.names(inputMatrix), select = "length")
mylengths=as.vector(mylengths$length)

# process transcripts
ncpm=1.0
ncov=100.0
cat("\nFiltering, cutoff: ", ncpm, ", COV: ", ncov, ", and normalizing input matrix...\n")
results = tappas_filterNormalize(data=inputMatrix, length=mylengths, factors=myfactors, filterCPM=ncpm, filterCOV=ncov)
cat("\nSaving results to file...")
for(row in 1:nrow(results)) {
  for(col in 1:ncol(results))
    results[row, col] <- round(results[row, col], 2)
}
write.table(results, dstfile, sep="\t", row.names=TRUE, quote=FALSE)