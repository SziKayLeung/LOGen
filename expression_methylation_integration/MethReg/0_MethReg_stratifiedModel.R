#' @title Fits linear models to triplet data (Target, TF, DNAm) for
#' samples with high DNAm or low DNAm separately, and annotates TF
#' (activator/repressor) and DNam effect over TF activity (attenuate, enhance).
#' @description Should be used after fitting \code{interaction_model}, and only
#' for triplet data with significant \code{TF*DNAm} interaction. This analysis
#' examines in more details on how TF activities differ in
#' samples with high DNAm or low DNAm values.
#' @param triplet Data frame with columns for
#' DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix or SummarizedExperiment
#' (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A matrix or SummarizedExperiment
#' (columns: samples in the same order as \code{dnam} matrix,
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @param cores Number of CPU cores to be used. Default 1.
#' @param tf.activity.es A matrix with normalized enrichment scores
#' for each TF across all samples
#' to be used in linear models instead of TF gene expression.
#' @param tf.dnam.classifier.pval.thld P-value threshold to consider
#' a linear model significant
#' of not. Default 0.001. This will be used to classify the TF role and DNAm
#' effect.
#' @param dnam.group.threshold DNA methylation threshold percentage to define samples 
#' in the low methylated group and high methylated group. For example, 
#' setting the threshold to 0.3 (30\%) will assign samples with the lowest 30\% 
#' methylation in the low group and the highest 30\% methylation in the high group. 
#' Default is 0.25 (25\%), accepted threshold range (0.0,0.5].
#' @return A data frame with \code{Region, TF, target, TF_symbol target_symbol},
#' results for
#' fitting linear models to samples with low methylation
#'  (\code{DNAmlow_pval_rna.tf}, \code{DNAmlow_estimate_rna.tf}),
#'  or samples with high methylation (\code{DNAmhigh_pval_rna.tf},
#' \code{DNAmhigh_pval_rna.tf.1}), annotations for TF (\code{class.TF})
#' and (\code{class.TF.DNAm}).
#'
#' @details This function fits linear model
#' \code{log2(RNA target) = log2(TF)}
#'
#' to samples with highest DNAm values (top 25 percent) or
#' lowest DNAm values (bottom 25 percent), separately.
#'
#' There are two implementations of these models, depending on whether there are an excessive
#' amount (i.e. more than 25 percent) of samples with zero counts in RNAseq data:
#'
#' \itemize{
#' \item When percent of zeros in RNAseq data is less than
#' 25 percent, robust linear models are implemented using \code{rlm}
#' function from \code{MASS} package. This
#' gives outlier gene expression values reduced weight. We used \code{"psi.bisqure"}
#' option in function \code{rlm} (bisquare weighting,
#' https://stats.idre.ucla.edu/r/dae/robust-regression/).
#'
#' \item When percent of zeros in RNAseq data is more than 25 percent,
#' zero inflated negative binomial models
#' are implemented using \code{zeroinfl} function from \code{pscl} package. This assumes there are
#' two processes that generated zeros (1) one where the counts are always zero
#' (2) another where the count follows a negative binomial distribution.
#'}
#'
#' To account for confounding effects from covariate variables,
#' first use the \code{get_residuals} function to obtain
#' RNA residual values which have covariate effects removed,
#' then fit interaction model. Note that no
#' log2 transformation is needed when \code{interaction_model}
#' is applied to residuals data.
#'
#' This function also provides annotations for TFs. A TF is annotated as
#' \code{activator} if
#' increasing amount of TF (higher TF gene expression) corresponds to
#' increased target gene expression. A TF
#' is annotated as \code{repressor} if increasing amount of TF
#' (higher TF gene expression) corresponds to
#' decrease in target gene expression.
#' A TF is annotated as \code{dual} if in the Q1 methylation group increasing
#' amount of TF (higher TF gene expression) corresponds to
#' increase in target gene expression, while in Q4 methylation group increasing
#' amount of TF (higher TF gene expression) corresponds to
#' decrease in target gene expression
#' (or the same but changing Q1 and Q4 in the previous sentence).
#'
#' In addition, a region/CpG is annotated as \code{enhancing} if more
#' TF regulation on gene transcription
#' is observed in samples with high DNAm. That is,  DNA methylation
#' enhances TF regulation on target gene expression.
#' On the other hand, a region/CpG is annotated as \code{attenuating}
#'  if more TF regulation on gene
#' transcription is observed in samples with low DNAm.
#' That is, DNA methylation reduces TF regulation
#' on target gene expression.
#'
#' @examples
#' library(dplyr)
#' dnam <- runif (20,min = 0,max = 1) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'
#' exp.target <-  runif (20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.target) <- c("ENSG00000232886")
#' colnames(exp.target) <- paste0("Samples",1:20)
#'
#' exp.tf <- runif (20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.tf) <- c("ENSG00000232888")
#' colnames(exp.tf) <- paste0("Samples",1:20)
#'
#' exp <- rbind(exp.tf, exp.target)
#'
#' triplet <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000232886",
#'    "TF" = "ENSG00000232888"
#')
#'
#' results <- stratified_model(
#'   triplet = triplet,
#'   dnam = dnam,
#'   exp = exp
#' )
#' @export
#' @importFrom tibble tibble
#' @importFrom rlang .data
stratified_model <- function(
  triplet,
  dnam,
  exp,
  cores = 1,
  tf.activity.es = NULL,
  tf.dnam.classifier.pval.thld = 0.001,
  dnam.group.threshold = 0.25
){
  
  if (missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
  if (missing(exp)) stop("Please set exp argument with gene expression matrix")
  
  if (is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
  if (is(exp,"SummarizedExperiment")) exp <- assay(exp)
  
  if (!all(grepl("ENSG", rownames(exp)))) {
    stop("exp must have the following row names as ENSEMBL IDs (i.e. ENSG00000239415)")
  }
  
  if (missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
  if (!all(c("regionID","TF","target") %in% colnames(triplet))) {
    stop("triplet must have the following columns names: regionID, TF, target")
  }
  
  message("Removing genes with RNA expression equal to 0 for all samples from triplets")
  exp <- filter_genes_zero_expression_all_samples(exp)
  
  message("Removing triplet with no DNA methylation information for more than 25% of the samples")
  regions.keep <- (rowSums(is.na(dnam)) < (ncol(dnam) * 0.75)) %>% which %>% names
  dnam <- dnam[regions.keep,,drop = FALSE]
  
  triplet <- triplet %>% dplyr::filter(
    .data$target %in% rownames(exp) &
      .data$regionID %in% rownames(dnam)
  )
  
  if(!"TF_symbol" %in% colnames(triplet)){
    triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)
  }
  
  if(!"target_symbol" %in% colnames(triplet)){
    triplet$target_symbol <- map_ensg_to_symbol(triplet$target)
  }
  
  # Remove cases where target is also the TF if it exists
  triplet <- triplet %>% dplyr::filter(
    .data$TF != .data$target
  )
  
  if (!is.null(tf.activity.es)) {
    
    if(any(is.na(rownames(tf.activity.es))))
      tf.activity.es <- tf.activity.es[!is.na(rownames(tf.activity.es)),]
    
    if (!all(grepl("^ENSG", rownames(tf.activity.es)))) {
      rownames(tf.activity.es) <- map_symbol_to_ensg(rownames(tf.activity.es))
    }
    
    triplet <- triplet %>% dplyr::filter(
      .data$TF %in% rownames(tf.activity.es)
    )
  } else {
    triplet <- triplet %>% dplyr::filter(
      .data$TF %in% rownames(exp)
    )
  }
  
  if (nrow(triplet) == 0) {
    stop("We were not able to find the same rows from triple in the data, please check the input.")
  }
  
  parallel <- register_cores(cores)
  
  out <- plyr::adply(
    .data = triplet,
    .margins = 1,
    .fun = function(row.triplet){
      
      data <- get_triplet_data(
        exp = exp,
        dnam = dnam,
        row.triplet = row.triplet,
        tf.es = tf.activity.es
      )
      interaction.significant <- TRUE
      if("RLM_DNAmGroup:TF_pvalue" %in% colnames(row.triplet)){
        interaction.significant <- ifelse(row.triplet[,"RLM_DNAmGroup:TF_pvalue"] < 0.05,TRUE,FALSE)
        if(is.na(interaction.significant)) interaction.significant <- FALSE
      }
      
      stratified_model_results(
        data = data, 
        tf.dnam.classifier.pval.thld = tf.dnam.classifier.pval.thld, 
        dnam.group.threshold = dnam.group.threshold,
        interaction.significant = interaction.significant
      )
    }, .progress = "time", .parallel = parallel, .inform = TRUE)
  
  if (!is.null(tf.activity.es)) {
    colnames(out) <- gsub("rna.tf","es.tf",colnames(out))
  }
  
  return(out)
}

stratified_model_results <- function(
  data,
  tf.dnam.classifier.pval.thld = 0.001,
  dnam.group.threshold = 0.25,
  interaction.significant = FALSE
){
  upper.cutoff <-  quantile(data$met,na.rm = TRUE,  1 - dnam.group.threshold)
  low.cutoff <-  quantile(data$met,na.rm = TRUE,  dnam.group.threshold)
  
  data.low <- data %>% dplyr::filter(.data$met <= low.cutoff)
  data.high <- data %>% dplyr::filter(.data$met >= upper.cutoff)
  
  results.low <- stratified_model_aux(data.low,"DNAmlow")
  results.low.pval <- results.low$pval
  results.low.estimate <- results.low$estimate
  
  results.high <- stratified_model_aux(data.high,"DNAmhigh")
  results.high.pval <- results.high$pval
  results.high.estimate <- results.high$estimate
  
  
  classification <- get_tf_dnam_classification(
    low.estimate = results.low.estimate,
    low.pval = results.low.pval,
    high.estimate = results.high.estimate,
    high.pval = results.high.pval,
    pvalue.threshold = tf.dnam.classifier.pval.thld
  )
  
  if(!interaction.significant){
    classification$DNAm.effect <- "ns"
  }
  
  tibble::tibble(
    "DNAm_low_RLM_target_vs_TF_pvalue" = results.low.pval %>% as.numeric(),
    "DNAm_low_RLM_target_vs_TF_estimate" = results.low.estimate %>% as.numeric(),
    "DNAm_high_RLM_target_vs_TF_pvalue" = results.high.pval %>% as.numeric(),
    "DNAm_high_RLM_target_vs_TF_estimate" = results.high.estimate %>% as.numeric(),
    "DNAm.effect" = classification$DNAm.effect,
    "TF.role" = classification$TF.role
  )
}

#' @importFrom MASS rlm psi.bisquare
#' @importFrom stats coef pt
stratified_model_aux <- function(data, prefix = ""){
  pct.zeros.samples <- sum(data$rna.target == 0, na.rm = TRUE) / nrow(data)
  
  if (pct.zeros.samples > 0.25) {
    results <-  tryCatch({
      pscl::zeroinfl(
        trunc(rna.target) ~ rna.tf | 1,
        data = data,
        dist = "negbin",
        EM = FALSE) %>% summary %>% coef
    }, error = function(e){
      # message("Binary model: ", e)
      return(NULL)
    })
    
    if (is.null(results)) return(stratified_model_aux_no_results(pct.zeros.samples))
    
    results <- results$count %>% data.frame
    
    results.pval <- results["rna.tf","Pr...z..",drop = FALSE] %>%
      t %>%
      as.data.frame()
    colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))
    
    results.estimate <- results["rna.tf","Estimate",drop = FALSE] %>%
      t %>%
      as.data.frame()
    colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))
    
  } else {
    
    results <- tryCatch({
      
      MASS::rlm(
        formula = as.formula("rna.target ~ rna.tf"),
        data = data,
        psi = psi.bisquare,
        maxit = 100) %>% summary %>% coef %>% data.frame
      
    }, error = function(e){
      # message("Binary model: ", e)
      return(NULL)
    })
    
    if (is.null(results)) return(stratified_model_aux_no_results(pct.zeros.samples))
    
    degrees.freedom.value <- nrow(data) - 2
    results$pval <- 2 * (1 - pt( abs(results$t.value), df = degrees.freedom.value) )
    
    results.pval <- results[-1,4,drop = FALSE] %>% t %>% as.data.frame()
    colnames(results.pval) <- paste0(prefix,"_pval_",colnames(results.pval))
    
    results.estimate <- results[-1,1,drop = FALSE] %>% t %>% as.data.frame()
    colnames(results.estimate) <- paste0(prefix,"_estimate_",colnames(results.estimate))
  }
  
  return(
    list(
      "estimate" = results.estimate,
      "pval" = results.pval,
      "Model" = ifelse(pct.zeros.samples > 0.25,
                       "Zero-inflated Negative Binomial Model",
                       "Robust Linear Model"),
      "percet_zero_target_genes" = paste0(round(pct.zeros.samples * 100, digits = 2)," %")
    )
  )
}
stratified_model_aux_no_results <- function(pct.zeros.samples){
  list("estimate" = NA,
       "pval" = NA,
       "Model" = "Robust Linear Model",
       "percent_zero_target_genes" = paste0(round(pct.zeros.samples * 100, digits = 2)," %")
  )
  
}

#' @title TF and DNAm roles classifier
#' @description
#' This function receives the pvalue and estimate
#' for the following linear models:
#' 1) Target ~ TF for  DNAm low  (Q1)
#' 2) Target ~ TF for  DNAm high (Q4)
#' Then it will classify the TF (TF.role) as:
#' - Repressor (Increase of TF decreases Target)
#' - Activator (Increase of TF increases Target)
#' - Invert (Increase of TF decreases Target for Q1 and
#' Increase of TF increases Target for Q1; or
#' Increase of TF decreases Target for Q4
#' )
#' And classify the DNA methylation effect (dnam.effect) as:
#' - Attenuating (Attenuates TF activity - Estimate Q4 < Estimate Q1)
#' - Enhancing (Enhances TF activity - Estimate Q4 > Estimate Q1)
#' @noRd
#' @examples
#' get_tf_dnam_classification(
#'   low.estimate = 0.2, low.pval = 0.05,
#'   high.estimate = 0.8, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = 0.8, low.pval = 0.05,
#'   high.estimate = 0.2, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = -0.8, low.pval = 0.05,
#'   high.estimate = -0.2, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = -0.2, low.pval = 0.05,
#'   high.estimate = -0.8, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = -0.8, low.pval = 0.05,
#'   high.estimate = 0.8, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = 0.8, low.pval = 0.05,
#'   high.estimate = -0.8, high.pval = 0.05
#' )
#' get_tf_dnam_classification(
#'   low.estimate = 0.8, low.pval = 1,
#'   high.estimate = -0.2, high.pval = 1
#' )
#' get_tf_dnam_classification(
#'   low.estimate = 0.8, low.pval = 1,
#'   high.estimate = 0.2, high.pval = 0.05
#' )
get_tf_dnam_classification <- function(
  low.estimate,
  low.pval,
  high.estimate,
  high.pval,
  pvalue.threshold = 0.001
){
  
  # output
  classification <- list("DNAm.effect" = NA, "TF.role" = NA)
  
  pval.vct <- c(low.pval %>% as.numeric, high.pval %>% as.numeric)
  pval.sig <- pval.vct <= pvalue.threshold
  pval.sig[is.na(pval.sig)] <- FALSE
  estimate.vector <- c(low.estimate %>% as.numeric, high.estimate %>% as.numeric)
  
  if (any(is.na(estimate.vector))) {
    return(classification)
  }
  
  # All estimates are not significant
  if (!any(pval.sig,na.rm = TRUE)) {
    return(classification)
  }
  
  
  # TF role classification
  # Activator
  # Repressor
  if (all(estimate.vector[pval.sig] > 0)) {
    classification$TF.role <- "Activator"
  } else if (all(estimate.vector[pval.sig] < 0)) {
    classification$TF.role <- "Repressor"
  } else {
    # One is + and the other -
    classification$TF.role <- "Dual"
  }
  
  if (is.na(low.pval)) {
    classification$DNAm.effect <- "Enhancing"
  } else if (is.na(high.pval)) {
    classification$DNAm.effect <- "Attenuating"
  } else if (low.pval <= pvalue.threshold & high.pval >= pvalue.threshold) {
    classification$DNAm.effect <- "Attenuating"
  } else if (high.pval <= pvalue.threshold & low.pval >= pvalue.threshold  ) {
    classification$DNAm.effect <- "Enhancing"
  } else if (high.pval <= pvalue.threshold &  low.pval <= pvalue.threshold & abs(low.estimate) > abs(high.estimate) ) {
    # in this case both are significant
    classification$DNAm.effect <- "Attenuating"
  } else if (high.pval <= pvalue.threshold &  low.pval <= pvalue.threshold & abs(low.estimate) < abs(high.estimate) ) {
    # in this case both are significant
    classification$DNAm.effect <- "Enhancing"
  }
  
  if (classification$TF.role == "Dual") {
    classification$DNAm.effect <- "Invert"
  }
  
  return(classification)
}
