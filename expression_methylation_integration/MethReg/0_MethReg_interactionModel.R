#' @title Fits linear models with interaction to triplet data (Target, TF, DNAm), where DNAm
#' is a binary variable (samples in Q1 or Q4)
#' @description Evaluates regulatory potential of DNA methylation (DNAm) on gene expression,
#' by fitting robust linear model or zero inflated negative binomial model to triplet data.
#' These models consist of terms to model direct effect of DNAm on target gene expression,
#' direct effect of TF on gene expression, as well as an interaction term that evaluates
#' the synergistic effect of DNAm and TF on gene expression.
#' @param triplet Data frame with columns for DNA methylation region (regionID), TF  (TF), and target gene  (target)
#' @param dnam DNA methylation matrix or SummarizedExperiment object
#' (columns: samples in the same order as \code{exp} matrix, rows: regions/probes)
#' @param exp A matrix or SummarizedExperiment object object
#'  (columns: samples in the same order as \code{dnam},
#' rows: genes represented by ensembl IDs (e.g. ENSG00000239415))
#' @param dnam.group.threshold DNA methylation threshold percentage to define samples 
#' in the low methylated group and high methylated group. For example, 
#' setting the threshold to 0.3 (30\%) will assign samples with the lowest 30\% 
#' methylation in the low group and the highest 30\% methylation in the high group. 
#' Default is 0.25 (25\%), accepted threshold range (0.0,0.5].
#' @param cores Number of CPU cores to be used. Default 1.
#' @param tf.activity.es A matrix with normalized enrichment scores for each TF across all samples
#' to be used in linear models instead of TF gene expression. See \code{\link{get_tf_ES}}.
#' @param sig.threshold Threshold to filter significant triplets.
#' Select if interaction.pval < 0.05 or pval.dnam < 0.05 or pval.tf < 0.05 in binary model
#' @param fdr Uses fdr when using sig.threshold.
#' Select if interaction.fdr < 0.05 or fdr.dnam < 0.05 or fdr.tf < 0.05 in binary model
#' @param filter.correlated.tf.exp.dnam  If wilcoxon test of TF expression Q1 and Q4 is significant (pvalue < 0.05),
#' triplet will be removed.
#' @param filter.correlated.target.exp.dnam  If wilcoxon test of target expression Q1 and Q4 is not significant (pvalue > 0.05),
#' triplet will be removed.
#' @param filter.triplet.by.sig.term Filter significant triplets ?
#' Select if interaction.pval < 0.05 or pval.dnam <0.05 or pval.tf < 0.05 in binary model
#' @param stage.wise.analysis A boolean indicating if stagewise analysis should be performed
#' to correct for multiple comparisons. If set to FALSE FDR analysis is performed.
#' @param verbose A logical argument indicating if
#' messages output should be provided.
#' @return A dataframe with \code{Region, TF, target, TF_symbo, target_symbol, estimates and P-values},
#' after fitting robust linear models or zero-inflated negative binomial models (see Details above).
#'
#' Model considering DNAm values as a binary variable generates \code{quant_pval_metGrp},
#' \code{quant_pval_rna.tf}, \code{quant_pval_metGrp.rna.tf},
#' \code{quant_estimates_metGrp}, \code{quant_estimates_rna.tf}, \code{quant_estimates_metGrp.rna.tf}.
#'
#' \code{Model.interaction} indicates which model (robust linear model or zero inflated model)
#' was used to fit Model 1, and \code{Model.quantile} indicates which model(robust linear model or zero
#' inflated model) was used to fit Model 2.
#'
#'@details This function fits the linear model
#'
#' \code{log2(RNA target) ~ log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' to triplet data as follow:
#'
#' Model by considering \code{DNAm} as a binary variable - we defined a binary group for
#' DNA methylation values (high = 1, low = 0). That is, samples with the highest
#' DNAm levels (top 25 percent) has high = 1, samples with lowest
#' DNAm levels (bottom 25 percent) has high = 0. Note that in this
#' implementation, only samples with DNAm values in the first and last quartiles
#' are considered.
#'
#' In these models, the term \code{log2(TF)} evaluates direct effect of TF on
#' target gene expression, \code{DNAm} evaluates direct effect of DNAm on target
#' gene expression, and \code{log2(TF)*DNAm} evaluates synergistic effect of DNAm
#' and TF, that is, if TF regulatory activity is modified by DNAm.
#'
#' There are two implementations of these models, depending on whether there are an excessive
#' amount (i.e. more than 25 percent) of samples with zero counts in RNAseq data:
#'
#' \itemize{
#' \item When percent of zeros in RNAseq data is less than
#' 25 percent, robust linear models are implemented using \code{rlm} function from \code{MASS} package. This
#' gives outlier gene expression values reduced weight. We used \code{"psi.bisqure"}
#' option in function \code{rlm} (bisquare weighting,
#' https://stats.idre.ucla.edu/r/dae/robust-regression/).
#'
#' \item When percent of zeros in RNAseq data is more than 25 percent, zero inflated negative binomial models
#' are implemented using \code{zeroinfl} function from \code{pscl} package. This assumes there are
#' two processes that generated zeros (1) one where the counts are always zero
#' (2) another where the count follows a negative binomial distribution.
#' }
#'
#' To account for confounding effects from covariate variables, first use the \code{get_residuals} function to obtain
#' RNA or DNAm residual values which have covariate effects removed, then fit interaction model. Note that no
#' log2 transformation is needed when \code{interaction_model} is applied to residuals data.
#'
#' Note that only triplets with TF expression not significantly different in high vs. low
#' methylation groups will be evaluated (Wilcoxon test, p > 0.05).
#'
#' @examples
#' library(dplyr)
#' dnam <- runif(20,min = 0,max = 1) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(dnam) <- c("chr3:203727581-203728580")
#' colnames(dnam) <- paste0("Samples",1:20)
#'
#' exp.target <-  runif(20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.target) <- c("ENSG00000252982")
#' colnames(exp.target) <- paste0("Samples",1:20)
#'
#' exp.tf <- runif(20,min = 0,max = 10) %>%
#'   matrix(ncol = 1) %>%  t
#' rownames(exp.tf) <- c("ENSG00000083937")
#' colnames(exp.tf) <- paste0("Samples",1:20)
#'
#' exp <- rbind(exp.tf, exp.target)
#'
#' triplet <- data.frame(
#'    "regionID" =  c("chr3:203727581-203728580"),
#'    "target" = "ENSG00000252982",
#'    "TF" = "ENSG00000083937"
#')
#' results <- interaction_model(
#'    triplet = triplet, 
#'    dnam = dnam, 
#'    exp = exp, 
#'     dnam.group.threshold = 0.25,
#'    stage.wise.analysis = FALSE, 
#'    sig.threshold = 1,
#'    filter.correlated.tf.exp.dnam = FALSE,
#'    filter.correlated.target.exp.dnam = FALSE,
#'    filter.triplet.by.sig.term = FALSE
#' )
#' @export
#' @importFrom rlang .data
#' @importFrom MASS rlm psi.bisquare
#' @importFrom plyr .
#' @importFrom stats wilcox.test
#' @importFrom dplyr group_by summarise filter_at contains vars any_vars pull filter
interaction_model <- function(
  triplet,
  dnam,
  exp,
  dnam.group.threshold = 0.25,
  cores = 1,
  tf.activity.es = NULL,
  sig.threshold = 0.05,
  fdr = TRUE,
  filter.correlated.tf.exp.dnam = TRUE,
  filter.correlated.target.exp.dnam = TRUE,
  filter.triplet.by.sig.term = TRUE,
  stage.wise.analysis = TRUE,
  verbose = FALSE
){
  
  if (stage.wise.analysis) library("stageR")
  
  if (missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
  if (missing(exp)) stop("Please set exp argument with gene expression matrix")
  
  if(!is(dnam.group.threshold,"numeric")) stop("dnam.group.threshold should be a value in the following interval (0,0.5]")
  if(dnam.group.threshold > 0.5) stop("dnam.group.threshold maximum valuee is 0.5")
  
  if (is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
  if (is(exp,"SummarizedExperiment")) exp <- assay(exp)
  
  #check_data(dnam, exp)
  
  if (!all(grepl("ENS", rownames(exp)))) {
    stop("exp must have the following row names as ENSEMBL IDs (i.e. ENSG00000239415)")
  }
  
  if (missing(triplet)) stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
  
  if (!all(c("regionID","TF","target") %in% colnames(triplet))) {
    stop("triplet must have the following columns names: regionID, TF, target")
  }
  
  if(verbose)  message("Removing genes with RNA expression equal to 0 for all samples from triplets")
  exp <- filter_genes_zero_expression_all_samples(exp)
  
  if(verbose)  message("Removing triplet with no DNA methylation information for more than 25% of the samples")
  regions.keep <- (rowSums(is.na(dnam)) < (ncol(dnam) * 0.75)) %>% which %>% names
  dnam <- dnam[regions.keep,,drop = FALSE]
  
  triplet <- triplet %>% dplyr::filter(
    .data$target %in% rownames(exp) &
      .data$regionID %in% rownames(dnam)
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
  
  # Remove cases where target is also the TF if it exists
  triplet <- triplet %>% dplyr::filter(
    .data$TF != .data$target
  )
  
  if (nrow(triplet) == 0) {
    stop("We were not able to find the same rows from triple in the data, please check the input.")
  }
  
  if(!"TF_symbol" %in% colnames(triplet))
    triplet$TF_symbol <- map_ensg_to_symbol(triplet$TF)
  
  if(!"target_symbol" %in% colnames(triplet))
    triplet$target_symbol <- map_ensg_to_symbol(triplet$target)
  
  if(!"target_region" %in% colnames(triplet))
    triplet$target_region <- map_ensg_to_region(triplet$target)
  
  if(!"distance_region_target_tss" %in% colnames(triplet)){
    triplet$distance_region_target_tss <- get_target_tss_to_region_distance(triplet$regionID,triplet$target)
  }
  
  if(verbose)  message("Evaluating ", nrow(triplet), " triplets")
  
  parallel <- register_cores(cores)
  
  ret <- plyr::adply(
    .data = triplet,
    .margins = 1,
    .fun = function(
      row.triplet
    ){
      
      data <- get_triplet_data(
        exp = exp,
        dnam = dnam,
        row.triplet = row.triplet,
        tf.es = tf.activity.es
      )
      
      upper.cutoff <-  quantile(data$met,na.rm = TRUE,  1 - dnam.group.threshold)
      low.cutoff <-  quantile(data$met,na.rm = TRUE,  dnam.group.threshold)
      
      quant.diff <- data.frame("met.IQR" = upper.cutoff - upper.cutoff)
      
      data.high.low <- data %>% filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
      data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)
      
      # pct.zeros.in.samples <- sum(data$rna.target == 0, na.rm = TRUE) / nrow(data)
      
      suppressWarnings({
        # Add information to filter TF if differenly expressed between DNAm high and DNAm low groups
        wilcoxon.tf.q4met.vs.q1met <- wilcox.test(
          data.high.low %>% filter(.data$metGrp == 1) %>% pull(.data$rna.tf),
          data.high.low %>% filter(.data$metGrp == 0) %>% pull(.data$rna.tf),
          exact = FALSE
        )$p.value
      })
      
      suppressWarnings({
        # Add information to filter Target if differently expressed between DNAm high and DNAm low groups
        wilcoxon.target.q4met.vs.q1met <- wilcox.test(
          data.high.low %>% filter(.data$metGrp == 1) %>% pull(.data$rna.target),
          data.high.low %>% filter(.data$metGrp == 0) %>% pull(.data$rna.target),
          exact = FALSE
        )$p.value
      })
      
      pct.zeros.in.quant.samples <- sum(
        data.high.low$rna.target == 0,
        na.rm = TRUE) / nrow(data.high.low)
      
      if (pct.zeros.in.quant.samples > 0.25) {
        itx.quant <- interaction_model_quant_zeroinfl(data.high.low)
      } else {
        itx.quant <- interaction_model_quant_rlm(data.high.low)
      }
      
      # Create output
      interaction_model_output(
        quant.diff,
        itx.quant,
        pct.zeros.in.quant.samples,
        wilcoxon.target.q4met.vs.q1met = wilcoxon.target.q4met.vs.q1met,
        wilcoxon.tf.q4met.vs.q1met = wilcoxon.tf.q4met.vs.q1met
      )
    },
    .progress = "time",
    .parallel = parallel,
    .inform = TRUE,
    .paropts = list(.errorhandling = 'pass')
  )
  
  if (stage.wise.analysis) {
    if(verbose)  message("Performing Stage wise correction for triplets")
    ret <- calculate_stage_wise_adjustment(ret)
  } else {
    if(verbose)  message("Performing FDR correction for triplets p-values per region")
    ret <- calculate_fdr_per_region_adjustment(ret)
  }
  
  if (filter.triplet.by.sig.term) {
    if(verbose)   message("Filtering results to have interaction, TF or DNAm significant")
    if(fdr){
      if (!stage.wise.analysis) pattern <- "fdr"
      if (stage.wise.analysis) pattern <- "triplet_stage_wise_adj_pvalue"
      
      ret <- ret %>% filter_at(vars(
        intersect(
          contains(pattern, ignore.case = TRUE),
          contains("RLM", ignore.case = TRUE)
        )
      ), 
      any_vars(. < sig.threshold)
      )
    } else { 
      ret <- ret %>% filter_at(vars(
        setdiff(
          intersect(
            contains("pvalue", ignore.case = TRUE),
            contains("RLM", ignore.case = TRUE)
          ),
          contains("stage_wise", ignore.case = TRUE)
        )
      ), 
      any_vars(. < sig.threshold)
      )
    }
  }
  
  if (filter.correlated.tf.exp.dnam) {
    if(verbose)  message("Filtering results to remove the significant in the wilcoxon test TF Q1 vs Q4")
    ret <- ret %>% 
      dplyr::filter(.data$TF_DNAm_high_vs_TF_DNAm_low_wilcoxon_pvalue > sig.threshold)
  }
  
  if (filter.correlated.target.exp.dnam) {
    if(verbose)  message("Filtering results to keep only the significant in the wilcoxon test target Q1 vs Q4")
    ret <- ret %>%
      dplyr::filter(.data$Target_gene_DNAm_high_vs_Target_gene_DNAm_low_wilcoxon_pvalue < sig.threshold)
  }
  
  # make the output more clear
  #colnames(ret) <- gsub(":","_x_",colnames(ret))
  colnames(ret) <- gsub("metGrp","DNAmGroup",colnames(ret))
  colnames(ret) <- gsub("rna.tf","TF",colnames(ret))
  # Since we used enrichment scores in the linear model
  # we will rename the output
  if (!is.null(tf.activity.es)) {
    colnames(ret) <- gsub("rna.tf","TF.es",colnames(ret))
  }
  colnames(ret) <- gsub("rna.tf","TF",colnames(ret))
  
  ret
}


interaction_all_model_no_results <- function(){
  cbind(
    "pval_met" = NA,
    "pval_rna.tf" = NA,
    "pval_met:rna.tf" = NA,
    "estimate_met" = NA,
    "estimate_rna.tf" = NA,
    "estimate_met:rna.tf" = NA
  ) %>% as.data.frame()
}

interaction_quant_model_no_results <- function(){
  cbind(
    "RLM_metGrp_pvalue" = NA,
    "RLM_rna.tf_pvalue" = NA,
    "RLM_metGrp:rna.tf_pvalue" = NA,
    "RLM_metGrp_estimate" = NA,
    "RLM_rna.tf_estimate" = NA,
    "RLM_metGrp:rna.tf_estimate" = NA
  ) %>% as.data.frame()
}

get_triplet_data <- function(
  exp,
  dnam,
  row.triplet,
  tf.es,
  add.groups = FALSE
){
  rna.target <- exp[which(rownames(exp) == row.triplet$target), , drop = FALSE]
  met <- dnam[which(rownames(dnam) == as.character(row.triplet$regionID)), , drop = FALSE]
  
  if (!is.null(tf.es)) {
    rna.tf <- tf.es[which(rownames(tf.es) == row.triplet$TF), , drop = FALSE]
  } else {
    rna.tf <- exp[which(rownames(exp) == row.triplet$TF), , drop = FALSE]
  }
  
  data <- data.frame(
    rna.target = rna.target %>% as.numeric,
    met = met %>% as.numeric,
    rna.tf = rna.tf %>% as.numeric
  )
  
  if(add.groups){
    quant.met <-  quantile(data$met,na.rm = TRUE)
    low.cutoff <- quant.met[2]
    upper.cutoff <- quant.met[4]
    data$metGrp <- NA
    data$metGrp[data$met <= low.cutoff] <- "DNAm low"
    data$metGrp[data$met >= upper.cutoff] <- "DNAm high"
    
    quant.rna.tf <-  quantile(data$rna.tf,na.rm = TRUE)
    low.cutoff <- quant.rna.tf[2]
    upper.cutoff <- quant.rna.tf[4]
    data$TF.group <- NA
    data$TF.group[data$rna.tf <= low.cutoff] <- "TF low"
    data$TF.group[data$rna.tf >= upper.cutoff] <- "TF high"
  }
  data
}

interaction_model_output <- function(
  # itx.all,
  #pct.zeros.in.samples,
  quant.diff,
  itx.quant,
  pct.zeros.in.quant.samples,
  wilcoxon.target.q4met.vs.q1met,
  wilcoxon.tf.q4met.vs.q1met
){
  if (is.null(itx.quant)) itx.quant <- interaction_quant_model_no_results()
  
  cbind(
    quant.diff,
    itx.quant,
    data.frame(
      "Model quantile" =
        ifelse(pct.zeros.in.quant.samples > 0.25,
               "Zero-inflated Negative Binomial Model",
               "Robust Linear Model"
        ),
      "Target_gene_DNAm_high_vs_Target_gene_DNAm_low_wilcoxon_pvalue" = wilcoxon.target.q4met.vs.q1met,
      "TF_DNAm_high_vs_TF_DNAm_low_wilcoxon_pvalue" = wilcoxon.tf.q4met.vs.q1met
    ),
    "% of target genes not expressed in DNAm_low and DNAm_high" = paste0(round(pct.zeros.in.quant.samples * 100,digits = 2)," %")
  )
}


interaction_model_rlm <- function(data){
  rlm.bisquare <- tryCatch({
    # 2) fit linear model: target RNA ~ DNAm + RNA TF
    rlm(
      rna.target ~ met + rna.tf + rna.tf * met,
      data = data,
      psi = MASS::psi.bisquare,
      maxit = 100
    ) %>% summary %>% coef %>% data.frame
  }, error = function(e) {
    # message("Continuous model: ", e)
    return(NULL)
  })
  
  # if (is.null(rlm.bisquare)) return(interaction_model_no_results())
  if (is.null(rlm.bisquare)) return(interaction_all_model_no_results())
  # if (is.null(rlm.bisquare)) return(NULL)
  
  if (!"met:rna.tf" %in% rownames(rlm.bisquare)) {
    rlm.bisquare <- rbind(
      rlm.bisquare,
      data.frame(
        row.names = "met:rna.tf",
        "Value" = NA,
        "Std..Error" = NA,
        "t.value" = NA
      )
    )
  }
  
  degrees.freedom.value <- nrow(data) - 4
  rlm.bisquare$pval <- 2 * (1 - pt( abs(rlm.bisquare$t.value), df = degrees.freedom.value) )
  
  all.pval <- rlm.bisquare[-1,4,drop = FALSE] %>% t %>% as.data.frame()
  colnames(all.pval) <- paste0(colnames(all.pval),"_pvalue")
  
  all.estimate <- rlm.bisquare[-1,1,drop = FALSE] %>% t %>% as.data.frame()
  colnames(all.estimate) <- paste0(colnames(all.estimate),"_estimate")
  return(cbind(all.pval, all.estimate))
}

#' @importFrom pscl zeroinfl
interaction_model_zeroinfl <- function(data){
  zinb <- tryCatch({
    suppressWarnings({
      pscl::zeroinfl(
        trunc(rna.target) ~ met + rna.tf + rna.tf * met | 1,
        data = data,
        dist = "negbin",
        EM = FALSE) %>% summary %>% coef
    })
  }, error = function(e) {
    # message("Continuous model: ", e)
    return(NULL)
  })
  if (is.null(zinb)) return(interaction_all_model_no_results())
  
  zinb <- zinb$count %>% data.frame
  
  all.pval <- zinb[c(-1,-5),4,drop = FALSE] %>% t %>% as.data.frame()
  colnames(all.pval) <- paste0(colnames(all.pval), "_pvalue")
  
  all.estimate <- zinb[c(-1,-5),1,drop = FALSE] %>% t %>% as.data.frame()
  colnames(all.estimate) <- paste0(colnames(all.estimate),"_estimate")
  return(cbind(all.pval, all.estimate))
}

interaction_model_quant_zeroinfl <- function(data){
  zinb.quant <- tryCatch({
    suppressWarnings({
      pscl::zeroinfl(
        trunc(rna.target) ~ metGrp + rna.tf + metGrp * rna.tf | 1,
        data = data,
        dist = "negbin",
        EM = FALSE
      ) %>% summary %>% coef
    })
  }, error = function(e) {
    # message("Continuous model: ", e)
    return(NULL)
  })
  if (is.null(zinb.quant)) return(interaction_quant_model_no_results())
  
  zinb.quant <- zinb.quant$count %>% data.frame
  quant.pval <- zinb.quant[c(-1,-5),4,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")
  
  quant.estimate <- zinb.quant[c(-1,-5),1,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")
  
  return(cbind(quant.pval, quant.estimate))
}


interaction_model_quant_rlm <- function(data){
  rlm.bisquare.quant <- tryCatch({
    suppressWarnings({
      rlm(
        rna.target ~ metGrp + rna.tf + metGrp * rna.tf,
        data = data,
        psi = MASS::psi.bisquare,
        maxit = 100
      ) %>% summary %>% coef %>% data.frame
    })
  }, error = function(e) {
    #message("Binary model: ", e)
    return(NULL)
  })
  
  if (is.null(rlm.bisquare.quant)) return(interaction_quant_model_no_results())
  
  # if the interaction is NA, it is removed from the data frame,
  # we have to re add it
  if (!"metGrp:rna.tf" %in% rownames(rlm.bisquare.quant)) {
    rlm.bisquare.quant <- rbind(
      rlm.bisquare.quant,
      data.frame(
        row.names = "metGrp:rna.tf",
        "Value" = NA,
        "Std..Error" = NA,
        "t.value" = NA
      )
    )
  }
  
  degrees.freedom.value <- nrow(data) - 4
  rlm.bisquare.quant$pval <- 2 * (1 - pt( abs(rlm.bisquare.quant$t.value),
                                          df = degrees.freedom.value) )
  
  quant.pval <- rlm.bisquare.quant[-1,4,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")
  
  quant.estimate <- rlm.bisquare.quant[-1,1,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")
  
  return(cbind(quant.pval, quant.estimate))
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
  
  interactiol.col <- "RLM_metGrp.rna.tf_pvalue"
  results <- stage_wise_adjustment(results, interactiol.col)
  
  dnam.col <- grep("metGrp_pvalue",colnames(results),value = TRUE)
  results <- stage_wise_adjustment(results, dnam.col)
  
  tf.col <- grep("tf_pvalue$",colnames(results),value = TRUE)
  tf.col <- grep("metGrp",tf.col,invert = TRUE,value = TRUE)
  results <- stage_wise_adjustment(results, tf.col)
  
  return(results)
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
#' results <- stage_wise_adjustment(results,"quant_pval_metGrp.rna.tf")
#' @noRd
stage_wise_adjustment <- function(
  results,
  col
){
  
  if(!"tripletID" %in% colnames(results)){
    results$tripletID <- create_triplet_ID(results)
  }
  
  
  min.pval <- results %>%
    dplyr::group_by(.data$regionID) %>%
    dplyr::summarise(min(.data[[col]]))
  
  # Preparing StageR input
  pScreen.pval <- min.pval[[2]]
  names(pScreen.pval) <- gsub("[[:punct:]]", "_", min.pval$regionID)
  
  pConfirmation <- results[,col] %>% as.matrix  ## LW: change to p-values
  rownames(pConfirmation) <-  results$tripletID
  colnames(pConfirmation) <- "transcript"
  
  triplet2region <- data.frame(
    row.names = results$tripletID,
    "transcript" = results$tripletID,
    "gene" = gsub("[[:punct:]]", "_", results$regionID)
  )
  
  if(all(is.na(pScreen.pval))) {
    stop("Stage wise error. Set it to FALSE")
    return(results)
  }
  pScreen.pval.stageRObj <- stageR::stageRTx(
    pScreen = pScreen.pval,
    pConfirmation = pConfirmation,
    pScreenAdjusted = FALSE,
    tx2gene = triplet2region
  )
  
  pScreen.pval.stageRObj <- stageR::stageWiseAdjustment(
    object = pScreen.pval.stageRObj,
    method = "dte",
    alpha = 0.05,
    allowNA = TRUE
  )
  
  padj <- stageR::getAdjustedPValues(
    pScreen.pval.stageRObj,
    onlySignificantGenes = FALSE,
    order = FALSE
  )
  # padj: geneID", "txID" ,"gene" ,"transcript"
  # equivalento to: region, triplet, region pval, triplet.pval
  
  results[[gsub("pvalue", "region_stage_wise_adj_pvalue", col)]] <- padj[[3]][match(results$tripletID, padj[[2]])]
  results[[gsub("pvalue", "triplet_stage_wise_adj_pvalue", col)]] <- padj[[4]][match(results$tripletID, padj[[2]])]
  results$tripletID <- NULL
  results <- results %>%
    relocate(
      c(
        gsub("pvalue","region_stage_wise_adj_pvalue",col),
        gsub("pvalue","triplet_stage_wise_adj_pvalue",col)
      ),
      .after = col
    )
  
  return(results)
}

#' @examples
#' \dontrun{
#' data("dna.met.chr21")
#' dna.met.chr21 <- make_se_from_dnam_probes(dna.met.chr21)
#' data("gene.exp.chr21.log2")
#' triplet <- data.frame(
#'     "regionID" = rownames(dna.met.chr21)[1:10],
#'     "TF" = rownames(gene.exp.chr21.log2)[11:20],
#'     "target" = rownames(gene.exp.chr21.log2)[1:10]
#' )
#' triplet$tripletID <- create_triplet_ID(triplet)
#' }
#' @noRd
create_triplet_ID <- function(df){
  paste0(gsub("[[:punct:]]", "_", df$regionID),"_TF_",df$TF,"_target_",df$target)
}

#' @title
#' Remove genes with gene expression level equal to 0 or NA in a all samples
#' @param exp Gene expression matrix or a Summarized Experiment object
#' @noRd
#' @examples
#' data("gene.exp.chr21.log2")
#' gene.exp.chr21.log2.filtered <- filter_genes_zero_expression_all_samples(
#'   gene.exp.chr21.log2
#' )
#' @return
#' A subset of the original matrix only with the rows
#' passing the filter threshold.
filter_genes_zero_expression_all_samples <- function(
  exp
){
  if(is(exp,"SummarizedExperiment")){
    exp <- assay(exp)
  }
  idx.all.zero <- rowSums(exp == 0, na.rm = TRUE) == ncol(exp)
  idx.all.na <- rowSums(is.na(exp)) == ncol(exp)
  
  # do not keep if it is all zero or all NA
  genes.keep <- rownames(exp)[!(idx.all.zero | idx.all.na)] %>% na.omit()
  if(length(genes.keep) < nrow(exp) & length(genes.keep) > 0){
    message(
      "Removing ", nrow(exp) - length(genes.keep),
      " out of ", nrow(exp), " genes"
    )
  }
  exp <- exp[genes.keep,,drop = FALSE]
  return(exp)
}
