#' @title Plot interaction model results
#' @description Create several plots to show interaction data
#' TF expression with target gene interaction using a linear model
#' \deqn{log2(RNA target) = log2(TF) + DNAm + log2(TF) * DNAm}
#'
#' To consider covariates, RNA can also be the residuals.
#' \deqn{log2(RNA target residuals) = log2(TF residual) + DNAm + log2(TF residual) * DNAm}
#'
#' @param triplet.results Output from function interaction_model
#' with Region ID, TF  (column name: TF),  and target gene  (column name: target),
#' p-values and estimates of interaction
#' @param dnam DNA methylation matrix or SummarizedExperiment object
#'  (columns: samples same order as met, rows: regions/probes)
#' @param exp gene expression matrix or a SummarizedExperiment object
#' (columns: samples same order as met, rows: genes)
#' @param metadata A data frame with samples as rownames and one columns that will be used to
#' color the samples
#' @param tf.activity.es A matrix with normalized enrichment scores for each TF across all samples
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
#' @param genome Genome of reference to be added to the plot as text
#' @param label.dnam Used for label text. Option "beta-value" and "residuals"
#' @param label.exp Used for label text. Option "expression" and "residuals"
#' @param add.tf.vs.exp.scatter.plot Add another row to the figure if the 
#' target gene expression vs TF expression stratified by DNA methylation groups
#' (DNAmLow - low quartile, DNAmHigh - high quartile)
#' @return A ggplot object, includes a table with results from fitting interaction model,
#' and the the following scatter plots: 1) TF vs DNAm, 2) Target vs DNAm,
#' 3) Target vs TF, 4) Target vs TF for samples in Q1 and Q4 for DNA methylation,
#' 5) Target vs DNAm for samples in Q1 and Q4 for the TF
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
#' plots <- plot_interaction_model(
#'     triplet.results = results,
#'     dnam = dnam,
#'     exp = exp
#' )
#' @export
#' @importFrom ggpubr ggplot ggarrange ggtexttable ttheme stat_cor
#' @importFrom ggplot2 xlab ylab geom_smooth
#' @importFrom tibble as_tibble
plot_interaction_model <-  function(
  triplet.results,
  dnam,
  exp,
  metadata,
  tf.activity.es = NULL,
  tf.dnam.classifier.pval.thld = 0.001,
  dnam.group.threshold = 0.25, 
  label.dnam = "beta-value",
  label.exp = "expression",
  genome = "hg38",
  add.tf.vs.exp.scatter.plot = FALSE
){
  
  if(!is(dnam.group.threshold,"numeric")) stop("dnam.group.threshold should be a value in the following interval (0,0.5]")
  if(dnam.group.threshold > 0.5) stop("dnam.group.threshold maximum valuee is 0.5")
  
  genome <- match.arg(genome, choices = c("hg38","hg19"))
  label.dnam <- match.arg(label.dnam, choices = c("beta-value","residuals"))
  label.exp <- match.arg(label.exp, choices = c("expression","residuals"))
  
  if (missing(dnam)) stop("Please set dnam argument with DNA methylation matrix")
  
  if (missing(exp)) stop("Please set exp argument with gene expression matrix")
  
  if (missing(triplet.results)) {
    stop("Please set triplet argument with interactors (region,TF, target gene) data frame")
  }
  
  if (!all(c("regionID","TF","target") %in% colnames(triplet.results))) {
    stop("triplet must have the following columns names: regionID, TF, target")
  }
  
  if (is(dnam,"SummarizedExperiment")) dnam <- assay(dnam)
  if (is(exp,"SummarizedExperiment")) exp <- assay(exp)
  
  
  triplet.results <- triplet.results %>% as.data.frame()
  
  out <- plyr::alply(
    .data = triplet.results,
    .margins = 1,
    .fun = function(row.triplet,metadata){
      
      
      df <- get_triplet_data(
        exp = exp,
        dnam = dnam,
        row.triplet = row.triplet,
        tf.es =  tf.activity.es
      )
      
      if(!any(grepl("DNAm.effect", colnames(row.triplet)))){
        stratified.results <- stratified_model_results(
          df,
          tf.dnam.classifier.pval.thld
        )
        if(!is.null(tf.activity.es)){
          colnames(stratified.results)[1:4] <- gsub("rna","es",colnames(stratified.results)[1:4])
        }
        
        row.triplet <- cbind(
          row.triplet,
          stratified.results
        )
      }
      
      color <- NULL
      if (!missing(metadata)) {
        df <- cbind(df,metadata)
        color <- colnames(metadata)[1]
      }
      
      # Reformat p-values for better looking on the plots
      idx <- grep("pval|fdr|value",colnames(row.triplet))
      idx <- idx[!is.na(as.numeric(row.triplet[idx]))]
      row.triplet[,idx] <- format.pval(row.triplet[,idx] %>% as.numeric(),digits = 3)
      idx <- grep("estimate|median|minus",colnames(row.triplet))
      idx <- idx[!is.na(as.numeric(row.triplet[idx]))]
      row.triplet[,idx] <- format(row.triplet[,idx] %>% as.numeric(),digits = 3)
      
      plots <- get_plot_results(
        df = df,
        row.triplet = row.triplet,
        color = color,
        use_tf_enrichment_scores = !is.null(tf.activity.es),
        label.dnam = label.dnam,
        label.exp = label.exp,
        dnam.group.threshold = dnam.group.threshold
      )
      
      table.plots <- get_table_plot(row.triplet, genome)
      
      # Arrange the plots on the same page
      suppressWarnings({
        suppressMessages({
          
          # Left side of the plot with all information from linear models
          plot.tables <- ggarrange(
            table.plots$table.plot.metadata,
            #table.plots$table.plot.wilcoxon,
            #table.plots$table.plot.lm.all,
            table.plots$table.plot.lm.quant,
            table.plots$table.plot.legend,
            #table.plots$table.plot.lm.dna.low,
            #table.plots$table.plot.lm.dna.high,
            heights = c(0.8,0.5,0.3),
            ncol = 1)
          
          # Right side of the plot with all scatter, box plots
          plots.top <- ggarrange(
            plots$tf.target,
            plots$dnam.target,
            #plots$dnam.tf,
            ncol = 2
          )
          
          if(add.tf.vs.exp.scatter.plot) {
            plots <- ggarrange( 
              plots.top,
              plots$tf.target.quantile,
              plots$dnam.target.quantile,
              nrow = 3
            )
          } else {
            plots <- ggarrange( 
              plots.top,
              plots$tf.target.quantile,
              nrow = 2
            )
          }
          
          # Tables on left, plots on right
          plot.all <- ggarrange(plot.tables, plots, ncol = 2,widths = c(1,2))
        })
      })
      plot.all + ggplot2::theme(plot.margin = grid::unit(c(1,2,1,2), "cm"))
    }, .progress = "time", metadata = metadata)
  attr(out,"split_type") <- NULL
  attr(out,"split_labels") <- NULL
  
  if (nrow(triplet.results) > 0) {
    names(out) <- paste0(
      triplet.results$regionID,
      "_TF_",
      triplet.results$TF,
      "_target_",
      triplet.results$target
    )
  }
  out
}

get_table_plot <- function(row.triplet, genome){
  
  row.triplet$genome <- genome
  columns <-  c(
    "genome",
    "regionID",
    "probeID",
    "target",
    "target_symbol",
    "target_region",
    "distance_region_target_tss",
    "TF",
    "TF_symbol",
    # "met.IQR",
    "TF.role",
    "DNAm.effect"
  )
  labels <-  c(
    "Genome of reference",
    "Region ID",
    "Probe ID",
    "Target gene ID",
    "Target gene Symbol",
    "Target region",
    "Distance Region to Target TSS",
    "TF gene ID",
    "TF gene Symbol",
    # "Diff. DNAm (Q4 - Q1)",
    "TF role",
    "DNAm effect on TF"
  )
  
  # Add probe ID if existant
  idx <- columns %in% colnames(row.triplet)
  columns <- columns[idx]
  labels <- labels[idx]
  
  base_size <- 9
  tab <- row.triplet %>%
    dplyr::select(
      columns
    ) %>% t() %>% as_tibble(rownames = "Variable")
  
  tab$Variable <- labels
  
  table.plot.metadata <- kable(tab)
  
  tab <- row.triplet %>%
    dplyr::select(
      c("TF_DNAm_high_vs_TF_DNAm_low_wilcoxon_pvalue")
    ) %>% t() %>% as_tibble(rownames = "Variable")
  
  tab$Variable <- c(
    "TF Q4 vs TF Q1"
  )
  #table.plot.wilcoxon <- ggtexttable(
  #  tab,
  #  rows = NULL,
  #  cols = c("Wilcoxon","P-Values"),
  #  theme = ttheme("mGreen", base_size = base_size)
  #)
  
  # Get results for linear model with all samples
  #table.plot.lm.all <- get_table_plot_results(row.triplet, type = "all")
  
  # Get results for linear model with all samples
  table.plot.lm.quant <- get_table_plot_results(row.triplet, type = "quantile")
  
  # Get results for linear model with DNAm low samples
  #table.plot.lm.dna.low <- get_table_plot_results(row.triplet, type = "DNAmlow")
  
  # Get results for linear model with DNAm high samples
  #table.plot.lm.dna.high <- get_table_plot_results(row.triplet, type = "DNAmhigh")
  
  table.plot.legend <- kable(data.frame(Legend = c("rlm: robust linear model", "ns: not significant")))
  
  table.plot.genome <- kable(data.frame("Genome of reference" = genome))
  
  table.plot.list <- list(
    "table.plot.metadata" = table.plot.metadata,
    #"table.plot.lm.all" = table.plot.lm.all,
    "table.plot.legend" = table.plot.legend,
    "table.plot.genome" = table.plot.genome,
    "table.plot.lm.quantile" = table.plot.lm.quant,
    #"table.plot.wilcoxon" = table.plot.wilcoxon,
    #"table.plot.lm.dna.low" = table.plot.lm.dna.low,
    #"table.plot.lm.dna.high" = table.plot.lm.dna.high
  )
  
  return(table.plot.list)
}

#' @importFrom stats quantile lm
#' @importFrom MASS rlm
get_plot_results <- function(
  df,
  row.triplet,
  color,
  use_tf_enrichment_scores = FALSE,
  label.dnam, 
  label.exp,
  dnam.group.threshold = 0.25
){
  
  target.lab <- bquote(atop("Target" ~.(row.triplet$target_symbol %>% as.character()), ~.(label.exp)))
  region.lab <- bquote(~.(row.triplet$regionID %>% as.character()) ~ "DNA methylation" ~.(label.dnam))
  if (use_tf_enrichment_scores) {
    tf.lab <- bquote("TF" ~.(row.triplet$TF_symbol %>% as.character()) ~" activity")
  } else {
    tf.lab <- bquote("TF" ~.(row.triplet$TF_symbol %>% as.character()) ~.(label.exp))
  }
  
  # quantile plots met
  quantile_upper_cutoff <-  quantile(df$met,na.rm = TRUE,  1 - dnam.group.threshold)
  quantile_lower_cutoff <-  quantile(df$met,na.rm = TRUE,  dnam.group.threshold)
  
  range1 <- paste0("[",paste(round(c(min(df$met,na.rm = TRUE),quantile_lower_cutoff), digits = 3),collapse = ","),"]")
  range2 <- paste0("[",paste(round(c(quantile_upper_cutoff,max(df$met,na.rm = TRUE)), digits = 3),collapse = ","),"]")
  
  df$DNAm.group <- NA
  df$DNAm.group[df$met >= quantile_upper_cutoff] <- paste0("DNAm high quartile ", range2)
  df$DNAm.group[df$met <= quantile_lower_cutoff] <- paste0("DNAm low quartile " , range1)
  
  df$DNAm.group <- factor(
    df$DNAm.group,
    levels = c(
      paste0("DNAm low quartile " , range1),
      paste0("DNAm high quartile ", range2)
    )
  )
  
  # quantile plots TF
  quantile_upper_cutoff <-  quantile(df$rna.tf,na.rm = TRUE,  1 - dnam.group.threshold)
  quantile_lower_cutoff <-  quantile(df$rna.tf,na.rm = TRUE,  dnam.group.threshold)
  
  range1 <- paste0("[",paste(round(c(min(df$rna.tf,na.rm = TRUE),quantile_lower_cutoff), digits = 3),collapse = ","),"]")
  range2 <- paste0("[",paste(round(c(quantile_upper_cutoff,max(df$rna.tf,na.rm = TRUE)), digits = 3),collapse = ","),"]")
  
  df$TF.group <- NA
  df$TF.group[df$rna.tf >= quantile_upper_cutoff] <- paste0("TF high quartile ", range2)
  df$TF.group[df$rna.tf <= quantile_lower_cutoff] <- paste0("TF low quartile ", range1)
  df$TF.group <- factor(
    df$TF.group,
    levels = c(
      paste0("TF low quartile " , range1),
      paste0("TF high quartile ", range2)
    )
  )
  
  
  tf.target.plot <- get_scatter_plot_results(
    df,
    x = "rna.tf",
    y = "rna.target",
    color = color,
    xlab = tf.lab,
    ylab = target.lab
  )
  
  #dnam.target.plot <- get_scatter_plot_results(
  #    df,
  #    x = "met",
  #    y = "rna.target",
  #    color = color,
  #    ylab = target.lab,
  #    xlab = region.lab
  #)
  
  dnam.target.plot <- get_box_plot_results(
    df,
    y = "rna.target",
    facet.by = "DNAm.group",
    ylab = target.lab
  )
  
  # dnam.target.plot <- get_histogram_plot_results(
  #     df,
  #    x = "rna.target",
  #    facet.by = "DNAm.group",
  #    xlab = target.lab
  # )
  
  dnam.tf.plot <- get_scatter_plot_results(
    df,
    x = "met",
    y = "rna.tf",
    color = color,
    ylab = tf.lab,
    xlab = region.lab
  )
  
  tf.target.quantile.plot <- get_scatter_plot_results(
    df[!is.na(df$DNAm.group),],
    x = "rna.tf",
    xlab = tf.lab,
    y = "rna.target",
    ylab = target.lab,
    facet.by = "DNAm.group",
    color = color # ifelse(is.null(color),"met",color)
  )
  
  dnam.target.quantile.plot <- get_scatter_plot_results(
    df[!is.na(df$TF.group),],
    x = "met",
    y = "rna.target",
    facet.by = "TF.group",
    color = color,
    ylab = target.lab,
    xlab = region.lab
  )
  
  return(
    list(
      "dnam.target.quantile" = dnam.target.quantile.plot,
      "tf.target.quantile" = tf.target.quantile.plot,
      "dnam.tf" = dnam.tf.plot,
      "dnam.target" = dnam.target.plot,
      "tf.target" = tf.target.plot
    )
  )
}

#' @noRd
#' @examples
#' df <- data.frame(
#'    x = runif(20),
#'    group = c(rep("low",10),rep("high",10))
#' )
#' get_histogram_plot_results(df = df, x =  "x",facet.by = "group", xlab = "expr")
get_histogram_plot_results <- function(
  df,
  x,
  facet.by,
  xlab
){
  
  df <- df[!is.na(df[[facet.by]]),]
  df[[facet.by]] <-  ifelse(grepl("low",df[[facet.by]]),"DNAm.low","DNAm.high")
  
  suppressWarnings({
    p <- ggplot(
      df,
      x = x,
      add = "mean",
      rug = TRUE,
      fill = facet.by,
      palette = c("#00AFBB", "#E7B800"),
      add_density = FALSE
    ) + geom_histogram() + labs(x = xlab) + ggplot2::theme(legend.title = ggplot2::element_blank())
  })
  p
}

#' @noRd
#' @examples
#' df <- data.frame(
#'    exp = runif(20),
#'    group = c(rep("high",10),rep("low",10))
#' )
#' get_box_plot_results(df = df, y =  "exp",facet.by = "group", ylab = "expr")
get_box_plot_results <- function(
  df,
  y,
  facet.by,
  ylab
){
  
  df <- df[!is.na(df[[facet.by]]),]
  df[[facet.by]] <- factor(
    ifelse(grepl("low",df[[facet.by]]),"DNAm.low","DNAm.high"),
    levels = c("DNAm.low","DNAm.high")
  )
  
  suppressWarnings({
    p <- ggplot(
      data = df,
      x = facet.by,
      y = y,
      add = "jitter",
      color = facet.by,
      palette = c("#477dbf", "#d04042"),
      add_density = FALSE
    ) + geom_boxplot() + labs(y = ylab) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      xlab("") +
      geom_signif(comparisons = list(c("group1", "group2")),
                  y_position = max(df[[y]]) * 1.1,
                  tip_length = 0) +
      ggplot2::guides( color = FALSE) +
      ggplot2::ylim(min(df[[y]]),max(df[[y]]) * 1.2)
  })
  p
}

#' @importFrom sfsmisc f.robftest
#' @importFrom stats as.formula
#' @noRd
#' @examples
#' df <- data.frame(
#' x = runif(20),
#' y = runif(20),
#' group = sample(c("high","low"),20,replace = TRUE)
#' )
#' get_scatter_plot_results(df, "x","y",NULL, xlab = "x", ylab = "y")
#' get_scatter_plot_results(df, "x","y",NULL, xlab = "x", ylab = "y", facet.by = "group")
#'
get_scatter_plot_results <- function(
  df,
  x,
  y,
  color,
  xlab,
  ylab,
  facet.by
){
  if (missing(facet.by)) {
    if (!is.null(color)) {
      p <- ggplot(
        df,
        x = x,
        y = y,
        color = color,
        size = 1
      )
    } else {
      p <- ggplot(
        df,
        x = x,
        y = y,
        size = 1
      )
    }
  } else {
    if (!is.null(color)) {
      p <- ggplot(
        df,
        x = x,
        y = y,
        facet.by = facet.by,
        color = color,
        size = 1
      ) #+ ggplot2::theme(legend.position = "right")
    } else {
      p <- ggplot(
        df,
        x = x,
        y = y,
        facet.by = facet.by,
        size = 1
      )
    }
  }
  
  p <- p + xlab(xlab) + ylab(ylab)
  
  suppressWarnings({
    suppressMessages({
      p <- p + geom_smooth(
        method = MASS::rlm,
        method.args = list(psi = "psi.bisquare"),
        se = FALSE
      )
    })
  })
  
  if (missing(facet.by)) {
    rlm.res <- get_rlm_val_pval(df, x, y)
    
    p <-  p + ggplot2::annotate(
      geom = "text",
      x = min(df[[x]], na.rm = TRUE),
      y = max(c(df[[y]] * 1.2, df[[y]] + 2), na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      color = 'blue',
      label = paste0(
        #gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
        "rlm estimate = ",
        formatC(
          rlm.res$rlm.val,
          digits = 3,
          format = ifelse(is.nan(rlm.res$rlm.val), "e",ifelse(abs(rlm.res$rlm.val) < 10^-3, "e","f"))
        ),
        "\nrlm p-value = ",
        formatC(
          rlm.res$rlm.p.value,
          digits = 3,
          format = ifelse(is.nan(rlm.res$rlm.p.value), "e",ifelse(rlm.res$rlm.p.value < 10^-3, "e","f"))
        )
      )
    )
  } else {
    # lower Annotation
    rlm.res.low <- get_rlm_val_pval(df %>% dplyr::filter(grepl("low",df[[facet.by]])),x , y)
    
    ann_text.low <- data.frame(
      x = min(df[[x]], na.rm = TRUE),
      y = max(c(df[[y]] * 1.2, df[[y]] + 2), na.rm = TRUE),
      facet.by = unique(factor(grep("low",df[[facet.by]],value = TRUE),levels = unique(df[[facet.by]])))
    )
    colnames(ann_text.low) <- c(x,y,facet.by)
    
    # higher Annotation
    p <- p + ggplot2::geom_text(
      #family = "Times New Roman",
      data = ann_text.low,
      hjust = 0,
      vjust = 1,
      color = 'blue',
      label = paste0(
        # gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
        "rlm estimate = ",
        formatC(
          rlm.res.low$rlm.val,
          digits = 3,
          format = ifelse(abs(rlm.res.low$rlm.val) < 10^-3, "e","f")
        ),
        "\nrlm p-value  = ",
        formatC(
          rlm.res.low$rlm.p.value,
          digits = 3,
          format = ifelse(is.nan(rlm.res.low$rlm.p.value), "e",ifelse(rlm.res.low$rlm.p.value < 10^-3, "e","f"))
        )
      )
    )
    
    rlm.res.high <- get_rlm_val_pval(df %>% dplyr::filter(grepl("high",df[[facet.by]])), x , y)
    ann_text.high <- data.frame(
      x = min(df[[x]], na.rm = TRUE),
      y = max(c(df[[y]] * 1.2, df[[y]] + 2), na.rm = TRUE),
      facet.by = unique(factor(grep("high",df[[facet.by]],value = TRUE),levels = unique(df[[facet.by]])))
    )
    colnames(ann_text.high) <- c(x,y,facet.by)
    
    # higher Annotation
    p <- p + ggplot2::geom_text(
      #family = "Times New Roman",
      data = ann_text.high,
      hjust = 0,
      vjust = 1,
      color = 'blue',
      label = paste0(
        #gsub("rna\\.","",y), " ~ ", gsub("rna\\.","",x),
        "rlm estimate = ",
        formatC(
          rlm.res.high$rlm.val,
          digits = 3,
          format =  ifelse(abs(rlm.res.high$rlm.val) < 10^-3, "e","f")
        ),
        "\nrlm p-value = ",
        formatC(
          rlm.res.high$rlm.p.value,
          digits = 3,
          format =  ifelse(is.nan(rlm.res.high$rlm.p.value), "e",ifelse(rlm.res.high$rlm.p.value < 10^-3, "e","f"))
        )
      )
    )
  }
  return(p)
  # stat_cor(method = "spearman",color = "blue")
}

get_rlm_val_pval <- function(df, x, y){
  suppressMessages({
    suppressWarnings({
      rls <- MASS::rlm(
        formula = as.formula(paste0(y, "~",x)),
        data = df,
        psi = psi.bisquare,
        maxit = 100
      )
      rlm <- rls %>% summary %>% coef %>% data.frame
      rlm.val <- rlm[-1,1]
    })
  })
  
  rlm.p.value <- tryCatch({
    degrees.freedom.value <- nrow(df) - 2
    pval <- 2 * (1 - pt( abs(rlm$t.value[-1]), df = degrees.freedom.value))
    pval
    #ftest <- sfsmisc::f.robftest(rls)
    #ftest$p.value
  }, error = function(e){
    #message(e);
    return("NA")
  })
  
  return(
    list(
      rlm.p.value = rlm.p.value,
      rlm.val = rlm.val
    )
  )
}

get_table_plot_results <- function(row.triplet, type){
  
  base.size <- 9
  
  if (type == "all") {
    pattern.estimate <- "^estimate"
    pattern.pval <- "^pval"
    title <- "Target ~ TF + DNAm +\n TF * DNAm"
    theme.color <- "mOrange"
  } else if (type == "quantile") {
    pattern.estimate <- "^RLM_.*estimate$"
    pattern.pval <- "^RLM_.*pvalue$"
    title <- "Target ~ TF + \nDNAm Quant. Group +\n TF * DNAm Quant. Group"
    theme.color <- "mGreen"
  } else if (type == "DNAmlow") {
    pattern.estimate <- "^DNAm_low_.*estimate"
    pattern.pval <- "^DNAm_low_.*pvalue"
    title <- "Target ~ TF\nDNAm low samples"
    theme.color <- "mBlue"
  } else if (type == "DNAmhigh") {
    pattern.estimate <-"^DNAm_high_.*estimate"
    pattern.pval <- "^DNAm_high_.*pvalue"
    title <- "Target ~ TF\nDNAm high samples"
    theme.color <- "mRed"
  }
  
  col.idx <- grep(pattern.estimate,colnames(row.triplet),value = TRUE)
  table.plot.estimate <- row.triplet[,col.idx,drop  = FALSE] %>%
    t() %>%
    as_tibble(rownames = "Variable")
  colnames(table.plot.estimate)[2] <- "Estimate"
  table.plot.estimate$Variable <- gsub(
    paste0("_estimate"),"", table.plot.estimate$Variable
  )
  
  table.plot.estimate$Variable[grep("DNAmGroup:TF",table.plot.estimate$Variable)] <- "Synergistic effect\n of DNAm and TF"
  table.plot.estimate$Variable[grep("DNAmGroup",table.plot.estimate$Variable)] <- "Direct effect of DNAm"
  table.plot.estimate$Variable[grep("_TF",table.plot.estimate$Variable)] <- "Direct effect of TF"
  
  
  col.idx <- gsub("estimate","pvalue",col.idx)
  table.plot.pval <- row.triplet[,col.idx, drop  = FALSE] %>%
    t() %>%
    as_tibble(rownames = "Variable")
  
  table.plot.pval$Variable <- gsub(
    pattern = "_pvalue",
    replacement = "",
    table.plot.pval$Variable
  )
  
  table.plot.pval$Variable[grep("DNAmGroup:TF",table.plot.pval$Variable)] <- "Synergistic effect\n of DNAm and TF"
  table.plot.pval$Variable[grep("DNAmGroup",table.plot.pval$Variable)] <- "Direct effect of DNAm"
  table.plot.pval$Variable[grep("_TF",table.plot.pval$Variable)] <- "Direct effect of TF"
  
  colnames(table.plot.pval)[2] <- "P-value"
  
  table.plot <- merge(
    table.plot.estimate,
    table.plot.pval,
    by = "Variable",
    sort = FALSE
  )
  
  table.plot.lm.all <- kable(data.frame(
    title = title,
    Estimate = table.plot$Estimate,
    `P-Values` = table.plot$`P-value`
  ))
  
  table.plot.lm.all
}
