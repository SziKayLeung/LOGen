## ---------- Manual Stats -----------------
# perform t-test/wilcox.test per gene between case and controls
gene_stats <- function(gene_matrix, phenotype_file){
  
  df <- loaded$iso$gene_matrix %>% rownames_to_column("gene") %>%
    melt(., id.vars = c("gene")) %>% dplyr::rename(sample = variable) %>%
    merge(., phenotype$iso, by = "sample")
  
  ttest_gene <- function(g){
    cat("******************", g,"***********\n")
    dat <- df %>% filter(gene == g) 
    t1 <- with(dat, shapiro.test(value[group == "CASE"]))
    t2 <- with(dat, shapiro.test(value[group == "CONTROL"]))
    t3 <- var.test(value ~ group, data = dat)
    
    if(t1$p.value > 0.05 & t2$p.value > 0.05 & t3$p.value > 0.05){
      print("Can assume normality, assumption met")
      res <- t.test(value ~ group, data = dat, var.equal = TRUE)
    }else{
      print("Cannot assume normality")
      res <- wilcox.test(value ~ group, data = dat, exact = FALSE)
    }
    print(res)
    return(res)
  }
  
  genetests <- lapply(unique(df$gene), function(x) ttest_gene(x))
  names(genetests) <- unique(df$gene)
  return(genetests)
}

## performing the t-test
# normalised_counts: normalised matrix for the input of expression values
# transcript: transcript of interest used to subset the normalised matrix
trans_stats <- function(normalised_counts, transcript){
  
  df <- normalised_counts[normalised_counts$isoform == transcript,] 
  
  # only perform the t-test if the transcript of interest is in the normalised matrix
  if(nrow(df) > 0){
    res <- t.test(value ~ group, data = df, var.equal = TRUE)
    # capture the transcript and p-value if significant
    if(res$p.value < 0.05){
      result <- list(transcript, res$p.value)
      return(result)
    }
  }else{
    exit()
  }
}

# running t-test to identify differentially expressed transcripts of target genes
run_trans_stats <- function(normalised_counts){
  
  normalised_counts <- normalised_counts %>% filter(associated_gene %in% TargetGene)
  
  # aggregate the results capturing the significant transcripts after performing t-test across all transcripts
  tdiff_stats <- lapply(unique(normalised_counts$isoform), function(x) trans_stats(normalised_counts, x)) %>% 
    .[lapply(.,length)>0] %>% do.call(rbind, .) %>% as.data.frame() %>% `colnames<-`(c("isoform", "P"))
  
  print(paste0("Number of differentially expressed transcripts:",nrow(tdiff_stats)))
  
  # plot the expression of significant transcripts
  tdiff_plots <- lapply(unlist(tdiff_stats$isoform), function(x) plot_trans_exp_individual(x, normalised_counts))
  
  # output
  output <- list(tdiff_stats, tdiff_plots)
  names(output) <- c("stats","plots")
  return(output)
}

