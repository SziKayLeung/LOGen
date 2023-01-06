# volcano plots of differentially expressed genes and transcripts
plot_volcano <- function(diff_results){
  
  #https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/#:~:text=A%20volcano%20plot%20is%20a,tools%20like%20EdgeR%20or%20DESeq2.
  diff_results <- diff_results %>% mutate(
    Expression = case_when(log2FC >= log(2) & `1-prob` <= 0.05 ~ "Up-regulated",
                           log2FC <= -log(2) & `1-prob` <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
  
  cat("Number of transcripts upregulated (red):", nrow(diff_results[diff_results$Expression == "Up-regulated",]))
  cat("Number of transcripts downregulated (blue):", nrow(diff_results[diff_results$Expression == "Down-regulated",]))
  
  top <- 10
  top_genes <- bind_rows(
    diff_results %>% 
      filter(Expression == 'Up-regulated') %>% 
      arrange(`1-prob`, desc(abs(log2FC))) %>% 
      head(top),
    diff_results %>% 
      filter(Expression == 'Down-regulated') %>% 
      arrange(`1-prob`, desc(abs(log2FC))) %>% 
      head(top)
  )
  
  p <- ggplot(diff_results, aes(log2FC, -log(1-prob,10))) + # -log10 conversion  
    geom_point(aes(color = Expression), size = 2/5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    ggrepel::geom_label_repel(data = top_genes,
                              mapping = aes(log2FC, -log(1-prob,10), label = associated_gene),
                              size = 2)
  
  output <- list(p, top_genes)
  names(output) <- c("p","top10")
  return(output)
  
}