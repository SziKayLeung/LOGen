
# number of transcripts that are filtered due to low expression as part of tappAS
num_tappas_filter <- function(normalised_counts, targeted.class.files){
  normalised_counts <- normalised_counts %>% rownames_to_column(., var = "isoform")
  retained_transcript <- intersect(targeted.class.files$isoform, normalised_counts$isoform)
  
  print(paste0("Number of targeted transcripts from SQANTI3 annotation: ", 
               length(targeted.class.files$isoform)))
  print(paste0("Number of retained targeted transcripts kept after filtering for low expression: ", 
               length(retained_transcript)))
  
  p1 <- targeted.class.files %>% mutate(Filter = ifelse(isoform %in% normalised_counts$isoform, "Retained","Filtered")) %>% 
    group_by(associated_gene, Filter) %>% tally() %>% 
    ggplot(., aes(x = associated_gene, y = n, fill= Filter)) + geom_bar(stat = "identity") + mytheme + 
    labs(x = "", y = "Number of Isoforms") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
  
  p2 <- targeted.class.files %>% mutate(Filter = ifelse(isoform %in% normalised_counts$isoform, "Retained","Filtered")) %>% 
    filter(Filter == "Filtered") %>%
    group_by(associated_gene, structural_category) %>% tally() %>%
    ggplot(., aes(x = associated_gene, y = n, fill = structural_category)) + geom_bar(stat = "identity") + mytheme + 
    labs(x = "", y = "Number of Isoforms") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") 
  
  return(list(p1,p2))
}

# manual means of filtering isoforms from the original input normalised dataframe
manual_filter_isoforms <- function(original_norm, norm_transcounts, gene, phenotype){
  dat <- norm_transcounts %>% filter(associated_gene == gene)
  
  control <- phenotype[phenotype$group == "CONTROL","sample"] 
  case <- phenotype[phenotype$group == "CASE","sample"] 
  control_exp <- original_norm[rownames(original_norm) %in% unique(dat$isoform),control]
  case_exp <- original_norm[rownames(original_norm) %in% unique(dat$isoform),case]
  
  exp_stats <- cbind(melt(rowSums(control_exp != 0)),
                     melt(apply(control_exp,1,mean)),
                     melt(apply(control_exp,1,sd)),
                     melt(rowSums(case_exp != 0)),
                     melt(apply(case_exp,1,mean)),
                     melt(apply(control_exp,1,sd))) %>% 
    `colnames<-`(c("num_no_0_ctrl", "ctrl_mean", "ctrl_sd","num_no_0_case","case_mean","case_sd")) %>% 
    mutate(ctrl_cov = ctrl_sd/ctrl_mean * 100, case_cov = case_sd/case_mean * 100) 
  
  exp_stats_filtered <- exp_stats %>% 
    filter(num_no_0_case != 0.00 & num_no_0_ctrl != 0.00) %>% 
    filter(ctrl_mean > 10 & case_mean > 10)
  
  output <- list(exp_stats, exp_stats_filtered)
  names(output) <- c("exp","exp_filtered")
  return(output)
}
