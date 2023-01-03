tabulate_structural_cate <- function(class.files){
  
  # group by structural_category and tally with proportion
  output <- class.files %>% group_by(structural_category) %>% tally() %>% mutate(perc = n/sum(n)*100)
  
  return(output)
}