generate_plots <- function(genelist,tappasinput,type,name,...){
  parms=list(...)
  
  if(type == "Gene"){
    plist <- lapply(lapply(genelist, function(gene) plot_gene_exp(gene, tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],"case_control",name)),ggplotGrob)
  }else if(type == "Gene_time"){
    plist <- lapply(lapply(genelist, function(gene) plot_gene_exp(gene, tappasinput[["GeneExp"]],tappasinput[["Norm_transcounts"]],"time_series",name)),ggplotGrob)
  }else if(type == "Transcript_Iso_Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],name,"time_series_sig","isoseq",...)),ggplotGrob)
  }else if(type == "Transcript_Rna_Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],name,"time_series_sig","rnaseq",...)),ggplotGrob)
  }else if(type == "Transcript Trajectory"){
    plist <- lapply(lapply(genelist, function(gene) plot_general_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],...,name)),ggplotGrob)
  }else if(type == "Transcript"){
    plist <- lapply(lapply(genelist, function(gene) plot_trans_exp(gene,tappasinput[["Norm_transcounts"]],"all",name)),ggplotGrob)
  }else if(type == "Top10_Transcript"){
    plist <- lapply(lapply(genelist, function(gene) plot_trans_exp(gene,tappasinput[["Norm_transcounts"]],"top10",name)),ggplotGrob)
  }else if(type == "2Cate_Transcript"){
    plist <- lapply(lapply(genelist, function(trans) twocate_plot_transexp(trans,tappasinput,name)),ggplotGrob)
  }else if(type == "2Cate_Transcript_overtime"){
    plist <- lapply(lapply(genelist, function(gene) twocate_plot_transexp_overtime(gene,tappasinput[["Norm_transcounts"]],name)),ggplotGrob)
  }else if (type == "Per_Transcript"){
    plist <- lapply(lapply(genelist, function(gene) plot_trans_exp_individual(gene,tappasinput[["Norm_transcounts"]])), ggplotGrob)
  }else if (type == "usage"){
    # <gene> <loaded$gene_transcripts> <normalised_matrix> <phenotype> <class.file>
    plist <- lapply(lapply(genelist, function(gene) IF_plot(as.character(gene),
                                                            tappasinput[["gene_transcripts"]], tappasinput[["input_normalized_matrix"]], 
                                                            parms[[1]], parms[[2]])), ggplotGrob)  
  }else if(type == "IF_time_series"){
    plist <- lapply(lapply(genelist, function(gene) 
      IF_plot_time_series(gene,tappasinput[["gene_transcripts"]],tappasinput[["input_normalized_matrix"]],name,...)[[3]]),ggplotGrob)
    
  }else{
    print("Type Required")
  }
  
  names(plist) = genelist
  print(genelist)
  return(plist)
}
