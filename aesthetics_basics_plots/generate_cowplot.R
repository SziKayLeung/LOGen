# ease cowplot functions
generate_cowplot <- function(...,num,ncol,nrow){
  if(num == "2"){label_input = c("A","B")}
  else if(num == "1"){label_input = c("")}
  else if(num == "4"){label_input = c("A","B","C","D")}
  else if(num == "1A"){label_input = c("A")}
  else if(num == "2BC"){label_input = c("B","C")}
  else if(num == "4D"){label_input = c("D")}
  else if(num == "3"){label_input = c("A","B","C")}
  else if(num == "3rw"){label_input = c("A","B","C")}
  else if(num == "6"){label_input = c("A","B","C","D","E","F")}
  else if(num == "7"){label_input = c("A","B","C","D","E","F","G")}
  else if(num == "8"){label_input = c("A","B","C","D","E","F","G","H")}
  else if(num == "3,-B"){label_input = c("A","","B","C")}
  else if(num == "5,-B"){label_input = c("A","","C","D","E")}
  else (print("num required"))
  
  p = plot_grid(...,labels = label_input, label_size = 30, label_fontfamily = "CM Roman", ncol = ncol, nrow = nrow, scale = 0.9)
  
  if(num == "rw"){p = plot_grid(...,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(2,1))}
  if(num == "3rw"){p = plot_grid(...,label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_heights = c(0.6,0.4,0.4), rel_widths = c(0.6,0.4))}
  return(p)
}


generate_cowplot <- function(...,num,ncol,nrow){
  if(num == "2"){label_input = c("A","B")}
  else if(num == "A"){label_input = c("A")}
  else if(num == "3"){label_input = c("A","B","C")}
  else if(num == "4"){label_input = c("A","B","C","D")}
  else if(num == "rw2"){label_input = c("B","C")}
  else if(num == "rw3"){label_input = c("B","C","D")}
  else if(num == "rw2D"){label_input = c("E","F","G","H")}
  else if(num == "rw4"){label_input = c("B","C","D","E")}
  else if(num == "rw8"){label_input = c("F","G","H","I")}
  else if(num == "rw"){label_input = c("A","B")}
  else if(num == "2DE"){label_input = c("D","E")}
  else if(num == "BCDE"){label_input = c("B","C","D","E")}
  else if(num == "BC"){label_input = c("B","C")}
  else if(num == "DE"){label_input = c("D","E")}
  else if(num == "BCD"){label_input = c("B","C","D")}
  else if(num == "rwBC"){label_input = c("B","C")}
  else if(num == "6"){label_input = c("A","B","C","D","E","F")}
  else if(num == "6B"){label_input = c("B","C","D","E","F","G")}
  else if(num == "7"){label_input = c("A","B","C","D","E","F","G")}
  else if(num == "5,-B"){label_input = c("A","","C","D","E")}
  else (print("num required"))
  
  p = plot_grid(...,labels = label_input, label_size = 30, label_fontfamily = "CM Roman", ncol = ncol, nrow = nrow, scale = 0.9)
  
  if(num == "rw"){p = plot_grid(...,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(2,1))}
  
  if(num == "rw2"){p = plot_grid(...,labels = c("B","C"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.8))}
  
  if(num == "rw3"){p = plot_grid(...,labels = c("B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.5))}
  
  if(num == "rw2D"){p = plot_grid(...,labels = c("D","E","F","G"), label_size = 30, label_fontfamily = "CM Roman", 
                                  ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rw4"){p = plot_grid(...,labels = c("A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rw8"){p = plot_grid(...,labels = c("E","F","G","H"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.2,0.3,0.2,0.3))}
  
  if(num == "rwBC"){p = plot_grid(...,labels = c("B","C"), label_size = 30, label_fontfamily = "CM Roman", 
                                  ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.6,0.4))}
  
  if(num == "BCD"){p = plot_grid(...,labels = c("B","C","D"), label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(0.55,0.225,0.225))}
  return(p)
}

generate_cowplot <- function(...,num,ncol,nrow){
  if(num == "2"){label_input = c("A","B")}
  else if(num == "1"){label_input = c("")}
  else if(num == "4"){label_input = c("A","B","C","D")}
  else if(num == "1A"){label_input = c("A")}
  else if(num == "2BC"){label_input = c("B","C")}
  else if(num == "4D"){label_input = c("D")}
  else if(num == "3"){label_input = c("A","B","C")}
  else if(num == "3rw"){label_input = c("A","B","C")}
  else if(num == "6"){label_input = c("A","B","C","D","E","F")}
  else if(num == "7"){label_input = c("A","B","C","D","E","F","G")}
  else if(num == "8"){label_input = c("A","B","C","D","E","F","G","H")}
  else if(num == "3,-B"){label_input = c("A","","B","C")}
  else if(num == "5,-B"){label_input = c("A","","C","D","E")}
  else (print("num required"))
  
  p = plot_grid(...,labels = label_input, label_size = 30, label_fontfamily = "CM Roman", ncol = ncol, nrow = nrow, scale = 0.9)
  
  if(num == "rw"){p = plot_grid(...,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", 
                                ncol = ncol, nrow = nrow, scale = 0.9,rel_widths = c(2,1))}
  if(num == "3rw"){p = plot_grid(...,label_size = 30, label_fontfamily = "CM Roman", 
                                 ncol = ncol, nrow = nrow, scale = 0.9,rel_heights = c(0.6,0.4,0.4), rel_widths = c(0.6,0.4))}
  return(p)
}

